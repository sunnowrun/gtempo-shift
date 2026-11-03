"""Publish a redacted user-visible transcript from one explicitly selected session.

Reads the source session without modifying it. Internal instructions, reasoning,
compaction state, and telemetry are excluded. Tools that inspect session internals
retain their timestamp/name but not their arguments or outputs.
"""
import argparse
import collections
import datetime
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

SECRET_PATTERNS = [
    (re.compile(r"-----BEGIN (?:[A-Z0-9 ]+ )?PRIVATE KEY-----[\s\S]*?-----END (?:[A-Z0-9 ]+ )?PRIVATE KEY-----"), '[REDACTED PRIVATE KEY]'),
    (re.compile(r'\bgh[pousr]_[A-Za-z0-9_]{20,}\b|\bgithub_pat_[A-Za-z0-9_]{20,}\b'), '[REDACTED GITHUB TOKEN]'),
    (re.compile(r'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{20,}\b'), '[REDACTED API KEY]'),
    (re.compile(r'(?i)(authorization\s*[:=]\s*["\x27]?bearer\s+)[A-Za-z0-9._~+/-]{12,}'), r'\1[REDACTED]'),
    (re.compile(r'(?i)(https?://)[^\s/:@]+:[^\s/@]+@'), r'\1[REDACTED CREDENTIALS]@'),
    (re.compile(r'(?i)(["\x27]?(?:api[_-]?key|access[_-]?token|refresh[_-]?token|password|passwd|client[_-]?secret)["\x27]?\s*[:=]\s*)["\x27]([^"\x27\n]+)["\x27]'), r'\1"[REDACTED]"'),
]
PRIVATE_TOOL = re.compile(r'rollout-|state_5\.sqlite|\.codex/(?:sessions|archived_sessions)|read_thread|ALL_TOOLS')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rollout', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT/'history')
    args = parser.parse_args()
    destination = args.output.resolve()
    if not destination.is_relative_to(ROOT):
        raise SystemExit('Output must stay inside this project.')
    destination.mkdir(parents=True, exist_ok=True)
    raw = args.rollout.read_bytes()
    replacements = collections.Counter()

    def redact(value):
        for number, (pattern, replacement) in enumerate(SECRET_PATTERNS):
            value, count = pattern.subn(replacement, value)
            replacements[str(number)] += count
        return value

    conversation, tools = [], []
    private_calls = set()
    excluded = collections.Counter()
    session_id, cutoff = None, None
    for line in raw.decode('utf-8').splitlines():
        record = json.loads(line)
        payload = record.get('payload', {})
        cutoff = record.get('timestamp', cutoff)
        if record.get('type') == 'session_meta':
            session_id = payload.get('id')
        if record.get('type') != 'response_item':
            continue
        kind = payload.get('type')
        if kind == 'message':
            role = payload.get('role')
            phase = payload.get('phase', payload.get('channel'))
            if role != 'user' and not (role == 'assistant' and phase in ('commentary', 'final_answer', 'final')):
                excluded['internal_message'] += 1
                continue
            parts = payload.get('content', [])
            body = '\n'.join(p.get('text', '') for p in parts if isinstance(p, dict) and p.get('text'))
            conversation.append(dict(timestamp=record.get('timestamp'),role=role,phase=phase,text=redact(body)))
        elif kind in ('custom_tool_call', 'function_call'):
            call_id = payload.get('call_id', payload.get('id'))
            body = payload.get('input', payload.get('arguments', ''))
            if not isinstance(body,str):body=json.dumps(body,ensure_ascii=False)
            private = bool(PRIVATE_TOOL.search(body))
            if private:private_calls.add(call_id)
            tools.append(dict(timestamp=record.get('timestamp'),type='tool_call',name=payload.get('name'),call_id=call_id,
                              input='[Session-internal inspection omitted]' if private else redact(body)))
        elif kind in ('custom_tool_call_output', 'function_call_output'):
            call_id = payload.get('call_id')
            body = payload.get('output','')
            if not isinstance(body,str):body=json.dumps(body,ensure_ascii=False)
            private = call_id in private_calls or bool(re.search(r'"type"\s*:\s*"(?:reasoning|compacted|retained_context)"',body))
            tools.append(dict(timestamp=record.get('timestamp'),type='tool_result',call_id=call_id,
                              output='[Session-internal inspection omitted]' if private else redact(body)))
            if private:excluded['session_internal_tool_output']+=1
        else:
            excluded[kind or 'unknown'] += 1

    for filename, rows in [('conversation.jsonl',conversation),('tools.jsonl',tools)]:
        (destination/filename).write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in rows))
    (destination/'conversation.txt').write_text('\n\n'.join(
        f'[{r["timestamp"]}] {r["role"]} ({r["phase"] or "message"})\n{r["text"]}' for r in conversation)+'\n')
    manifest=dict(session_id=session_id,source_rollout_sha256=hashlib.sha256(raw).hexdigest(),
                  snapshot_cutoff_utc=cutoff,exported_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  message_count=len(conversation),tool_record_count=len(tools),
                  excluded_response_item_counts=dict(excluded),credential_redactions_by_rule=dict(replacements),
                  scope='Current project session only; user-visible messages and redacted tool records; no internal reasoning or instructions.')
    (destination/'export_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('session_id','message_count','tool_record_count','snapshot_cutoff_utc')},ensure_ascii=False))


if __name__ == '__main__':
    main()
