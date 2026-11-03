"""Verify recorded study artifacts and portable foundation snapshots offline."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    study = ROOT/'theory_reassessment'
    recorded = json.loads((study/'results/artifact_manifest.json').read_text())
    for name, expected in recorded.items():
        assert sha(study/name) == expected, name
    sources = json.loads((ROOT/'foundations/sources.json').read_text())
    for source in sources:
        assert sha(ROOT/source['snapshot']) == source['sha256'], source['snapshot']
    conversation = [json.loads(line) for line in (ROOT/'history/conversation.jsonl').read_text().splitlines()]
    assert all(row['role'] == 'user' or (row['role'] == 'assistant' and row['phase'] in ('commentary','final_answer','final')) for row in conversation)
    meta = json.loads((ROOT/'history/export_manifest.json').read_text())
    assert len(conversation) == meta['message_count']
    print(f'Verified {len(recorded)} study files, {len(sources)} foundation snapshots and {len(conversation)} conversation messages.')


if __name__ == '__main__':
    main()
