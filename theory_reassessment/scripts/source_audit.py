"""Read-only provenance of external theoretical foundations; run locally."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/zjsun/Desktop/Agent/try_sth/GTEMPO_theory')
W=Path('/home/zjsun/.codex/worktrees/4288/GTEMPO_theory/order2')
A=Path('/home/zjsun/.codex/worktrees/a0d0/GTEMPO_theory/order2/weyl_tti/if_constant_audit')
sources=[BASE/'order1_实时GTEMPO_从Hamiltonian到代码目标的完整推导.md',
 BASE/'GTEMPO_order1_虚时间逐步推导.md',BASE/'GTEMPO_理论与离散阶数审计.md',
 W/'weyl_tti/推导.md',W/'weyl_tti/验证报告.md',W/'weyl_tti/TTI审计.md',
 A/'理论审查.md',W/'core.py',W/'weyl_tti/weyl.py']
records=[]
for p in sources:
    records.append(dict(path=str(p),exists=p.exists(),sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None))
for name,p in [('core.py',W/'core.py'),('weyl.py',W/'weyl_tti/weyl.py')]:
    assert (ROOT/'reference'/name).read_bytes()==p.read_bytes()
out=dict(sources=records,read_only_external_sources=True,
 historical_threads=['01a159e9-e6c7-7300-932b-d89936c2ba49','01a15f48-b981-7a10-985f-7b95b5dada96','01a1644b-70cf-7340-aa9c-bd5b98dd9e43'],
 historical_session_lookup='SQLite state_5.sqlite opened with mode=ro; final_answer messages in rollout_path JSONL',
 note='The older general order audit contains a retraction; use the later complete first-order derivation.',
 versions=(ROOT.parent/'results/source_versions.txt').read_text(),
 upstream_HEAD_rechecked='2026-11-03; GTEMPO and QuAPI git ls-remote HEAD matched the recorded commits')
(ROOT/'results/source_audit.json').write_text(json.dumps(out,indent=2,ensure_ascii=False))
print('External source hashes recorded; copied reference utilities match exactly.')
