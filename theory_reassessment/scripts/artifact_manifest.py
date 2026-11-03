"""Checksum the deliverable; no numerical simulation or external writes."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
target=ROOT/'results/artifact_manifest.json'
records={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
         for p in sorted(ROOT.rglob('*')) if p.is_file() and p!=target}
target.write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n')
print('Checksummed',len(records),'deliverable files')
