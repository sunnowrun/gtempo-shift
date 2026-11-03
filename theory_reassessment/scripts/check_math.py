"""Explicitly check math delimiters in authored Markdown, not vendor docs."""
from pathlib import Path
root=Path(__file__).resolve().parents[1]
paths=[root/'REPORT.md',root/'README.md',*sorted((root/'docs').glob('*.md'))]
for p in paths:
    if not p.exists():
        continue
    state=False
    for n,line in enumerate(p.read_text().splitlines(),1):
        i=0
        while i<len(line):
            if line.startswith('$$',i):
                state=not state
                i+=2
                continue
            if line[i]=='$':
                state=not state
            if line[i]=='\\' and not state:
                raise SystemExit(f'{p}:{n}: command outside math delimiters')
            i+=1
    if state:
        raise SystemExit(f'{p}: unclosed math delimiter')
    print(f'{p.relative_to(root)}: math delimiters explicitly checked')
