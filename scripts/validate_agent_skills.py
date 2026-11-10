#!/usr/bin/env python3
"""Lightweight stdlib structural check of the local Agent Skills package.
Not a replacement for the agentskills.io reference validator.
"""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent.parent
SKILLS = ROOT / '.agents' / 'skills'
errors = []
count = 0
for directory in sorted(SKILLS.iterdir() if SKILLS.is_dir() else []):
    if not directory.is_dir():
        continue
    skill = directory / 'SKILL.md'
    count += 1
    if not skill.exists():
        errors.append(f'{directory}: missing SKILL.md')
        continue
    data = skill.read_text(encoding='utf-8')
    match = re.match(r'\A---\n(.*?)\n---\n', data, flags=re.S)
    if not match:
        errors.append(f'{skill}: missing YAML frontmatter')
        continue
    fm = match.group(1)
    fields = {}
    for row in fm.splitlines():
        if ':' in row:
            name, value = row.split(':',1)
            fields[name.strip()] = value.strip()
    name = fields.get('name','').strip('"\'')
    desc = fields.get('description','').strip('"\'')
    if name != directory.name or not re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', name) or len(name)>64:
        errors.append(f'{skill}: invalid or mismatched name')
    if not desc or len(desc)>1024:
        errors.append(f'{skill}: missing/oversize description')
    for target in re.findall(r'\]\((references/[^)#]+)\)', data):
        if not (directory / target).is_file():
            errors.append(f'{skill}: missing linked reference {target}')
    print(f'checked: {directory.name}')

if count == 0:
    errors.append(f'No skills found in {SKILLS}')
if errors:
    print('\n'.join('ERROR: '+e for e in errors), file=sys.stderr)
    sys.exit(1)
print(f'OK: {count} skill(s) with metadata and references checked.')
