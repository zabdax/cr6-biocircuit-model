"""Quick structural validation of the manuscript tex file.

Resolves the manuscript path relative to this script so the check works
from any working directory (and on any platform).
"""
import os
import re
from collections import Counter

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX_FILE = os.path.join(REPO_ROOT, "Journal_Manuscript_2026.tex")

with open(TEX_FILE, encoding='utf-8') as f:
    text = f.read()

# Brace balance
opens = text.count('{')
closes = text.count('}')
print(f'Braces: open={opens} close={closes} balance={opens - closes}')

# Environment balance
begins = re.findall(r'\\begin\{(\w+)\*?\}', text)
ends = re.findall(r'\\end\{(\w+)\*?\}', text)
b_counts = Counter(begins)
e_counts = Counter(ends)
all_envs = set(list(b_counts) + list(e_counts))
print(f'Begin envs: {len(begins)}; End envs: {len(ends)}')
mismatch = False
for env in sorted(all_envs):
    if b_counts[env] != e_counts[env]:
        print(f'  MISMATCH: {env}: {b_counts[env]} begin vs {e_counts[env]} end')
        mismatch = True
if not mismatch:
    print('  All environments balanced.')

# Sections
sections = re.findall(r'\\section\{([^}]+)\}', text)
print()
print('Top-level sections:')
for s in sections:
    print(f'  - {s}')

# Subsections
subs = re.findall(r'\\subsection\{([^}]+)\}', text)
print()
print('Subsections:')
for s in subs:
    print(f'  - {s}')

# Citations
cites = re.findall(r'\\cite\{([^}]+)\}', text)
all_keys = set()
for c in cites:
    for k in c.split(','):
        all_keys.add(k.strip())
print()
print(f'Citation references: {len(cites)}')
print(f'Unique citation keys used: {len(all_keys)}')

# Bibliography
bibs = re.findall(r'\\bibitem\{([^}]+)\}', text)
print(f'Bibliography entries: {len(bibs)}; unique: {len(set(bibs))}')

missing = all_keys - set(bibs)
extra = set(bibs) - all_keys
if missing:
    print(f'Cited but not in bib: {sorted(missing)}')
if extra:
    print(f'In bib but not cited: {sorted(extra)}')
if not missing and not extra:
    print('All citations match bibliography.')
