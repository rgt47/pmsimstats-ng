#!/usr/bin/env python3
"""Corrections to the compendium bibliographies, verified 2026-10-07
against PubMed. Entries are matched by title, so every copy is fixed
whatever its citation key.

Usage: python3 tools/fix-bib-errors-2026-10-07.py BACKUP_DIR FILE.bib [...]
Each changed file is first copied to BACKUP_DIR, then rewritten via a
temporary file in BACKUP_DIR and moved into place (no in-place editing
of cloud-synced paths).
"""
import os
import re
import shutil
import sys

DOIS = {
    'predictive approaches to treatment effect heterogeneity': '10.7326/M18-3667',
    'trial of prazosin for post-traumatic stress disorder in military veterans':
        '10.1056/NEJMoa1507598',
    'understanding variation in sets of n-of-1 trials': '10.1371/journal.pone.0167167',
    'combining single patient': '10.1016/S0895-4356(96)00429-5',
    'power and design issues in crossover-based': '10.3390/healthcare7030084',
    'overt versus covert treatment': '10.1016/S1474-4422(04)00908-1',
}


def norm(s):
    return re.sub(r'\s+', ' ', re.sub(r'[{}]', '', s)).strip().lower()


def entries(text):
    """Yield (start, end) spans of entries: from '@type{' to the next
    line that starts an entry or the end of the file."""
    starts = [m.start() for m in re.finditer(r'(?m)^@\w+\s*\{', text)]
    for i, s in enumerate(starts):
        yield s, starts[i + 1] if i + 1 < len(starts) else len(text)


def title_of(entry):
    m = re.search(r'(?is)\btitle\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|"[^"]*")', entry)
    return norm(m.group(1)) if m else ''


def fix_entry(entry, log):
    t = title_of(entry)
    new = entry
    new = re.sub(r'Zucker, David R\b\.?', 'Zucker, Deborah R.', new)
    if 'power and design issues in crossover-based' in t:
        new = re.sub(r'Wang, Yan\b(?!pin)', 'Wang, Yanpin', new)
    if 'a comparison of four methods for the analysis of n-of-1' in t:
        new = re.sub(r'\s+and\s+Xia,\s*Yinglin', '', new)
    if 'power analysis for idiographic' in t:
        new = re.sub(r'(\byear\s*=\s*)\{2022\}', r'\g<1>{2023}', new)
    if 'individual (n-of-1) trials can be combined' in t:
        new = re.sub(r'author\s*=\s*\{[^}]*\}',
                     'author = {Zucker, Deborah R. and Ruthazer, Robin and Schmid, Christopher H.}',
                     new, count=1)
    for key, doi in DOIS.items():
        if key in t and not re.search(r'(?i)\bdoi\s*=', new):
            new = re.sub(r'^(@\w+\s*\{[^,\n]*,)', r'\1\n  doi     = {' + doi + '},', new, count=1)
    if new != entry:
        head = entry.split('\n', 1)[0]
        log.append(head.strip())
    return new


def main():
    backup = sys.argv[1]
    os.makedirs(backup, exist_ok=True)
    for path in sys.argv[2:]:
        with open(path, encoding='utf-8') as f:
            text = f.read()
        log = []
        out, last = [], 0
        for s, e in entries(text):
            out.append(text[last:s])
            out.append(fix_entry(text[s:e], log))
            last = e
        out.append(text[last:])
        new = ''.join(out)
        if new != text:
            tag = path.replace('/', '_')
            shutil.copy2(path, os.path.join(backup, tag))
            tmp = os.path.join(backup, tag + '.new')
            with open(tmp, 'w', encoding='utf-8') as f:
                f.write(new)
            shutil.move(tmp, path)
            print(f'{path}: {len(log)} entries changed')
            for h in log:
                print(f'    {h}')


if __name__ == '__main__':
    main()
