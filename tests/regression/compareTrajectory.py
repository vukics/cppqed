#!/usr/bin/env python3
# Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
"""Regression check of a trajectory executable against a stored reference.

usage: compareTrajectory.py EXECUTABLE REFERENCE [--regenerate]

The reference file's first line is `# args: …`, the command-line arguments of the run.
The remaining lines are the streamed data lines and the `# OUTRO:` line of that run.
Data are compared numerically (complex numbers as `(re,im)`), and the jump trajectory
from the OUTRO line must agree in channels and, within tolerance, in times.
Header lines (banner, commit, paths) are ignored because they legitimately change.
"""
import json
import math
import shlex
import subprocess
import sys

RTOL, ATOL = 1e-8, 1e-12


def run(executable, args):
    out = subprocess.run([executable, *args], check=True, capture_output=True, text=True).stdout
    lines = [l for l in out.splitlines() if l.strip() and (not l.startswith('#') or l.startswith('# OUTRO: '))]
    return lines


def numbers(line):
    res = []
    for tok in line.split():
        if tok.startswith('('):
            re_, im = tok.strip('()').split(',')
            res += [float(re_), float(im)]
        else:
            res.append(float(tok))
    return res


def close(a, b):
    return math.isclose(a, b, rel_tol=RTOL, abs_tol=ATOL)


def compare(got, ref):
    errors = []
    data_got = [l for l in got if not l.startswith('#')]
    data_ref = [l for l in ref if not l.startswith('#')]
    if len(data_got) != len(data_ref):
        errors.append(f'{len(data_got)} data lines, reference has {len(data_ref)}')
    for i, (g, r) in enumerate(zip(data_got, data_ref)):
        ng, nr = numbers(g), numbers(r)
        if len(ng) != len(nr) or not all(map(close, ng, nr)):
            errors.append(f'data line {i}:\n  got {g}\n  ref {r}')
    outro = lambda lines: [json.loads(l[len('# OUTRO: '):]) for l in lines if l.startswith('# OUTRO: ')]
    og, or_ = outro(got), outro(ref)
    if len(og) != 1 or len(or_) != 1:
        errors.append('missing OUTRO line')
    else:
        jg, jr = og[0]['QJMC']['Jump trajectory'], or_[0]['QJMC']['Jump trajectory']
        if len(jg) != len(jr) or any(a[1] != b[1] or not close(a[0], b[0]) for a, b in zip(jg, jr)):
            errors.append(f'jump trajectory:\n  got {jg}\n  ref {jr}')
    return errors


def main():
    executable, reference = sys.argv[1], sys.argv[2]
    with open(reference) as f:
        header, *ref = f.read().splitlines()
    assert header.startswith('# args: '), 'reference must start with "# args: …"'
    args = shlex.split(header[len('# args: '):])
    got = run(executable, args)
    if '--regenerate' in sys.argv[3:]:
        with open(reference, 'w') as f:
            f.write('\n'.join([header, *got]) + '\n')
        print(f'regenerated {reference}')
        return 0
    errors = compare(got, ref)
    for e in errors:
        print(e, file=sys.stderr)
    print(f'{len(ref)} reference lines: ' + ('FAILED' if errors else 'agree'))
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
