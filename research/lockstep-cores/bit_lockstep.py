#!/usr/bin/env python3
"""Bit side (h = 23) with PR129's shared cores, executed in lockstep.

Inputs, pinned by SHA-256 (all in this branch):
  PR128's class-split profile of its minimal-V bit word (one child per residual, PR104's rule), per core and x 2v;
  PR129's partition of the 1771 triples into 161 sunflower classes of 11;
  PR129's consecutive shared bit profile (only as a cross-check).
Checks: the partition covers every triple once, and two triples of one class meet in exactly one point, so their
outer lines are orthogonal for H = 9I - J (t_T^T H t_S = 9|T cap S| - 9 = 0, t_T^T H t_T = 18). The consecutive
schedule rebuilt here reproduces PR129's shared profile exactly.
Lockstep: the 11 cores of a class run the same inner word together. On a shared scratch stream, their frame steps
U -> U' are the orthogonal residuals (U' - U) x <t_b>, b in the class; together they are one projector residual on
(U' - U) x S_J of rank 11 d, hence one child of width 11 d (PR104's rule). Copied-centre transforms merge the same
way. Source and target data streams, data bridges and endpoint copies are separate streams per core and stay as
they are. Each group-role keeps PR129's complement of rank 529 - 11 (23 - dim sigma).
Writes lockstep-bit-profile.json.  Usage: python3 research/lockstep-cores/bit_lockstep.py
Prepared by Avi Eisenberg (ikeboy) with Anthropic Claude assistance. Apache-2.0.
"""
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BS = ROOT / 'research' / 'cyclic-deferred' / 'bit-sharing'
PINS = {BS / 'inputs' / 'research' / 'round8-network-opposite-profile.json':
        '8403eb58f658999584de3f1e1af56e0ddb90ce0521ae6a68f9f57864da0dcdd1',
        BS / 'groups161.json': '220dc207ab0774212949a57301391506ead8bf9f6eadd60a6a82ec2de7e0bc61'}
SHARED = ROOT / 'research' / 'cyclic-deferred' / 'shared-bit-profile.json'


def require(ok, msg):
    if not ok:
        raise SystemExit('FAIL: ' + msg)


def main():
    for path, digest in PINS.items():
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'pinned hash of ' + path.name)
    base = json.loads((BS / 'inputs' / 'research' / 'round8-network-opposite-profile.json').read_text())
    groups = [[tuple(T) for T in G] for G in json.loads((BS / 'groups161.json').read_text())]
    h, m, v, R, N = base['h'], base['m'], base['v'], base['R'], base['N']
    require((h, m, v, R, N) == (23, 529, 1771, 28866, 1771 ** 2), 'bit header')
    trip = list(combinations(range(h), 3))
    require(len(groups) == 161 and all(len(G) == 11 for G in groups), '161 classes of 11')
    flat = [tuple(sorted(T)) for G in groups for T in G]
    require(len(flat) == len(set(flat)) == len(trip) and set(flat) == set(trip), 'partition of all triples')
    for G in groups:
        for T, S in combinations(G, 2):
            require(9 * len(set(T) & set(S)) - 9 == 0, 'H-orthogonal outer lines in a class')
    K, Gc = 11, len(groups)
    cls = {k: Counter({int(t): c for t, c in val['hist'].items()}) for k, val in base['classes'].items()}
    require(set(cls) == {'auxiliary_steps', 'source_steps', 'target_steps', 'copied_centres', 'auxiliary_exteriors',
                         'data_bridges', 'endpoint_copies'}, 'class names')
    sig = Counter()
    for t, c in cls['auxiliary_exteriors'].items():
        require(c % (2 * v) == 0, 'exterior multiplicity'); sig[t - (m - h)] += c // (2 * v)
    require(sum(sig.values()) == R, 'one exterior per role and core')
    L = 2 * v * h * (h - 1); W = 2 * N + 2 * Gc * R

    def profile(lock):
        z = Counter()
        for name in ('auxiliary_steps', 'copied_centres'):
            for t, c in cls[name].items():
                require(c % (2 * v) == 0, 'per-core multiplicity')
                if lock: z[K * t] += 2 * Gc * (c // (2 * v))
                else: z[t] += c
        for name in ('source_steps', 'target_steps', 'data_bridges', 'endpoint_copies'):
            z.update(cls[name])
        for f, c in sig.items(): z[m - K * (h - f)] += 2 * Gc * c
        z = Counter({t: c for t, c in z.items() if c})
        require(sum(t * c for t, c in z.items()) == W * m - N + L, 'rank mass')
        require(all(0 < t < m for t in z), 'children below m')
        return z

    seq = profile(False)
    shared = json.loads(SHARED.read_text())
    require({int(t): c for t, c in shared['child_multiplicities'].items()} == dict(seq) and shared['W'] == W,
            'consecutive schedule reproduces PR129 shared-bit-profile.json')
    lock = profile(True)
    out = dict(h=h, m=m, v=v, N=N, R=R, L=L, W=W, groups=Gc, group_size=K, total_rank=W * m - N + L,
               deficit=N - L, maxchild=max(lock), schedule='lockstep: the 11 cores of a class run the same inner word together',
               merged_classes=['auxiliary_steps', 'copied_centres'],
               per_core_classes=['source_steps', 'target_steps', 'data_bridges', 'endpoint_copies'],
               sigma_dims=dict(sorted(sig.items())), consecutive_matches_pr129=True,
               child_multiplicities=dict(sorted(lock.items())))
    (HERE / 'lockstep-bit-profile.json').write_text(json.dumps(out, indent=1) + '\n')
    print('PASS bit: 161 classes checked, PR129 consecutive profile reproduced, lockstep W=%d maxchild %d' % (W, max(lock)))


if __name__ == '__main__':
    main()
