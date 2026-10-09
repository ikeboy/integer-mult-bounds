# Lockstep shared cores

Conditional κ = 1028671/5000000000 = **2.057342 × 10^-4** (+59.3% over PR #129's 1.29186205798 × 10^-4).

PR #128 and PR #129 let the cores of an orthogonal group share scratch streams, and run those cores one after
another. The cores of a group run the same inner word on disjoint fibres and orthogonal slices, so they can also run
in lockstep. Then every frame step of a shared stream is one projector residual of rank |J|·d instead of |J|
residuals of rank d, and under PR #104's one-child rule it is one child of width |J|·d. Nothing physical changes:
the same gates, frames, streams and scalar work, in a different order.

| | PR #129 | here |
|---|---:|---:|
| stopped bit saving (binds) | 1.292196e-4 | 2.058190e-4 |
| complex saving | 1.306965e-4 | 3.083799e-4 |
| κ | 1.291862e-4 | 2.057342e-4 |

See [PROOF.md](PROOF.md) for the schedule, the merge argument, what is checked and the inherited assumptions.

## Reproduce

From the repository root (standard-library Python 3; about two minutes):

```sh
python3 research/lockstep-cores/complex_shrunk.py
python3 research/lockstep-cores/bit_lockstep.py
python3 research/lockstep-cores/certificate.py
```

## Files

| file | content |
|---|---|
| `complex_shrunk.py` | PR #117's DAG, inherited matcher, role compile, lifted and shrunk frames, replay, exact F2 checks, per-core row histogram |
| `bit_lockstep.py` | PR #129's 161 classes and PR #128's class-split bit profile (pinned); reproduces PR #129's shared profile; lockstep bit profile |
| `certificate.py` | lockstep complex profile (PR #128's 87 groups), exact moments, stopped bit saving, finite bridge, 47-constraint assembly, κ |
| `complex-row.json`, `lockstep-bit-profile.json`, `certificate.json` | generated outputs |
