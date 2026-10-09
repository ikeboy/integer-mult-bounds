# Lockstep shared cores

## 1. Claim

Under every interface that PR #128 and PR #129 retain (PR #104's stopped product-ring interchange and one-child
rule, the transparent dirty-core contract, translated endpoints, copied centres, the common generic basis, the
stopped-bit transfer, routing, recovery, the fixed-tape interfaces and the analytic estimates),

κ = 1028671/5000000000 = 2.057342 × 10^-4,

which is +59.3% over PR #129 (1.29186205798 × 10^-4). The stopped bit saving binds. The next κ grid point at
denominator 10^10 is rejected. This is a finite conditional witness.

| | PR #129 (cores run one after another) | here (cores run in lockstep) |
|---|---:|---:|
| bit: coarse saving | 1.293104e-4 | 2.059865e-4 |
| bit: stopped saving | 1.292196e-4 | 2.058190e-4 |
| complex saving | 1.306965e-4 | 3.083799e-4 |
| κ | 1.291862e-4 | 2.057342e-4 |

The only new ingredient is the order in which the cores of a shared group execute. No physical word, frame, gate,
stream allocation, partition or scalar operation is added or removed.

## 2. Setting: completed-core sharing

In each stage there are v cores, one per outer triple label b. Core b runs the same inner word (the same DAG, roles,
frames, gates and readout coefficients) on its own data fibre. PR #128 (complex, h = 24) and PR #129 (bit, h = 23)
partition the labels into groups J whose outer lines are orthogonal: PR #128's 83 groups of 24 and 4 groups of 8,
with binary Gram matrix I, and PR #129's 161 sunflower classes of 11, with t_T^T H t_S = 9|T ∩ S| − 9 = 0 for
H = 9I − J. Each role u of a group gets one scratch stream. Core b uses the slice X ⊗ ⟨t_b⟩ of that stream's address
space, and the slices of a group are orthogonal. PR #128 and PR #129 run the cores of a group one after another.

## 3. The lockstep schedule

Run the cores of each group together, gate by gate. Because every core executes the same inner word, at every point
of the word the cores agree on which gate comes next, which streams it touches and at which frames.

* A gate of core b acts on the slice of core b in each shared stream, and on core b's own data and target streams,
  at exactly the common frames it has in the consecutive schedule. Per slice, nothing changes.
* The cores of a group touch disjoint data fibres and disjoint (orthogonal) slices of the shared streams, so their
  gates commute and any interleaving computes the same values. Lockstep is one such interleaving.
* The joint core is the tensor product of identical transparent words on disjoint fibres. Its action is L_b on the
  data of every core b, and on each shared stream it leaves the same total residual that the consecutive schedule
  leaves after the last core of the group (C_{G_J} in PR #128's notation). Group complements are therefore unchanged.

## 4. Why a merged frame step is one child

In lockstep, all cores of a group raise the frame of a shared stream at the same time, from U ⊗ ⟨t_b⟩ to
U' ⊗ ⟨t_b⟩ for every b in J. The residuals (U' − U) ⊗ ⟨t_b⟩ lie on mutually orthogonal subspaces:

* complex side: C_X C_Y = C_{X⊕Y} for orthogonal nondegenerate X, Y, with the exact Gaussian phase (PR #128, §2);
* bit side: the rank-one projectors of the triples of a class annihilate each other under H (PR #129).

So the product of the |J| residuals is the single nested projector residual from U ⊗ S_J to U' ⊗ S_J, where S_J is
the span of the group's outer lines. Its rank is |J|·d, with d = dim U' − dim U. Under PR #104's rule it is one child
of width |J|·d, by the same argument as `C_r ⊗ C_f = C_{rf}` in `notes/stopped-product-factorization.tex`. In the
consecutive schedule the same operator is applied as |J| separate residuals at different times, hence |J| children
of width d.

Copied-centre transforms are temporaries on their own streams and merge in the same way.

## 5. What stays per core

Source and target data streams belong to one core each. Their frame chains, the data bridges, the data fronts and
the endpoint copies stay exactly as in PR #128 and PR #129. Each group-role keeps its complement: width 384 for the
four partial complex groups (PR #128) and width 529 − 11(23 − dim σ) on the bit side (PR #129). The stream count W
and the total rank are unchanged: merging keeps every residual's rank and only changes how the ranks are grouped.

## 6. The finite profiles

**Bit side** (`bit_lockstep.py`). PR #128's class-split profile of its minimal-V word (one child per residual) and
PR #129's 161 classes, both pinned by SHA-256. The script checks that the classes partition all 1771 triples and are
pairwise H-orthogonal, and it rebuilds PR #129's consecutive shared profile exactly (`shared-bit-profile.json`). In
lockstep the per-core auxiliary steps and copied-centre transforms of the 11 cores of a class merge into one child
of width 11·r. W = 2N + 2·161·R = 15,567,734; the largest child is 518 < 529.

**Complex side** (`complex_shrunk.py`). PR #117's DAG (PR #129's pinned copy, the same bytes as PR #117's), PR #104's
inherited matcher (71,185 links, R = 28,705), PR #110's explicit role compile and lifted frames, and PR #125's frame
shrinking rule without deferral: every role climbs from frame 0 to F in every core, as PR #128's contract requires.
The per-core row histogram has rank 24R + 552. In lockstep, the row children of a group merge into children of width
24·r (full groups) or 8·r (partial groups). W = 2N + 2·87·R = 13,187,822; the largest child is 552 < 576. The
complex side does not bind.

## 7. What is checked

* `complex_shrunk.py`: every root support of PR #117's DAG; the matching size and R; with zero scratch the compiled
  word computes every root; the transparent word restores arbitrary scratch and adds exactly x_S to every target
  (two seeds, mod 2^61 − 1); exactly over F2, every role chain 0 → F is nested and nondegenerate and every label lies
  in its frame; the row rank 24R + 552.
* `bit_lockstep.py`: the pins, the partition and its H-orthogonality, the exact reproduction of PR #129's consecutive
  shared profile, and the rank mass W·m − N + L of the lockstep profile.
* `certificate.py`: the largest coarse bit and complex savings on the 10^-10 grid with PR #104's exact moment (next
  grid points rejected), the stopped bit saving (1 − 1/1000)·coarse + (1/1000)·384599/10^10, the finite bridge
  (its scalar term counts every DAG addition and every readout multiply-add) and the 47 strict constraints and 7
  margins of PR #104's assembly at the largest κ on the grid (next κ rejected).

Not machine-checked here: the lockstep schedule itself (§3–§4) is stated as a proof, not replayed literally; the
consecutive schedule's literal replays belong to PR #128 and PR #129.

## 8. Scope

All interfaces retained by PR #128 and PR #129 remain assumptions, including the opposite-bank factorization for
every residual in the finite list. That list now contains the merged residuals (U' − U) ⊗ S_J, so the common generic
basis is chosen for the enlarged finite family, as PR #129 already does for its projector family. The larger
children (bit 518, complex 552) change only the halving degrees, the row stock and the eventual cutoffs, not κ.
No global optimality is claimed.

## 9. Credits

* **an664**: completed-core sharing, the signed resolvable partition and the minimal-V bit word (PR #128).
* **eumemic**: the 161-class bit partition and shared bit cores (PR #129); the PR #117 complex DAG.
* **icekylinx**: the stopped product-ring interchange, one-child rule and rational centres (PR #104).
* **GamingPuzzled**: the frame shrinking rule (PR #125).
* **Swapnil Jain**: the round-seven bit word; **Zhihao Chen (jacklightChen)**: PR #97's ledgers and gauges.
* Everyone credited in the NOTICE and SOURCES files of PR #104, PR #128 and PR #129.

The lockstep schedule, its accounting and these scripts were prepared by Avi Eisenberg (ikeboy) with Anthropic
Claude assistance.
