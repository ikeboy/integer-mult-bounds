#!/usr/bin/env python3
"""PR117's h=24 complex DAG with lifted and shrunk binary frames: the per-core row histogram for lockstep shared cores.

Reads PR117's DAG from PR129's pinned copy (SHA-256 below, the same bytes as PR117's), checks every root support exactly, rebuilds the maximum
carrier matching under PR104's inherited adjacency (71185 links, R = 28705) and PR110's explicit role compile, then
chooses the frame of every addition:
  lifted: U_x = common F2 kernel of the functionals of every root reachable along slot successors (PR110);
  shrunk: U_x is replaced by Lb_x = label(x) + (frames of every role arriving at x) when Lb_x is nondegenerate,
          strictly smaller, and the one-child cost d*ln(m/d) of the chain steps through x drops (PR125's rule, here
          without deferral: every role starts at frame 0 and ends at F inside each core).
Checks: B. the compiled word computes every root value; C. the transparent word (readouts at 0, V, L, centre and
output reads, L^-1, V^-1) restores arbitrary scratch and adds exactly x_S to every target (two seeds, mod 2^61-1);
D. exactly over F2, every role chain 0 -> F is nested and every frame on it nondegenerate, every label lies in its
frame, every root frame contains the last node frame; E. the per-core row histogram has rank 24 R + 552.
Writes complex-row.json next to this file.  Usage: python3 research/lockstep-cores/complex_shrunk.py
Prepared by Avi Eisenberg (ikeboy) with Anthropic Claude assistance. Apache-2.0.
"""
import array
import gzip
import hashlib
import math
import json
import random
import struct
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
DAG = ROOT / 'research' / 'cyclic-deferred' / 'inputs' / 'complex-dag.json.gz'
DAG_SHA256 = '3c034d0aae388ef567a454826f4f48b26fd8a94c71e8ffed4835271b349a783b'

P = (1 << 61) - 1


def require(ok, msg):
    if not ok:
        raise SystemExit('FAIL: ' + msg)


def load117(path):
    raw = Path(path).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == DAG_SHA256, 'pinned PR117 DAG hash')
    d = json.loads(gzip.decompress(raw))
    h, v = d['h'], d['v']; trip = list(combinations(range(h), 3)); fa = d['args']; c = len(fa) // 2
    require((h, v, c, d['central_disjoint'], d['center_denominator']) == (24, 2024, 91770, 24, 21), 'PR117 header')
    n = v + 1 + c
    args = [(0, 0)] * n; sup = [0] * n; core = [0] * n; cover = [0] * n
    for i, T in enumerate(trip):
        sup[i + 1] = 1 << i; core[i + 1] = cover[i + 1] = sum(1 << p for p in T)
    for i in range(c):
        a, b = fa[2 * i], fa[2 * i + 1]; x = v + 1 + i
        require(0 < a < x and 0 < b < x and not sup[a] & sup[b], 'disjoint addition %d' % x)
        args[x] = (a, b); sup[x] = sup[a] | sup[b]; core[x] = core[a] & core[b]; cover[x] = cover[a] | cover[b]
    pm = [sum(1 << j for j, t in enumerate(trip) if a in t) for a in range(h)]; full = (1 << v) - 1
    lab = {}
    for j, r in enumerate(d['D']):
        T = trip[j]; require(sup[r] == full & ~(pm[T[0]] | pm[T[1]] | pm[T[2]]), 'disjoint root %d' % j)
        lab[('A',) + T] = r
    for r in d['P']:
        s_ = sup[r]; pts = [a for a in range(h) if s_ & pm[a] == s_]
        require(len(pts) == 2, 'pair root core')
        a, b = pts; miss = [i for i in range(h) if i not in (a, b) and not s_ & pm[i] & pm[a] & pm[b]]
        require(len(miss) == 1 and s_ == pm[a] & pm[b] & ~pm[miss[0]], 'pair root support')
        lab[('P', a, b, miss[0])] = r
    for i, r in enumerate(d['A']):
        require(sup[r] == full & ~pm[i], 'centre root %d' % i); lab[('C', i)] = r
    roots = [lab[('A',) + T] for T in trip]
    for a, b in combinations(range(h), 2):
        for i in sorted((i for i in range(h) if i not in (a, b)), key=lambda i: ((i ^ 1) in (a, b), i)):
            roots.append(lab[('P', a, b, i)])
    roots += [lab[('C', i)] for i in range(h)]
    require(len(roots) == len(lab) == 8120, 'root count')
    kind = [0] * (8120 - h) + [1] * h
    types = [0] * n; ranks = [0] * n
    for x in range(1, n):
        if x <= v: types[x], ranks[x] = 1, 1
        elif bin(core[x]).count('1') >= 2: types[x], ranks[x] = 1, bin(sup[x]).count('1')
        else: types[x], ranks[x] = 2, bin(cover[x]).count('1')
    active = [0] * n; st = list(roots)
    while st:
        x = st.pop()
        if x and not active[x]:
            active[x] = 1
            if args[x][0]: st.extend(args[x])
    for x in range(v + 1, n):                      # PR104 label nesting: every operand label inside its node label
        if active[x]:
            for y in args[x]:
                if types[y] == 1 and types[x] == 1: require(not (core[x] & ~core[y] or cover[y] & ~cover[x]), 'nest %d' % x)
                else: require(types[x] == 2 and not cover[y] & ~cover[x], 'nest %d' % x)
    return h, v, n, len(roots), args, core, cover, roots, kind, active, ranks, types


def load(prefix):
    def rd(fh, typ, cnt):
        a = array.array(typ); a.fromfile(fh, cnt); return a
    with open(str(prefix) + '.bin', 'rb') as f:
        h, v, n, q = struct.unpack('<4I', f.read(16))
        flat = rd(f, 'I', 2 * n); core = rd(f, 'Q', n); cover = rd(f, 'Q', n)
        roots = rd(f, 'I', q); kind = rd(f, 'I', q); active = rd(f, 'B', n)
    with open(str(prefix) + '.labels', 'rb') as f:
        ranks = rd(f, 'I', n); types = rd(f, 'B', n)
    args = [(flat[2 * i], flat[2 * i + 1]) for i in range(n)]
    return h, v, n, q, args, core, cover, roots, kind, active, ranks, types


# ---------------------------------------------------------------- F2 linear algebra on bit vectors
def dot(a, b): return bin(a & b).count('1') & 1


def reduce(basis):
    rows = []; piv = []
    for vec in basis:
        for p, r in zip(piv, rows):
            if vec >> p & 1: vec ^= r
        if vec:
            p = vec.bit_length() - 1
            for i in range(len(rows)):
                if rows[i] >> p & 1: rows[i] ^= vec
            rows.append(vec); piv.append(p)
    return rows


def kernel(funcs, h):
    K = reduce(funcs); piv = [r.bit_length() - 1 for r in K]
    out = []
    for c in (i for i in range(h) if i not in piv):
        vec = 1 << c
        for r, p in zip(K, piv):
            if r >> c & 1: vec |= 1 << p
        out.append(vec)
    return out


def restrict(basis, funcs):
    B = reduce(basis)
    for f in funcs:
        if not B: break
        vals = [dot(b, f) for b in B]
        if not any(vals): continue
        k = vals.index(1); pv = B[k]
        B = [b ^ pv if val else b for b, val in zip(B, vals)]; B.pop(k)
    return reduce(B)


def contains(A, B):
    RB = reduce(B)
    return len(reduce(RB + list(A))) == len(RB)


def cap(A, B):
    A = reduce(A); B = reduce(B)
    if not A or not B: return []
    ker = []; red = []
    for i, vec in enumerate(A + B):
        tag = 1 << i
        for rv, rt in red:
            if vec >> (rv.bit_length() - 1) & 1: vec ^= rv; tag ^= rt
        if vec: red.append((vec, tag))
        else:
            comb_ = 0
            for k in range(len(A)):
                if tag >> k & 1: comb_ ^= A[k]
            ker.append(comb_)
    return reduce(ker)


def nondeg(B):
    B = reduce(B); k = len(B)
    if not k: return True
    M = [sum(dot(B[i], B[j]) << j for j in range(k)) for i in range(k)]
    return len(reduce(M)) == k


def hopcroft_karp(left, adj):
    """maximum bipartite matching; adj[x] = sorted list of right vertices."""
    INF = 1 << 60
    mate_l = {x: None for x in left}; mate_r = {}
    while True:
        dist = {}; queue = []
        for x in left:
            if mate_l[x] is None: dist[x] = 0; queue.append(x)
            else: dist[x] = INF
        found = INF; head = 0
        while head < len(queue):
            x = queue[head]; head += 1
            if dist[x] >= found: continue
            for r in adj[x]:
                y = mate_r.get(r)
                if y is None: found = min(found, dist[x] + 1)
                elif dist[y] == INF: dist[y] = dist[x] + 1; queue.append(y)
        if found == INF: break
        gained = 0
        for x0 in left:
            if mate_l[x0] is not None: continue
            stack = [(x0, iter(adj[x0]))]; path = []
            while stack:
                x, it = stack[-1]
                advanced = False
                for r in it:
                    y = mate_r.get(r)
                    if (y is None and dist[x] + 1 == found) or (y is not None and dist[y] == dist[x] + 1):
                        path.append((x, r))
                        if y is None:
                            for xx, rr in path: mate_l[xx] = rr; mate_r[rr] = xx
                            gained += 1; stack = []; break
                        stack.append((y, iter(adj[y]))); advanced = True; break
                else:
                    dist[x] = INF; stack.pop()
                    if path: path.pop()
                    continue
                if not advanced: break
        if not gained: break
    return mate_l


def main():
    h, v, n, q, args, core, cover, roots, kind, active, ranks, types = load117(DAG)
    trip = list(combinations(range(h), 3)); require(len(trip) == v, 'triple count')
    m = h * h; N = v * v; FULL = (1 << h) - 1
    tmask = [sum(1 << p for p in T) for T in trip]
    tid = {T: i for i, T in enumerate(trip)}

    # ------------------------------------------------------------ uses and matching (inherited adjacency)
    degree = [0] * n; uses = defaultdict(list); c_add = 0
    for x in range(1, n):
        if active[x] and args[x][0]:
            c_add += 1
            for pos, y in enumerate(args[x]): degree[y] += 1; uses[y].append(('gate', x, pos))
    for j, x in enumerate(roots): degree[x] += 1; uses[x].append(('root', j))
    def unode(u): return roots[u[1]] if u[0] == 'root' else u[1]
    def before(x, u):
        y = unode(u)
        return (ranks[x], x) < (ranks[y], (n + u[1]) if u[0] == 'root' else y)
    def incl(x, y):
        tx, ty = types[x], types[y]
        if tx == 1 and ty == 1: return not (core[y] & ~core[x]) and not (cover[x] & ~cover[y])
        if tx in (1, 2) and ty == 2: return not (cover[x] & ~cover[y])
        if tx == 1 and ty == 3: return bool(core[x] & core[y])
        if tx == 2 and ty == 3: return not (cover[x] & ~core[y])
        if tx == 3 and ty == 2: return cover[y] == FULL
        if tx == 3 and ty == 3: return core[x] == core[y]
        return False
    left = []; adj = {}
    for x in range(1, n):
        if not (active[x] and args[x][0]): continue
        rs = sorted({(y, k) for y in args[x] for k, u in enumerate(uses[y]) if before(x, u) and incl(x, unode(u))})
        if rs: left.append(x); adj[x] = rs
    mate = hopcroft_karp(left, adj)
    links = {x: r for x, r in mate.items() if r is not None}
    R = c_add + q - len(links)
    require(len(links) == 71185 and R == 28705, 'matching size %d, roles %d' % (len(links), R))
    linked_use = {(y, k): x for x, (y, k) in links.items()}

    # ------------------------------------------------------------ root targets, read coefficients, root functionals
    target = [None] * q
    for j in range(v): target[j] = j
    idx = v
    for a, b in combinations(range(h), 2):
        others = sorted((i for i in range(h) if i not in (a, b)), key=lambda i: ((i ^ 1) in (a, b), i))
        for i in others: target[idx] = tid[tuple(sorted((a, b, i)))]; idx += 1
    centre_of = {}
    for j in range(q):
        if kind[j]:
            miss = FULL & ~cover[roots[j]]; require(bin(miss).count('1') == 1, 'centre cover')
            centre_of[j] = miss.bit_length() - 1
    require(idx + len(centre_of) == q, 'root order')
    rootfun = [(1 << centre_of[j]) if kind[j] else tmask[target[j]] for j in range(q)]

    # ------------------------------------------------------------ lifted binary frames, monotone
    succ = defaultdict(set); droot = defaultdict(set)
    for x in range(1, n):
        if not active[x]: continue
        for k, u in enumerate(uses[x]):
            if (x, k) in linked_use: continue
            (succ[x].add(u[1]) if u[0] == 'gate' else droot[x].add(u[1]))
    for x, (y, k) in links.items():
        u = uses[y][k]; (succ[x].add(u[1]) if u[0] == 'gate' else droot[x].add(u[1]))
    order_desc = sorted((x for x in range(1, n) if active[x]), key=lambda x: (ranks[x], x), reverse=True)
    K = {}
    for x in order_desc:
        rows = [rootfun[j] for j in droot[x]]
        for t in succ[x]: rows += K[t]
        K[x] = reduce(rows)
    def envelope(x):
        if not args[x][0]: return [tmask[x - 1]]
        if types[x] == 2: return [1 << i for i in range(h) if cover[x] >> i & 1]
        require(types[x] == 1, 'only ordinary and common-pair labels occur')
        ks = [i for i in range(h) if (cover[x] & ~core[x]) >> i & 1]
        require(len(ks) == ranks[x], 'common-pair support')
        return [core[x] | (1 << k) for k in ks]
    U = {}
    for x in range(1, n):
        if not active[x]: continue
        if not args[x][0]: U[x] = envelope(x); continue
        B = kernel(K[x], h)
        require(contains(envelope(x), B), 'label outside lifted frame')
        U[x] = B if nondeg(B) else envelope(x)
    changed = True
    while changed:
        changed = False
        for x in order_desc:
            if args[x][0] and any(not contains(U[x], U[t]) for t in succ[x]):
                env = envelope(x)
                if len(reduce(U[x])) != len(reduce(env)) or not contains(U[x], env): U[x] = env; changed = True
    dimU = {x: len(reduce(B)) for x, B in U.items()}

    # ------------------------------------------------------------ explicit role compile
    order = sorted((x for x in range(1, n) if active[x]), key=lambda x: (ranks[x], x))
    holds = []; first_node = []; ops = []; edge_key = {}; role_root = {}
    def new_role(node): holds.append([node]); first_node.append(node); return len(holds) - 1
    def serve(s, y, k):
        u = uses[y][k]
        if u[0] == 'gate': edge_key[(y, k)] = s
        else: role_root[s] = u[1]
    for x in order:
        free = [k for k in range(len(uses[x])) if (x, k) not in linked_use]
        if not args[x][0]:
            for k in free:
                s = new_role(x); ops.append(('src', s, x)); serve(s, x, k)
            continue
        a, b = args[x]
        ka = next(k for k, u in enumerate(uses[a]) if u == ('gate', x, 0))
        kb = next(k for k, u in enumerate(uses[b]) if u == ('gate', x, 1))
        sa, sb = edge_key.pop((a, ka)), edge_key.pop((b, kb))
        if x in links and links[x][0] == a: sa, sb = sb, sa
        piv, oth = sa, sb
        ops.append(('add', piv, oth, x)); holds[piv].append(x); holds[oth].append(x)
        if x in links: serve(oth, *links[x])
        require(free, 'every use linked at node %d' % x)
        serve(piv, x, free[0])
        for k in free[1:]:
            f = new_role(x); ops.append(('copy', piv, f, x)); serve(f, x, k)
    require(not edge_key and len(holds) == R, 'compile')
    Rr = R
    def fcost(d): return d * math.log(m / d) if d > 0 else 0.0
    through = defaultdict(list)
    for s in range(Rr):
        for i, x in enumerate(holds[s]):
            if args[x][0]: through[x].append((s, i))
    def nextdim(s, i):
        return dimU[holds[s][i + 1]] if i + 1 < len(holds[s]) else (h - 1 if s in role_root else h)
    nshrink = 0
    for x in order:
        if not args[x][0]: continue
        Lb = reduce([vec for s, i in through[x] if i > 0 for vec in U[holds[s][i - 1]]] + envelope(x))
        du = dimU[x]; dl = len(Lb)
        if dl >= du or not nondeg(Lb): continue
        old = new = 0.0
        for s, i in through[x]:
            dp = dimU[holds[s][i - 1]] if i > 0 else 0
            dn = nextdim(s, i)
            old += fcost(du - dp) + fcost(dn - du); new += fcost(dl - dp) + fcost(dn - dl)
        if new < old - 1e-9:
            U[x] = Lb; dimU[x] = dl; nshrink += 1

    # ------------------------------------------------------------ B. the word computes every root value
    rng = random.Random(20261008)
    xs = [rng.randrange(P) for _ in range(v)]
    val = [0] * n
    for i in range(v): val[i + 1] = xs[i]
    for x in order:
        if args[x][0]: val[x] = (val[args[x][0]] + val[args[x][1]]) % P
    leaf_of = {s: first_node[s] for s in range(Rr) if not args[first_node[s]][0]}
    a_ = [0] * Rr
    for s, leaf in leaf_of.items(): a_[s] = xs[leaf - 1]
    for o in ops:
        if o[0] == 'add': a_[o[1]] = (a_[o[1]] + a_[o[2]]) % P
        elif o[0] == 'copy': a_[o[2]] = (a_[o[2]] + a_[o[1]]) % P
    require(all(a_[s] == val[roots[j]] for s, j in role_root.items()), 'root values')

    # ------------------------------------------------------------ phase one and garbage reach
    def touch(o): return [o[1], o[2]] if o[0] in ('add', 'copy') else []
    prev = {}; pred = defaultdict(list); last = {}
    for i, o in enumerate(ops):
        for s in touch(o):
            if s in prev: pred[i].append(prev[s])
            prev[s] = i; last[s] = i
    centre_roles = [s for s, j in role_root.items() if kind[j]]
    Anc = set(); st = [last[s] for s in centre_roles]
    while st:
        i = st.pop()
        if i not in Anc: Anc.add(i); st.extend(pred[i])
    touched = set(centre_roles)
    for i in Anc: touched.update(touch(ops[i]))
    require(all(last[s] in Anc for s in centre_roles), 'centres complete in phase one')
    reach = [set() for _ in range(Rr)]; reach_all = [False] * Rr
    for s, j in role_root.items():
        if kind[j]: reach_all[s] = True
        else: reach[s].add(target[j])
    for o in reversed(ops):
        if o[0] == 'add': reach[o[2]] |= reach[o[1]]; reach_all[o[2]] = reach_all[o[2]] or reach_all[o[1]]
        elif o[0] == 'copy': reach[o[1]] |= reach[o[2]]; reach_all[o[1]] = reach_all[o[1]] or reach_all[o[2]]

    # ------------------------------------------------------------ no deferral inside a shared core
    def F0(s): return U[first_node[s]]
    placed = {}; deferred = []; dset = set()

    # ------------------------------------------------------------ C. replay with arbitrary scratch and data
    inv = lambda a: pow(a % P, P - 2, P); HALF = inv(2); I21 = inv(21)
    cvec = [None] * Rr; dpart = [dict() for _ in range(Rr)]
    for s, j in role_root.items():
        if kind[j]: c = [0] * h; c[centre_of[j]] = 1; cvec[s] = c
        else: dpart[s][target[j]] = (P - HALF) if j >= v else HALF
    seed_c = {s: (None if cvec[s] is None else list(cvec[s])) for s in role_root}
    seed_d = {s: dict(dpart[s]) for s in role_root}
    def addc(dst, src):
        if cvec[src] is not None:
            cvec[dst] = list(cvec[src]) if cvec[dst] is None else [(p + q_) % P for p, q_ in zip(cvec[dst], cvec[src])]
        for t, c in dpart[src].items(): dpart[dst][t] = (dpart[dst].get(t, 0) + c) % P
    for o in reversed(ops):
        if o[0] == 'add': addc(o[2], o[1])
        elif o[0] == 'copy': addc(o[1], o[2])
    readout_terms = sum(len(dpart[s]) + (v if cvec[s] is not None and any(cvec[s]) else 0) for s in range(Rr))
    readout_terms += sum(v if kind[j] else 1 for s, j in role_root.items())
    def readout(y, s, value, sign, seed=False):
        cv = seed_c[s] if seed else cvec[s]; dp = seed_d[s] if seed else dpart[s]
        if cv is not None and any(cv):
            tot = sum(cv) % P
            for t, (i0, i1, i2) in enumerate(trip):
                y[t] = (y[t] + sign * value * (I21 * tot - HALF * (cv[i0] + cv[i1] + cv[i2]))) % P
        for t, c in dp.items(): y[t] = (y[t] + sign * c * value) % P
    phase1 = sorted(Anc); rest = [i for i in range(len(ops)) if i not in Anc]
    def replay(seed):
        rng = random.Random(seed)
        x = [rng.randrange(P) for _ in range(v)]; z = [rng.randrange(P) for _ in range(Rr)]
        y0 = [rng.randrange(P) for _ in range(v)]; a = list(z); y = list(y0)
        for s in range(Rr):
            if s not in dset: readout(y, s, a[s], -1)
        for s, leaf in leaf_of.items():
            if s not in dset: a[s] = (a[s] + x[leaf - 1]) % P
        def run(i, sign=1):
            o = ops[i]
            if o[0] == 'add': a[o[1]] = (a[o[1]] + sign * a[o[2]]) % P
            elif o[0] == 'copy': a[o[2]] = (a[o[2]] + sign * a[o[1]]) % P
        for i in phase1: run(i)
        for s in centre_roles: readout(y, s, a[s], +1, True)
        for s in deferred: readout(y, s, a[s], -1)
        for s in deferred:
            if s in leaf_of: a[s] = (a[s] + x[leaf_of[s] - 1]) % P
        for i in rest: run(i)
        for s, j in role_root.items():
            if not kind[j]: readout(y, s, a[s], +1, True)
        for i in reversed(range(len(ops))): run(i, -1)
        for s, leaf in leaf_of.items(): a[s] = (a[s] - x[leaf - 1]) % P
        return a == z, all((y[t] - y0[t] - x[t]) % P == 0 for t in range(v))
    rep = [replay(seed) for seed in (1, 2)]
    require(all(r == (True, True) for r in rep), 'replay %s' % rep)

    # ------------------------------------------------------------ D. exact F2 frame facts
    root_frame = {}
    for s, j in role_root.items():
        root_frame[s] = [1 << i for i in range(h) if i != centre_of[j]] if kind[j] else kernel([tmask[target[j]]], h)
    FULLB = [1 << i for i in range(h)]
    chain_dims = []
    for s in range(Rr):
        seq = [placed.get(s, []), F0(s)] + [U[x] for x in holds[s][1:]]
        if s in root_frame: seq.append(root_frame[s])
        seq.append(FULLB)
        require(all(contains(A, B) for A, B in zip(seq, seq[1:])), 'role chain nesting %d' % s)
        require(all(nondeg(B) for B in seq[1:-1] if B), 'degenerate frame on role %d' % s)
        chain_dims.append([len(reduce(B)) for B in seq])
    for x in range(1, n):
        if active[x] and args[x][0]:
            require(contains(envelope(x), U[x]) and nondeg(U[x]), 'label or frame at node %d' % x)
    for s, X in placed.items():
        require(contains(X, F0(s)) and nondeg(X), 'deferral frame %d' % s)
        require(all(not any(dot(xv, tmask[t]) for xv in X) for t in reach[s]), 'target frame %d' % s)
    byT2 = defaultdict(list)
    for s in deferred:
        for t in reach[s]: byT2[t].append(s)
    for t, ss in byT2.items():
        ss.sort(key=lambda s: (len(placed[s]), s))
        require(all(contains(placed[p], placed[q_]) for p, q_ in zip(ss, ss[1:])), 'target chain %d' % t)

    # ------------------------------------------------------------ E. per-core row histogram
    H = Counter()
    for s in range(Rr):
        ds = chain_dims[s]
        for a, b in zip(ds, ds[1:]):
            if b > a: H[b - a] += 1                                            # one child per strict increase
        if s in role_root and kind[role_root[s]]: H[h - 1] += 1                # copied centre transform
    require(sum(t * c for t, c in H.items()) == h * R + h * (h - 1), 'row rank 24R + 552')
    out = dict(h=h, v=v, additions=c_add, roots=q, links=len(links), R=R, m=m, N=N, shrunk_frames=nshrink,
               readout_terms=readout_terms, dag_sha256=DAG_SHA256, row_rank=sum(t * c for t, c in H.items()),
               replay=dict(seeds=[1, 2], scratch_restored=True, y_plus_x=True, field='Z/(2^61-1)'),
               row_histogram=dict(sorted(H.items())))
    (HERE / 'complex-row.json').write_text(json.dumps(out, indent=1) + '\n')
    print('PASS complex row: R=%d, links %d, shrunk frames %d, row children %d, row rank %d' % (
        R, len(links), nshrink, sum(H.values()), out['row_rank']))


if __name__ == '__main__':
    main()
