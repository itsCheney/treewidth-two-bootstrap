"""Independent H-bootstrap percolation simulator (no shared code with the Lean project).

A round adds every missing pair {x,y} such that for some edge uv of H there is an
injective map f: V(H) -> V(G) with f(u)=x, f(v)=y and every other edge of H mapped
onto an edge of the *current* graph G (non-induced copies, all read the old graph).
"""
import itertools


def adj_from_edges(n, edges):
    A = [set() for _ in range(n)]
    for a, b in edges:
        A[a].add(b); A[b].add(a)
    return A


def order_for(H, n_h, u, v):
    """Vertex order starting with u, v then BFS-ish by max connectivity to placed vertices."""
    order = [u, v]
    rest = set(range(n_h)) - {u, v}
    while rest:
        best = max(rest, key=lambda w: (len(H[w] & set(order)), len(H[w])))
        order.append(best); rest.remove(best)
    return order


def new_edges(H, n_h, G, n_g, typed=None):
    """Return set of missing pairs witnessed in one round. typed: optional (htype, gtype) lists
    requiring type preservation."""
    found = set()
    edgesH = [(a, b) for a in range(n_h) for b in H[a] if a < b]
    for (u, v) in edgesH:
        order = order_for(H, n_h, u, v)
        # constraints: for position k, earlier neighbours in H-uv
        pos = {w: i for i, w in enumerate(order)}
        back = []
        for i, w in enumerate(order):
            nb = [x for x in H[w] if pos[x] < i and not ({w, x} == {u, v})]
            back.append(nb)
        for x in range(n_g):
            if typed and typed[1][x] != typed[0][u]:
                continue
            for y in range(n_g):
                if y == x or y in G[x]:
                    continue
                key = (min(x, y), max(x, y))
                if key in found:
                    continue
                if typed and typed[1][y] != typed[0][v]:
                    continue
                f = {u: x, v: y}
                used = {x, y}
                if _extend(H, G, order, back, 2, f, used, n_g, typed):
                    found.add(key)
    return found


def _extend(H, G, order, back, i, f, used, n_g, typed):
    if i == len(order):
        return True
    w = order[i]
    nb = back[i]
    if nb:
        cands = set(G[f[nb[0]]])
        for x in nb[1:]:
            cands &= G[f[x]]
    else:
        cands = range(n_g)
    for c in cands:
        if c in used:
            continue
        if typed and typed[1][c] != typed[0][w]:
            continue
        f[w] = c; used.add(c)
        if _extend(H, G, order, back, i + 1, f, used, n_g, typed):
            return True
        del f[w]; used.discard(c)
    return False


def run(H, n_h, G, n_g, typed=None, maxrounds=10**6, trace=False):
    G = [set(s) for s in G]
    t = 0
    hist = []
    while t < maxrounds:
        add = new_edges(H, n_h, G, n_g, typed)
        if not add:
            return t, G, hist
        for a, b in add:
            G[a].add(b); G[b].add(a)
        if trace:
            hist.append(add)
        t += 1
    raise RuntimeError("no stabilization")


def theta(a, b, c):
    """Theta graph with poles 0,1 and arms of lengths a,b,c."""
    edges = []
    n = 2
    for L in (a, b, c):
        prev = 0
        for _ in range(L - 1):
            edges.append((prev, n)); prev = n; n += 1
        edges.append((prev, 1))
    return n, edges
