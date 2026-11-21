# Independent reconstruction of the ordinary 166,464-vertex rule from the Lean definitions.
from collections import deque
roles=[('x','e0'),('x','e1'),('y','e0'),('y','e1')]+[(s,k) for s in 'xy' for k in ['u0','v0','u1','v1']]
kindId={('x','e'):0,('y','e'):1,('x','u0'):2,('x','v0'):3,('x','u1'):4,('x','v1'):5,('y','u0'):6,('y','v0'):7,('y','u1'):8,('y','v1'):9}
def kind(r): return (r[0], 'e' if r[1][0]=='e' else r[1])
code={r:kindId[kind(r)] for r in roles}
RE=[]
for s in 'xy':
    for k in '01':
        RE += [((s,'e0'),(s,'u'+k)),((s,'u'+k),(s,'v'+k)),((s,'v'+k),(s,'e1'))]
RE += [(('x','e0'),('y','e0')),(('x','e1'),('y','e1'))]
radj={r:set() for r in roles}
for a,b in RE: radj[a].add(b); radj[b].add(a)
order={('x','e0'):0,('x','e1'):1,('y','e0'):2,('y','e1'):3,('x','u0'):4,('x','v0'):5,('x','u1'):6,('x','v1'):7,('y','u0'):8,('y','v0'):9,('y','u1'):10,('y','v1'):11}
def selected(r): return min(radj[r], key=lambda s: order[s])
def port(r,s): return 0 if s==selected(r) else 1
V={}; E=set()
def vid(x):
    if x not in V: V[x]=len(V)
    return V[x]
def add(a,b):
    a,b=vid(a),vid(b); assert a!=b; E.add((min(a,b),max(a,b)))
for r in roles:
    for blk,s in (('P',20+code[r]),('S',50-code[r])):
        for a in range(s): vid((r,blk,'c',a))
        cedges=[(0,1)]+[(0,j) for j in range(2,s)]+[(1,j) for j in range(2,s)]
        assert len(cedges)==2*s-3
        pg=0
        for (a,d) in cedges:
            add((r,blk,'c',a),(r,blk,'c',d))
            for _ in range(103):
                add((r,blk,'p',pg),(r,blk,'c',a)); add((r,blk,'p',pg),(r,blk,'c',d)); pg+=1
    for k in (0,1): add((r,'P','c',k),(r,'S','c',k))
for a,b in RE: add((a,'P','c',port(a,b)),(b,'P','c',port(b,a)))
n=len(V); print('vertices',n,'edges',len(E))
adj=[set() for _ in range(n)]
for a,b in E: adj[a].add(b); adj[b].add(a)
print('min degree',min(len(s) for s in adj),'max degree',max(len(s) for s in adj))
# 2-connectivity via articulation points (iterative Tarjan)
import sys
disc=[-1]*n; low=[0]*n; t=0; art=set()
for root in [0]:
    stack=[(root,-1,iter(adj[root]))]; disc[root]=low[root]=t; t+=1; rootch=0
    while stack:
        v,p,it=stack[-1]
        for w in it:
            if disc[w]==-1:
                disc[w]=low[w]=t; t+=1; stack.append((w,v,iter(adj[w])))
                if v==root: rootch+=1
                break
            elif w!=p: low[v]=min(low[v],disc[w])
        else:
            stack.pop()
            if stack:
                u=stack[-1][0]; low[u]=min(low[u],low[v])
                if u!=root and low[v]>=disc[u]: art.add(u)
    if rootch>1: art.add(root)
print('connected',all(d>=0 for d in disc),'articulation points',len(art))
# treewidth<=2 via series-parallel reduction
A=[set(s) for s in adj]; alive=n
q=deque(v for v in range(n) if len(A[v])<=2); dead=[False]*n
while q:
    v=q.popleft()
    if dead[v] or len(A[v])>2: continue
    nb=list(A[v]); dead[v]=True; alive-=1
    for w in nb: A[w].discard(v)
    if len(nb)==2:
        a,b=nb; A[a].add(b); A[b].add(a)
    for w in nb:
        if len(A[w])<=2: q.append(w)
print('remaining after SP reduction',alive,'=> treewidth<=2' if alive==0 else '=> treewidth>=3')
# triangle exists?
tri=any(len(adj[a]&adj[b])>0 for a,b in E)
print('has triangle (tw>=2)',tri)
