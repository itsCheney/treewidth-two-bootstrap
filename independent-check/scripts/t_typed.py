import random, math, sys
from hbp import *
KINDS=['e','u0','u1','v0','v1']
# rule
R={}; 
for s in 'xy':
    for r in ['e0','e1','u0','u1','v0','v1']: R[(s,r)]=len(R)
RE=[]
for s in 'xy':
    for k in '01':
        RE += [(R[(s,'e0')],R[(s,'u'+k)]),(R[(s,'u'+k)],R[(s,'v'+k)]),(R[(s,'v'+k)],R[(s,'e1')])]
RE += [(R[('x','e0')],R[('y','e0')]),(R[('x','e1')],R[('y','e1')])]
assert len(RE)==14
Hn=12; H=adj_from_edges(Hn,RE)
htype=[None]*12
for (s,r),i in R.items(): htype[i]=(s, 'e' if r[0]=='e' else r)
def host(p,q,E):
    V={}
    for s,m in (('x',p),('y',q)):
        for i in range(m):
            for k in KINDS: V[(s,k,i)]=len(V)
    ed=[]
    for s,m in (('x',p),('y',q)):
        for i in range(m):
            for k in '01':
                ed += [(V[(s,'e',i)],V[(s,'u'+k,i)]),(V[(s,'u'+k,i)],V[(s,'v'+k,i)]),(V[(s,'v'+k,i)],V[(s,'e',(i+1)%m)])]
    ed += [(V[('x','e',i)],V[('y','e',j)]) for (i,j) in E]
    gt=[None]*len(V)
    for (s,k,i),v in V.items(): gt[v]=(s,k)
    return V, adj_from_edges(len(V),ed), gt
def plus(E,p,q):
    return E | {((i+1)%p,(j+1)%q) for i,j in E} | {((i-1)%p,(j-1)%q) for i,j in E}
random.seed(1); fails=0; tests=0
for p in range(3,7):
  for q in range(3,7):
    for trial in range(4):
        E={(i,j) for i in range(p) for j in range(q) if random.random()<0.15} or {(0,0)}
        V,G,gt=host(p,q,E)
        inv={v:k for k,v in V.items()}
        add=new_edges(H,Hn,G,len(V),typed=(htype,gt))
        got=set()
        ok=True
        for a,b in add:
            A,B=inv[a],inv[b]
            if A[0]=='y': A,B=B,A
            if not(A[0]=='x' and B[0]=='y' and A[1]=='e' and B[1]=='e'): ok=False; break
            got.add((A[2],B[2]))
        ok = ok and (E|got)==plus(E,p,q)
        tests+=1; fails+= not ok
        if not ok: print('FAIL',p,q,sorted(E),sorted(got))
print('typed one-round tests',tests,'fails',fails); sys.stdout.flush()
# full run from single edge, typed
for p,q in [(3,3),(3,4),(4,5),(3,5),(5,6),(4,6),(6,7)]:
    V,G,gt=host(p,q,{(0,0)})
    t,_,_=run(H,Hn,G,len(V),typed=(htype,gt))
    print('typed run',p,q,'rounds',t,'lcm/2',math.lcm(p,q)//2, 'OK' if t==math.lcm(p,q)//2 else 'FAIL'); sys.stdout.flush()
# untyped one round, does the bare 12-vertex rule also give E+?
for p,q in [(3,4),(4,5),(5,6)]:
    E={(0,0)}
    V,G,gt=host(p,q,E)
    add=new_edges(H,Hn,G,len(V))
    inv={v:k for k,v in V.items()}
    print('untyped',p,q,'new edges',len(add),'sample',[ (inv[a],inv[b]) for a,b in list(add)[:4]])
