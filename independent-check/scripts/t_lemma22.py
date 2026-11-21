# Random test of Lemma II.2: W full copy of H; P (|P|=h) all adjacent to f(x), Q all adjacent to f(y);
# random extra host edges allowed. Claim: after 2h rounds all P-Q pairs are present.
import random, subprocess, itertools
from hbp import theta
def rules():
    yield 'theta(2,2,3)', *theta(2,2,3)
    yield 'theta(2,3,3)', *theta(2,3,3)
    yield 'theta(2,2,4)', *theta(2,2,4)
    yield 'C5', 5, [(i,(i+1)%5) for i in range(5)]
    yield 'C4', 4, [(i,(i+1)%4) for i in range(4)]
    # K_{2,3} with a pendant two-edge path from a pole (leaf is a base vertex; bridge)
    n,E=theta(2,2,2); yield 'K23+pendant P2', n+2, E+[(0,n),(n,n+1)]
    # two triangles-ish: bowtie subdivided: centre 0, two returning routes of length 3 (triangles)
    yield 'bowtie', 5, [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
    # short route plus long route between different base vertices: K_{1,3}-like with subdivisions
    yield 'theta(2,2,5)', *theta(2,2,5)
random.seed(11)
tot=bad=0
for name,h,E in rules():
    A=[set() for _ in range(h)]
    for a,b in E: A[a].add(b); A[b].add(a)
    worst=0
    for (x,y) in E:
        for orient in (0,1):
            if orient: x,y=y,x
            for trial in range(6):
                n=3*h
                if n>64: continue
                P=list(range(h,2*h)); Q=list(range(2*h,3*h))
                edges=set((min(a,b),max(a,b)) for a,b in E)
                edges|={(x,p) for p in P}|{(y,q) for q in Q}
                dens=random.choice([0,0.03,0.08])
                for a,b in itertools.combinations(range(n),2):
                    if random.random()<dens: edges.add((a,b))
                # Lemma assumes nothing about P-Q edges; keep random ones too
                inp=f"{h} {len(E)}\n"+"".join(f"{a} {b}\n" for a,b in E)+f"{n} {len(edges)}\n"+"".join(f"{a} {b}\n" for a,b in edges)
                out=subprocess.run(['./simr',str(2*h)],input=inp,capture_output=True,text=True).stdout.split()
                G=[int(r) for r in out[1:]]
                ok=all((G[p]>>q)&1 for p in P for q in Q)
                # also find the first round at which P-Q is complete
                tot+=1; bad+=not ok
                if not ok: print('FAIL',name,x,y,dens)
    print(name,'h=',h,'tests so far',tot,'failures',bad,flush=True)
print('TOTAL',tot,'failures',bad)
