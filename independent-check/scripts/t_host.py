# Full-scale independent checks of Lemmas III.4/III.5 premises and of the lifting direction,
# on the real 166,464-vertex rule H and the real ordinary host G(E) (built from the preprint text).
import random, itertools, sys
exec(open('t_rule.py').read().split("n=len(V); print")[0])   # builds rule: V (name->id), E (edges), roles, RE, code, port
Hnames={i:nm for nm,i in V.items()}; Hedges=E; Hcode=code; Hport=port
Hadj=[set() for _ in range(len(V))]
for a,b in Hedges: Hadj[a].add(b); Hadj[b].add(a)
KINDS=['e','u0','v0','u1','v1']
kindId={('x','e'):0,('y','e'):1,('x','u0'):2,('x','v0'):3,('x','u1'):4,('x','v1'):5,('y','u0'):6,('y','v0'):7,('y','u1'):8,('y','v1'):9}
def size(c,B): return 20+c if B=='P' else 50-c
class Host:
    def __init__(s,p,q,Ecross):
        s.p,s.q,s.E=p,q,set(Ecross)
        s.roles=[(sd,k,i) for sd,m in (('x',p),('y',q)) for i in range(m) for k in KINDS]
        s.code={h:kindId[(h[0],h[1])] for h in s.roles}
        s.F={h:set() for h in s.roles}
        for sd,m in (('x',p),('y',q)):
            for i in range(m):
                for k in '01':
                    for a,b in [((sd,'e',i),(sd,'u'+k,i)),((sd,'u'+k,i),(sd,'v'+k,i)),((sd,'v'+k,i),(sd,'e',(i+1)%m))]:
                        s.F[a].add(b); s.F[b].add(a)
        for i,j in s.E: s.F[('x','e',i)].add(('y','e',j)); s.F[('y','e',j)].add(('x','e',i))
    def desig(s,X,Y):  # designated clique adjacency, X=(h,B)
        (h,B),(g,C)=X,Y
        if h==g: return B!=C
        return B=='P' and C=='P' and g in s.F[h]
    def adj(s,x,y):
        if x==y: return False
        if x[2]=='p' and y[2]=='p': return False
        if x[2]=='p': x,y=y,x
        if y[2]=='p': return (x[0],x[1])==(y[0],y[1])
        return (x[0],x[1])==(y[0],y[1]) or s.desig((x[0],x[1]),(y[0],y[1]))
    def nbrs(s,x):
        h,B=x[0],x[1]; out=set()
        blocks=[(h,B)]
        if x[2]=='c':
            blocks += [(h,'S' if B=='P' else 'P')]
            if B=='P': blocks += [(g,'P') for g in s.F[h]]
            for (g,C) in blocks:
                for a in range(size(s.code[g],C)): out.add((g,C,'c',a))
            sz=size(s.code[h],B)
            for k in range(103*(2*sz-3)): out.add((h,B,'p',k))
        else:
            for a in range(size(s.code[h],B)): out.add((h,B,'c',a))
        out.discard(x); return out
    def order(s): return sum(size(s.code[h],B)+103*(2*size(s.code[h],B)-3) for h in s.roles for B in 'PS')

# ---- (1) rule side: every core edge of H has >=103 common neighbours, each page serves one core edge
core_edges=[(a,b) for a,b in Hedges if Hnames[a][2]=='c' and Hnames[b][2]=='c' and Hnames[a][:2]==Hnames[b][:2]]
mincn=min(len(Hadj[a]&Hadj[b]) for a,b in core_edges)
pages=[v for v,nm in Hnames.items() if nm[2]=='p']
print('rule: core edges',len(core_edges),'min common nbrs',mincn,
      'every page has degree 2 on a core edge:',all(len(Hadj[v])==2 and Hnames[min(Hadj[v])][:2]==Hnames[max(Hadj[v])][:2] for v in pages))
print('rule: book core minus any one core edge stays connected (s>=3 star-fan structure) -- checked below')
def core_conn_after_delete(s):
    ce=[(0,1)]+[(0,j) for j in range(2,s)]+[(1,j) for j in range(2,s)]
    for d in ce:
        A={i:set() for i in range(s)}
        for a,b in ce:
            if (a,b)!=d: A[a].add(b); A[b].add(a)
        seen={0}; st=[0]
        while st:
            x=st.pop()
            for y in A[x]:
                if y not in seen: seen.add(y); st.append(y)
        if len(seen)!=s: return False
    return True
print('  ', all(core_conn_after_delete(s) for s in list(range(20,30))+list(range(41,51))))

# ---- (2) host side: common-neighbour bounds by edge class, for initial and for dense E
for (p,q,Ec,label) in [(3,3,{(0,0)},'initial p=q=3'),(3,4,{(i,j) for i in range(3) for j in range(4)},'complete E p=3,q=4'),(4,5,{(0,0),(1,1),(3,2)},'sparse E p=4,q=5')]:
    G=Host(p,q,Ec)
    worst={'private':0,'cross-clique':0,'same-clique-min':10**9}
    for h in G.roles:
        for B in 'PS':
            x=(h,B,'c',0); Nx=G.nbrs(x)
            for y in [(h,B,'c',1),(h,B,'p',0)]+[(h,'S' if B=='P' else 'P','c',0)]+([(g,'P','c',0) for g in G.F[h]] if B=='P' else []):
                assert G.adj(x,y)
                cn=len(Nx & G.nbrs(y))
                if y[2]=='p': worst['private']=max(worst['private'],cn)
                elif (y[0],y[1])!=(h,B): worst['cross-clique']=max(worst['cross-clique'],cn)
                else: worst['same-clique-min']=min(worst['same-clique-min'],cn)
    # triangle-freeness of F(E)
    tf=not any(G.F[a]&G.F[b] for a in G.roles for b in G.F[a])
    print(label,'| host order',G.order(),'=69360(p+q)?',G.order()==69360*(p+q),'| F(E) triangle-free',tf,'|',worst)
    sys.stdout.flush()

# ---- (3) max degree and connectivity of G_{p,q}
for p,q in [(3,3),(3,4),(5,7)]:
    G=Host(p,q,{(0,0)})
    md=0
    for h in G.roles:
        for B in 'PS':
            md=max(md,len(G.nbrs((h,B,'c',0))),len(G.nbrs((h,B,'p',0))))
    print('G_{%d,%d} max degree'%(p,q),md,'<=14016',md<=14016)

# ---- (4) explicit lifted witnesses, checked edge by edge against the host oracle
def typed_witness(G,target):
    """Brute-force a type-preserving injective map of the 12-vertex rule L minus target edge into F(E)."""
    Lr=roles; Ladj={r:set() for r in Lr}
    for a,b in RE: Ladj[a].add(b); Ladj[b].add(a)
    def kind(r): return (r[0],'e' if r[1][0]=='e' else r[1])
    sols=[]
    order=[('x','e0'),('x','u0'),('x','v0'),('x','e1'),('x','u1'),('x','v1'),('y','e0'),('y','u0'),('y','v0'),('y','e1'),('y','u1'),('y','v1')]
    def ok(r,h,phi):
        if kindId[(h[0],h[1])]!=kindId[kind(r)] or h in phi.values(): return False
        for s in Ladj[r]:
            if s in phi and {r,s}!=set(target) and phi[s] not in G.F[h]: return False
        return True
    def rec(i,phi):
        if i==12: sols.append(dict(phi)); return
        r=order[i]
        for h in G.roles:
            if ok(r,h,phi):
                phi[r]=h; rec(i+1,phi); del phi[r]
    rec(0,{}); return sols
random.seed(7)
checked=0
for p,q,Ec in [(3,3,{(0,0)}),(4,5,{(0,0),(2,3)}),(5,6,{(1,4)})]:
    G=Host(p,q,Ec)
    for tgt in [(('x','e0'),('y','e0')),(('x','e1'),('y','e1'))]:
        for phi in typed_witness(G,tgt):
            a,b=tgt
            if phi[b] in G.F[phi[a]]: continue   # pair already present
            # prescribed endpoints: random physical points in the two primary cliques
            pa=random.randrange(size(G.code[phi[a]],'P')); pb=random.randrange(size(G.code[phi[b]],'P'))
            mp={}
            for r in roles:
                h=phi[r]
                for B in 'PS':
                    s=size(Hcode[r],B); perm=list(range(s)); random.shuffle(perm)
                    if B=='P' and r in tgt:
                        prt=Hport(a,b) if r==a else Hport(b,a); want=pa if r==a else pb
                        j=perm.index(want); perm[j],perm[prt]=perm[prt],perm[j]
                    for c in range(s): mp[V[(r,B,'c',c)]]=(h,B,'c',perm[c])
                    pg=list(range(103*(2*s-3))); random.shuffle(pg)
                    for k in range(len(pg)): mp[V[(r,B,'p',k)]]=(h,B,'p',pg[k])
            assert len(mp)==len(V) and len(set(mp.values()))==len(V), 'not injective'
            ta=V[(a,'P','c',Hport(a,b))]; tb=V[(b,'P','c',Hport(b,a))]
            assert mp[ta]==(phi[a],'P','c',pa) and mp[tb]==(phi[b],'P','c',pb)
            bad=[(u,v) for u,v in Hedges if {u,v}!={ta,tb} and not G.adj(mp[u],mp[v])]
            assert not bad, ('missing',len(bad))
            assert not G.adj(mp[ta],mp[tb])
            checked+=1
    print('p,q',p,q,'lifted witnesses verified so far',checked); sys.stdout.flush()
print('ALL lifted witnesses verified:',checked)
