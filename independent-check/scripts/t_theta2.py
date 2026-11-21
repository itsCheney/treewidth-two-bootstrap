import itertools, subprocess, sys
sys.path.insert(0,'.')
from hbp import theta
bad=cnt=0
for a,b,c in itertools.combinations_with_replacement(range(2,12),3):
    n,E=theta(a,b,c)
    if n>18: continue
    inp=f"{n} {len(E)}\n"+"".join(f"{x} {y}\n" for x,y in E)
    out=subprocess.run(['./sim'],input=inp+inp,capture_output=True,text=True,timeout=3600).stdout.split()
    t=int(out[0]); G=[int(r) for r in out[1:]]
    col={0:0}; st=[0]; A=[set() for _ in range(n)]
    for x,y in E: A[x].add(y); A[y].add(x)
    while st:
        x=st.pop()
        for y in A[x]:
            if y not in col: col[y]=1-col[x]; st.append(y)
    full=all((G[x]>>y)&1 for x in range(n) for y in range(n) if x!=y)
    bip=all(((G[x]>>y)&1)==(col[x]!=col[y]) for x in range(n) for y in range(n) if x!=y)
    if len({a%2,b%2,c%2})==2: ok=full; pred='clique'
    elif a==b==c: ok=(t==0); pred='stable'
    else: ok=bip; pred='bipartite'
    cnt+=1; bad+=not ok
    print(a,b,c,'n',n,'rounds',t,pred,'OK' if ok else 'FAIL',flush=True)
print('cases',cnt,'failures',bad)
