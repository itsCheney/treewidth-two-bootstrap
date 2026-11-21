import itertools, sys
from hbp import *
def colour(n, A):
    col={0:0}; st=[0]
    while st:
        x=st.pop()
        for y in A[x]:
            if y not in col: col[y]=1-col[x]; st.append(y)
            elif col[y]==col[x]: return None
    return col
bad=0; cnt=0
for a,b,c in itertools.combinations_with_replacement(range(2,8),3):
    n,E=theta(a,b,c)
    if n>14: continue
    H=adj_from_edges(n,E)
    t,G,_=run(H,n,H,n)
    m=sum(len(s) for s in G)//2
    col=colour(n,H)
    if len({a%2,b%2,c%2})==2: pred='clique'; ok = m==n*(n-1)//2
    elif a==b==c: pred='stable'; ok = (t==0)
    else:
        pred='bipartite'
        ok = all((y in G[x])==(col[x]!=col[y]) for x in range(n) for y in range(n) if x!=y)
    cnt+=1; bad+= not ok
    print(a,b,c,'n=',n,'rounds=',t,'edges',m,pred,'OK' if ok else 'FAIL'); sys.stdout.flush()
print('cases',cnt,'failures',bad)
