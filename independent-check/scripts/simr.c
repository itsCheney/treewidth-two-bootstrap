// Bitset H-bootstrap simulator (n<=64): new edge xy iff some edge uv of H and injective f with
// f(u)=x,f(v)=y and all other H-edges mapped into the current graph G.
#include <stdio.h>
#include <stdint.h>
#include <string.h>
typedef uint64_t u64;
int hn, gn; u64 H[64], G[64], NEW[64];
int ord[64], back[64][64], nback[64];
int u0,v0; int f[64]; u64 used;
static int ext(int i){
  if(i==hn) return 1;
  int w=ord[i]; u64 c = ~used & ((gn==64)?~0ULL:((1ULL<<gn)-1));
  for(int k=0;k<nback[i];k++) c &= G[f[back[i][k]]];
  while(c){ int z=__builtin_ctzll(c); c&=c-1; f[w]=z; used|=1ULL<<z; if(ext(i+1)) return 1; used&=~(1ULL<<z);} return 0;
}
int round_(void){
  memset(NEW,0,sizeof NEW); int any=0;
  for(int u=0;u<hn;u++) for(int v=0;v<hn;v++) if((H[u]>>v&1)){
    // order: u, v, then greedy
    int pos[64]; for(int i=0;i<hn;i++) pos[i]=-1; ord[0]=u; ord[1]=v; pos[u]=0; pos[v]=1;
    for(int i=2;i<hn;i++){ int best=-1,bs=-1; for(int w=0;w<hn;w++) if(pos[w]<0){ int s=0; for(int x=0;x<hn;x++) if(pos[x]>=0 && (H[w]>>x&1)) s++; if(s>bs){bs=s;best=w;} } ord[i]=best; pos[best]=i; }
    for(int i=0;i<hn;i++){ nback[i]=0; int w=ord[i]; for(int x=0;x<hn;x++) if((H[w]>>x&1) && pos[x]<i && !((w==u&&x==v)||(w==v&&x==u))) back[i][nback[i]++]=x; }
    for(int x=0;x<gn;x++) for(int y=0;y<gn;y++){ if(x==y || (G[x]>>y&1) || (NEW[x]>>y&1)) continue;
      f[u]=x; f[v]=y; used=(1ULL<<x)|(1ULL<<y);
      if(ext(2)){ NEW[x]|=1ULL<<y; NEW[y]|=1ULL<<x; any=1; } }
  }
  for(int x=0;x<gn;x++) G[x]|=NEW[x]; return any;
}
#include <stdlib.h>
int main(int argc,char**argv){
  // read: hn m edges ; then gn m edges ; output rounds and final adjacency count / matrix rows
  int m; if(scanf("%d %d",&hn,&m)!=2) return 1; for(int i=0;i<m;i++){int a,b;scanf("%d %d",&a,&b);H[a]|=1ULL<<b;H[b]|=1ULL<<a;}
  scanf("%d %d",&gn,&m); for(int i=0;i<m;i++){int a,b;scanf("%d %d",&a,&b);G[a]|=1ULL<<b;G[b]|=1ULL<<a;}
  int lim=argc>1?atoi(argv[1]):1000000; int t=0; while(t<lim && round_()) t++;
  printf("%d\n",t); for(int x=0;x<gn;x++) printf("%llu\n",(unsigned long long)G[x]); return 0;
}
