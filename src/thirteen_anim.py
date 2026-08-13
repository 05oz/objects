"""One object, thirteen times.

Each of the thirteen rigid {I3,TT4}-free witnesses of k(3,4) = 21 is brought into its
closest found agreement with the published witness W by relabelling (simulated annealing
over the 20! relabellings; 82-98% of arcs match). Vertices are then fixed at W's spectral
embedding and only the arcs move.

What stays is what the family shares. What blinks is what makes them thirteen objects and
not one -- the residual no relabelling could reconcile, which is exactly the obstruction to
isomorphism that rigidity guarantees cannot be smoothed away.

Hue encodes persistence: how many of the thirteen an arc survives in.
"""
import json, glob, math, os
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__))
SRC="/Users/kirt/Documents/reserch math/certify-repo/k34add-certificates"
S=1500; R=S*0.40; C=S/2
W=[json.load(open(f)) for f in sorted(glob.glob(f"{SRC}/w*.json"))]
P=json.load(open("/private/tmp/claude-501/-Applications/20ab6041-f5e3-44f3-9953-5c9b386f5721/scratchpad/perms.json"))["perms"]
N=W[0]["N"]

def embed(arcs):
    A=np.zeros((N,N))
    for u,v in arcs: A[u][v]=1.0
    K=A-A.T
    w,V=np.linalg.eig(K)
    v=V[:,np.argsort(-np.abs(w.imag))[0]]
    x,y=v.real.copy(),v.imag.copy()
    x-=x.mean(); y-=y.mean()
    s=max(np.abs(x).max(),np.abs(y).max()) or 1
    return x/s,y/s

x,y=embed([tuple(a) for a in W[0]["arcs"]])
pos=[(C+x[i]*R, C+y[i]*R) for i in range(N)]

frames=[]
for k in range(13):
    p=P[k]
    inv=[0]*N
    for i,pi in enumerate(p): inv[pi]=i     # the annealer scores B's (p[u],p[v]) against A's (u,v)
    frames.append(set((inv[u],inv[v]) for u,v in map(tuple,W[k]["arcs"])))
union=sorted(set().union(*frames))
persist={a:sum(a in f for f in frames) for a in union}

def hsv(h,s,v):
    h=(h%360)/60.0; i=int(h)%6; f=h-int(h)
    p,q,t=v*(1-s),v*(1-s*f),v*(1-s*(1-f))
    r,g,b=((v,t,p),(q,v,p),(p,v,t),(p,q,v),(t,p,v),(v,p,q))[i]
    return "#%02x%02x%02x"%(int(r*255),int(g*255),int(b*255))

DUR=6.5                       # 0.5 s per witness
out=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{S}" height="{S}" viewBox="0 0 {S} {S}">',
     f'<rect width="{S}" height="{S}" fill="#0a0d12"/>',
     '<g fill="none" stroke-linecap="round">']
for a in union:
    u,v=a; n=persist[a]
    x1,y1=pos[u]; x2,y2=pos[v]
    dx,dy=x2-x1,y2-y1; L=math.hypot(dx,dy) or 1
    mx,my=(x1+x2)/2,(y1+y2)/2
    qx,qy = mx-dy/L*L*0.14, my+dx/L*L*0.14
    t=(n-1)/12.0                                  # 0 = appears once, 1 = in all thirteen
    col=hsv(198-166*t, 0.40+0.26*t, 0.72+0.20*t)
    on=[("0.80" if a in f else "0.07") for f in frames]
    # hold each state for 72% of its slot, then flick to the next: alive, not strobing
    vals=";".join(v for x in on for v in (x,x))+";"+on[0]
    kt=";".join(f"{k:.4f}" for i in range(13) for k in (i/13.0,(i+0.72)/13.0))+";1"
    wid=1.1+2.4*t
    base="0.80" if a in frames[0] else "0.07"
    out.append(f'<path d="M{x1:.1f} {y1:.1f}Q{qx:.1f} {qy:.1f} {x2:.1f} {y2:.1f}" stroke="{col}" '
               f'stroke-width="{wid:.2f}" opacity="0.045">'
               f'<animate attributeName="opacity" values="{vals}" keyTimes="{kt}" '
               f'dur="{DUR}s" calcMode="linear" repeatCount="indefinite"/></path>')
out.append('</g><g>')
for i,(px,py) in enumerate(pos):
    out.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="6.0" fill="#f2ece0" opacity="0.97"/>')
out.append('</g></svg>')
open(os.path.join(HERE,"thirteen.svg"),"w").write("".join(out))
hist={}
for a,n in persist.items(): hist[n]=hist.get(n,0)+1
print("union arcs:",len(union))
print("persistence histogram (arcs present in n of 13):",dict(sorted(hist.items())))
print("in all 13:",hist.get(13,0),"| in 1 only:",hist.get(1,0))
