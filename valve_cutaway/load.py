import numpy as np, sys, re
from pathlib import Path
def load(f):
    txt=Path(f).read_text()
    a=txt.index('[Data]'); b=txt.index('[Lines]')
    data=txt[a:b].splitlines()[2:]
    pts=np.loadtxt([l for l in data if l.strip()],delimiter=',')
    lines=np.loadtxt([l for l in txt[b:].splitlines()[1:] if l.strip()],delimiter=',',dtype=np.int64)
    return pts[:,1:4], pts[:,4], pts[:,5], lines
def polylines(lines):
    # chain segments into polylines
    nxt={}; has_prev=set()
    for a,b in lines: nxt[a]=b; has_prev.add(b)
    starts=[a for a in nxt if a not in has_prev]
    out=[]
    for s in starts:
        p=[s]
        while p[-1] in nxt: p.append(nxt[p[-1]])
        out.append(np.array(p))
    return out
if __name__=="__main__":
    import os
    base=Path(os.environ.get("CSV_DIR", Path(__file__).resolve().parent.parent/"35cc_inputs"/"120BPM"))
    for f,ph in [("35cc 120BPM 16.67% Updated.csv",16.67),("35cc 120BPM 33.33%.csv",33.33),("35cc 120BPM 50.00%.csv",50.0),("35cc 120BPM 66.67%.csv",66.67),("35cc 120BPM 83.33%.csv",83.33)]:
        X,tau,U,L=load(base/f); pl=polylines(L)
        off=np.cumsum([0]+[len(p) for p in pl]); idx=np.concatenate(pl)
        np.savez_compressed(f"ph{int(ph)}.npz",X=X[idx],tau=tau[idx],U=U[idx],off=off)
        print(ph,len(X),len(pl),np.median([len(p) for p in pl]), tau.max(), U.max())
