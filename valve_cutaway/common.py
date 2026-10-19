import numpy as np, pickle, pyvista as pv
DISC={47,48,49,50,51}
STRUT=set(range(4,29))|{38,39,40,41,42,43,44}
PH={16:16.67,33:33.33,50:50.0,66:66.67,83:83.33}
def facemesh(fa, ids):
    Vs=[];Fs=[];o=0
    for t in ids:
        v=fa[t]; Vs.append(v['V']); Fs.append(np.c_[np.full(len(v['F']),3),v['F']+o]); o+=len(v['V'])
    return pv.PolyData(np.vstack(Vs), np.vstack(Fs).ravel()).clean()
def geom(k):
    fa,ed=pickle.load(open(f"geom{k}.pkl","rb"))
    return fa,ed
def edges_mesh(ed, sel=None):
    pts=[];segs=[];o=0
    for e in ed:
        if sel is not None and not (set(e['faces'])&sel): continue
        pts.append(e['V']); segs.append(np.c_[np.full(len(e['S']),2),e['S']+o]); o+=len(e['V'])
    return pv.PolyData(np.vstack(pts), lines=np.vstack(segs).ravel())
def lines(k):
    d=np.load(f"ph{k}.npz"); return d['X'],d['tau'],d['U'],d['off']
def plane_normal(fa,t):
    V=fa[t]['V']; c=V.mean(0); _,_,vt=np.linalg.svd(V-c); return c,vt[2]
