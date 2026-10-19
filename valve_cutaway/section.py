import sys; sys.path.insert(0,'.')
from common import *
A=np.array([0.427,0.82,0.382]); A/=np.linalg.norm(A)
H=np.array([0.0,-0.423,0.906]); H-= (H@A)*A; H/=np.linalg.norm(H)
E2=np.cross(H,A)
O=np.array([8.7,11.1,4.7])/1000
def to2d(P): Q=P-O; return np.c_[Q@A, Q@E2]*1000, (Q@H)*1000
def section(fa, ids):
    m=facemesh(fa,ids); s=m.slice(normal=H, origin=O)
    if s.n_points==0: return []
    s=s.strip()  # polylines
    out=[]
    for i in range(s.n_cells):
        c=s.get_cell(i); out.append(to2d(np.asarray(c.points))[0])
    return out
def tangents(X,off):
    T=np.zeros_like(X)
    for s,e in zip(off[:-1],off[1:]):
        x=X[s:e]
        if len(x)<2: continue
        g=np.gradient(x,axis=0); T[s:e]=g
    n=np.linalg.norm(T,axis=1,keepdims=True); return T/np.maximum(n,1e-12)
