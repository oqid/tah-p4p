import sys, json, pickle; sys.path.insert(0,'.')
from render_v3 import *
from features import feats3d, largest
import vtk
C_BF='#2E6DB4'
def project(p,P):
    out=[]
    for pt in P:
        c=vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(*map(float,pt))
        x,y=c.GetComputedDoubleDisplayValue(p.renderer); out.append((x,p.window_size[1]-y))
    return out
def render_v4(k,tag,camera,n=800,window=(1400,1400),seed=1,amin=-14e-3,ghost=True):
    c=classify(k); off=c['off']; clip=clipper(amin,None)
    idx=np.random.default_rng(seed).choice(len(off)-1,size=min(n,len(off)-1),replace=False)
    pl=clip(polys(c['X'],off,idx,np.log10(np.clip(c['tau'],1e-2,1e2))))
    p=base_scene(k,clip,window)
    p.add_mesh(pl.tube(radius=0.00008,n_sides=8),scalars='s',cmap='magma',clim=(-2,2),show_scalar_bar=False,smooth_shading=True,ambient=0.3)
    ghost_mesh=None
    if ghost:
        g=pv.read(f"vol{k}.vtk"); rv=g.contour([0.03],scalars='rev')
        rv=rv.smooth_taubin(n_iter=40) if rv.n_points else rv
        rv=clip(rv)
        if rv.n_points:
            rv=rv.connectivity('largest') if rv.n_points else rv
            rv=pv.PolyData(rv.extract_surface(algorithm='dataset_surface')) if not isinstance(rv,pv.PolyData) else rv
            p.add_mesh(rv,color=C_BF,opacity=0.16,smooth_shading=True,ambient=0.6)
            ghost_mesh=rv
    p.camera_position=camera; p.enable_anti_aliasing('ssaa')
    fn=f"v4_{tag}_{k}.png"; p.screenshot(fn)
    F=feats3d(k)
    F['recirc']=F.get('recirc',{}) ; 
    proj={kk:project(p,[v['p']])[0] for kk,v in F.items() if 'p' in v}
    # also an in-image unit length for the recirculation glyph: project centre +/- 2 mm along A
    if 'p' in F['recirc']:
        a0,a1=project(p,[F['recirc']['p']-0.002*A,F['recirc']['p']+0.002*A]); proj['recirc_scale']=float(np.hypot(a1[0]-a0[0],a1[1]-a0[1]))
        e0,e1=project(p,[F['recirc']['p'],F['recirc']['p']+0.002*E2]); proj['e2_dir']=(e1[0]-e0[0],e1[1]-e0[1])
        a0,a1=project(p,[F['recirc']['p'],F['recirc']['p']+0.002*A]); proj['a_dir']=(a1[0]-a0[0],a1[1]-a0[1])
    if ghost_mesh is not None:
        sil=vtk.vtkPolyDataSilhouette(); sil.SetInputData(ghost_mesh); sil.SetCamera(p.camera); sil.SetEnableFeatureAngle(False); sil.BorderEdgesOn(); sil.Update()
        sl=pv.wrap(sil.GetOutput()).strip()
        lines=[]
        for i in range(sl.n_cells):
            pts=np.asarray(sl.get_cell(i).points)
            if len(pts)>3: lines.append(project(p,pts))
        proj['backflow_outline']=lines
    json.dump(proj,open(f"v4_{tag}_{k}.json","w"))
    p.close()
if __name__=="__main__":
    k=int(sys.argv[1]); tag=sys.argv[2]
    camera=view('front_valve')[0] if tag=='front' else swivel(*[float(x) for x in tag.replace('sw','').split('e')])
    render_v4(k,tag,camera)
