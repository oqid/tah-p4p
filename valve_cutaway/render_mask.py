"""Second pass from the same camera: a black-on-white mask of the backflow volume, for a multiply-blend shade."""
import sys; sys.path.insert(0,'.')
from render_v4 import *
def render_mask(k,tag,camera,window=(1400,1400),amin=-14e-3):
    clip=clipper(amin,None)
    g=pv.read(f"vol{k}.vtk"); rv=g.contour([0.03],scalars='rev')
    p=pv.Plotter(off_screen=True,window_size=window); p.set_background('white')
    if rv.n_points:
        rv=clip(rv.smooth_taubin(n_iter=40)).connectivity('largest')
        p.add_mesh(rv,color='black',lighting=False)
    p.camera_position=camera; p.screenshot(f"mask_{tag}_{k}.png"); p.close()
if __name__=="__main__":
    k=int(sys.argv[1]); tag=sys.argv[2]
    camera=view('front_valve')[0] if tag=='front' else swivel(*[float(x) for x in tag.replace('sw','').split('e')])
    render_mask(k,tag,camera)
    if len(sys.argv)>3: render_v4(k,tag,camera,ghost=False)
