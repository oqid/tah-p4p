import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import haem5
import haem5_paper_comparison as p

class PaperComparisonTests(unittest.TestCase):
    def test_uniform_stress_hi3_matches_power_law(self):
        c=haem5.POWER_LAW_CONSTANTS["GW"]
        _, hi, _=p.coefficients_on_path([0,100,100],[0,.04,.06],c)
        self.assertAlmostEqual(hi,c["A"]*100**c["alpha"]*.1**c["beta"],places=12)

    def test_hi3_segment_refinement_and_order(self):
        c=haem5.POWER_LAW_CONSTANTS["HO"]
        a=p.coefficients_on_path([10,100],[.2,.05],c)
        b=p.coefficients_on_path([100,10],[.05,.2],c)
        refined=p.coefficients_on_path([10,10,100,100],[.1,.1,.025,.025],c)
        self.assertAlmostEqual(a[1],b[1],places=12)
        self.assertAlmostEqual(a[1],refined[1],places=12)
        self.assertNotAlmostEqual(a[0],b[0],places=8)

    def test_weighting_is_only_an_aggregation(self):
        hi=np.array([1.,.5])
        mean,weighted=p.aggregate_hi(hi,[1.,2.])
        self.assertEqual(mean,.75)
        self.assertAlmostEqual(weighted,2/3)
        np.testing.assert_array_equal(hi,[1.,.5])
        self.assertEqual(p.aggregate_hi(hi,[2.,2.]),(.75,.75))
        self.assertTrue(np.isnan(p.aggregate_hi(hi,[0.,0.])[1]))

    def test_linear_density_normalizes_and_handles_zero(self):
        x=np.linspace(0,10,5000)
        density=p.reflection_kde(np.array([0.,1.,2.,3.]),x,.2)
        self.assertAlmostEqual(np.trapezoid(density,x),1.,places=5)

    def test_top_share_and_invalid_inputs(self):
        x,y=p.top_contribution([7.,1.,1.,1.])
        self.assertEqual(x[1],25)
        self.assertEqual(y[1],70)
        self.assertEqual(y[-1],100)
        with self.assertRaises(ValueError):
            p.aggregate_hi([1.,2.],[1.,-1.])

    def test_all_figures_and_paired_device_metrics(self):
        cases=[]
        rng=np.random.default_rng(4)
        for group in range(2):
            for phase in (25.,75.):
                n=40
                base=np.exp(rng.normal(-8,.5,n))*(group+1)
                d=pd.DataFrame({"SA_dyne_s_cm2":np.exp(rng.normal(0,.3,n)),
                                "transit_ms":rng.uniform(50,500,n),
                                "max_stress_Pa":rng.uniform(20,150,n),
                                "mean_stress_Pa":rng.uniform(.1,5,n),
                                "path_mean_speed_m_s":rng.uniform(.1,2,n),
                                "mean_inlet_velocity":rng.uniform(.1,1,n)})
                for model in p.MODELS:
                    d[model+"_HI2_percent"]=base*.5
                    d[model+"_HI3_percent"]=base
                a={f"p{q}":np.array([[0,0,0],[.1,100,q*1e-6],[.2,1,q*2e-6]]) for q in (50,95,99)}
                a.update(cloud=np.column_stack((rng.random((80,3))*.05,rng.uniform(.1,100,80))),
                         time_bands_ms=np.array([10,20,30,40,50]),dose_bands_percent=np.array([1,4,10,30,55]))
                cases.append(dict(label=f"device {group}",group=group,phase=phase,bpm=120,cardiac_output_l_min=None,
                                  n_paths=n,data=d,arrays=a,color=p.COLORS[int(phase==75)],
                                  style="-" if group==0 else "--",display=f"device {group} | {phase}",
                                  path_id_p50=1,path_id_p95=2,path_id_p99=3))
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            p.render(cases,out,{"comparison_context":"Synthetic validation"})
            self.assertEqual(len(list(out.glob("*.png"))),12)
            ratios=pd.read_csv(out/"matched_phase_ratios.csv")
            self.assertEqual(ratios.phase.tolist(),[25.,75.])
            self.assertTrue((out/"gallery.html").exists())

if __name__=="__main__":
    unittest.main()
