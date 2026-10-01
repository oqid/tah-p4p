"""Regression checks for geometry-dependent streamline plot framing."""
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import proj3d

from haem5 import plot_streamlines_3d, _select_plot_streamlines


class StreamlineFramingTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_different_geometries_fit_with_labels_and_equal_scale(self):
        sizes = []
        for spans, origin in [
            ((1, 1, 1), (0, 0, 0)),
            ((30, 1, 2), (0, 0, 0)),
            ((1, 30, 2), (0, 0, 0)),
            ((1, 2, 30), (0, 0, 0)),
            ((1e-6, 2e-6, 3e-6), (0, 0, 0)),
            ((1, 2, 3), (1000, -2000, 3000)),
        ]:
            with self.subTest(spans=spans, origin=origin):
                t = np.linspace(0, 1, 40)
                points = np.column_stack((t, t ** 2, t ** 3)) * spans + origin
                df = pd.DataFrame(points, columns=["x", "y", "z"])
                df["shear"] = np.linspace(1, 100, len(df))
                before = df.copy(deep=True)
                fig = plot_streamlines_3d([df])
                ax, cax = fig.axes
                fig.canvas.draw()
                renderer = fig.canvas.get_renderer()
                sizes.append(tuple(fig.get_size_inches()))
                limits = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
                scales = ax.get_box_aspect() / np.diff(limits, axis=1).ravel()
                np.testing.assert_allclose(scales, np.full(3, scales[0]))
                self.assertEqual((ax.elev, ax.azim), (15, -80))
                pd.testing.assert_frame_equal(df, before)
                # Retain the currently accepted +90-degree X rotation.
                expected = np.column_stack((points[:, 0], -points[:, 2], points[:, 1]))
                np.testing.assert_allclose(ax.collections[0]._segments3d[0], expected[:2])
                corners = np.array(np.meshgrid(*limits)).reshape(3, -1)
                projected = np.column_stack(proj3d.proj_transform(*corners, ax.get_proj())[:2])
                pixels = ax.transData.transform(projected)
                plot_bounds = Bbox.union([
                    Bbox.from_extents(*pixels.min(axis=0), *pixels.max(axis=0)),
                    *(axis.get_tightbbox(renderer)
                      for axis in (ax.xaxis, ax.yaxis, ax.zaxis)),
                    *(axis.label.get_window_extent(renderer)
                      for axis in (ax.xaxis, ax.yaxis, ax.zaxis)),
                ])
                for bounds in (plot_bounds, cax.get_tightbbox(renderer)):
                    self.assertGreaterEqual(bounds.x0, 0)
                    self.assertGreaterEqual(bounds.y0, 0)
                    self.assertLessEqual(bounds.x1, fig.bbox.width)
                    self.assertLessEqual(bounds.y1, fig.bbox.height)
                self.assertLess(plot_bounds.x1, cax.get_tightbbox(renderer).x0)
                # At least one dimension should use most of the reserved area.
                self.assertGreater(max(plot_bounds.width / (0.795 * fig.bbox.width),
                                       plot_bounds.height / (0.88 * fig.bbox.height)), 0.85)
        self.assertGreater(len(set(sizes)), 1)

    def test_planar_plot_returns_a_figure(self):
        df = pd.DataFrame({"x": [0, 1, 2], "y": [0, 1, 0],
                           "z": [0, 0, 0], "shear": [1, 10, 100]})
        fig = plot_streamlines_3d([df])
        self.assertIsNotNone(fig)
        self.assertEqual(fig.axes[0].name, "rectilinear")

    def test_sparse_selection_is_unique_deterministic_and_keeps_extremes(self):
        dfs = [pd.DataFrame({"x": [0., 1.], "y": [0., 1.], "z": [0., 1.],
                             "shear": [float(i + 1)] * 2, "dt": [0., 1.]})
               for i in range(300)]
        dfs[10]["shear"] = [1., 1000.]
        dfs[20]["dt"] = [0., 100.]
        indices, records = _select_plot_streamlines(dfs, 8)
        self.assertEqual(len(set(indices)), 8)
        self.assertIn(10, indices)
        self.assertIn(20, indices)
        self.assertEqual(indices, _select_plot_streamlines(dfs, 8)[0])
        self.assertEqual(len(records), 8)
        sparse = plot_streamlines_3d(dfs, max_lines=8)
        default = plot_streamlines_3d(dfs)
        full = plot_streamlines_3d(dfs, max_lines=None)
        self.assertEqual(len(sparse.axes[0].collections), 8)
        self.assertEqual(len(default.axes[0].collections), 250)
        self.assertEqual(len(full.axes[0].collections), 300)
        np.testing.assert_allclose(sparse.axes[0].get_xlim(), full.axes[0].get_xlim())
        self.assertEqual(sparse.axes[0].collections[0].get_alpha(), 0.4)


if __name__ == "__main__":
    unittest.main()
