"""STEP import, units, rotation, and optional overlay regression checks."""
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fluid_domain import (FluidDomainMesh, load_fluid_domain, camera_border_collection,
                          domain_surface_collection)
from haem5 import plot_streamlines_3d

try:
    import gmsh
except (ImportError, OSError):
    gmsh = None


class DomainTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_camera_contours_follow_view_and_ignore_face_winding(self):
        vertices = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0],
                             [0, -1., 0], [0, 0, 1.], [0, 0, -1.]])
        triangles = np.array([[4, 0, 2], [4, 2, 1], [4, 1, 3], [4, 3, 0],
                              [5, 2, 0], [5, 1, 2], [5, 3, 1], [5, 0, 3]])
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.set_proj_type("ortho")
        empty = np.empty((0, 2, 3))
        border = camera_border_collection(vertices, triangles, empty)
        reversed_border = camera_border_collection(vertices, triangles[:, ::-1], empty)
        ax.add_collection3d(border, autolim=False)
        ax.add_collection3d(reversed_border, autolim=False)
        ax.view_init(elev=20, azim=30)
        fig.canvas.draw()
        first = border._segments3d.copy()
        self.assertGreater(len(first), 0)
        self.assertLess(len(first), 12)  # Interior tessellation edges stay hidden.
        np.testing.assert_allclose(first, reversed_border._segments3d)
        ax.view_init(elev=20, azim=140)
        fig.canvas.draw()
        self.assertFalse(np.array_equal(first, border._segments3d))

    def test_shading_is_finite_preserves_opacity_and_tracks_camera(self):
        vertices = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [0, 0, 1.]])
        triangles = np.array([[0, 1, 2], [0, 1, 3]])
        for mode in ("matte", "glossy"):
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
            shell = domain_surface_collection(vertices, triangles, shading=mode, alpha=0.18)
            ax.add_collection3d(shell)
            fig.canvas.draw()
            colors = np.array(shell.get_facecolor())
            self.assertTrue(np.isfinite(colors).all())
            np.testing.assert_allclose(colors[:, 3], 0.18)
            ax.view_init(elev=60, azim=120)
            fig.canvas.draw()
            self.assertFalse(np.allclose(colors, shell.get_facecolor()))
        with self.assertRaises(ValueError):
            domain_surface_collection(vertices, triangles, shading="invalid")

    def test_overlay_rotates_border_and_frames_full_domain(self):
        vertices = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.], [0., 0., 6.]])
        mesh = FluidDomainMesh(vertices, np.array([[0, 1, 2], [0, 1, 3]]),
                               vertices[np.array([[0, 1], [0, 2], [0, 3]])])
        # A planar streamline still uses 3D when a domain is supplied.
        df = pd.DataFrame({"x": [0., 1.], "y": [0., 1.], "z": [0., 0.], "shear": [1., 10.]})
        fig = plot_streamlines_3d([df], fluid_domain=mesh)
        ax = fig.axes[0]
        self.assertEqual(ax.name, "3d")
        self.assertEqual(len(ax.collections), 3)
        self.assertEqual(ax.collections[0].get_alpha(), 0.06)
        np.testing.assert_allclose(ax.collections[1]._segments3d[2], [[0, 0, 0], [0, -6, 0]])
        self.assertGreater(ax.get_xlim()[1], 4)
        self.assertLess(ax.get_ylim()[0], -6)
        self.assertGreater(ax.get_zlim()[1], 5)
        self.assertEqual((ax.elev, ax.azim), (15, -80))
        fig_no_border = plot_streamlines_3d([df], fluid_domain=mesh, domain_border=False)
        self.assertEqual(len(fig_no_border.axes[0].collections), 2)
        with self.assertRaises(ValueError):
            plot_streamlines_3d([df], fluid_domain=mesh, domain_opacity=2)

    @unittest.skipIf(gmsh is None, "Optional Gmsh dependency not installed")
    def test_real_step_import_units_scale_and_session_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "box.stp"
            gmsh.initialize()
            try:
                gmsh.option.setNumber("General.Terminal", 0)
                gmsh.model.add("source_box")
                gmsh.model.occ.addBox(0, 0, 0, 50, 60, 70)
                gmsh.model.occ.synchronize()
                gmsh.write(str(path))
            finally:
                gmsh.finalize()
            mesh = load_fluid_domain(path)
            np.testing.assert_allclose(np.ptp(mesh.vertices, axis=0), [0.05, 0.06, 0.07])
            self.assertGreater(len(mesh.triangles), 0)
            self.assertGreater(len(mesh.border_segments), 0)
            self.assertFalse(gmsh.isInitialized())
            df = pd.DataFrame({"x": [0.01, 0.03], "y": [0.01, 0.04],
                               "z": [0.01, 0.05], "shear": [1., 10.]})
            output = Path(directory) / "streamlines_3d_step_test.png"
            plot_streamlines_3d([df], fluid_domain=path, save_path=output)
            self.assertGreater(output.stat().st_size, 0)
            self.assertTrue(output.with_name("streamline_selection_step_test.csv").is_file())
            gmsh.initialize()
            try:
                gmsh.model.add("existing_model")
                previous_size = gmsh.option.getNumber("Mesh.MeshSizeMax")
                scaled = load_fluid_domain(path, scale=2)
                np.testing.assert_allclose(np.ptp(scaled.vertices, axis=0), [0.1, 0.12, 0.14])
                self.assertEqual(gmsh.model.getCurrent(), "existing_model")
                self.assertEqual(gmsh.option.getNumber("Mesh.MeshSizeMax"), previous_size)
                self.assertEqual(gmsh.model.list(), ["", "existing_model"])
            finally:
                gmsh.finalize()

    def test_missing_file_and_unsupported_extension(self):
        with self.assertRaises(FileNotFoundError):
            load_fluid_domain("nonexistent_fluid_domain.stp")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.txt"
            path.write_text("bad")
            with self.assertRaises(ValueError):
                load_fluid_domain(path)


if __name__ == "__main__":
    unittest.main()
