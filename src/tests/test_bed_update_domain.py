import unittest
import os
import tempfile

import h5py
import numpy as np

from ICESEE.src.utils.inference_plugin import (
    apply_bed_domain_gate_global,
    apply_bed_observation_anchor_global,
    resolve_analysis_cycle_time,
    resolve_bed_update_schedule,
)
from ICESEE.src.utils.stable_bed_inference import (
    apply_bed_regularized_correction,
)


class BedUpdateDomainTests(unittest.TestCase):
    def test_bed_schedule_is_resolved_without_private_runtime_state(self):
        active, is_snapshot, mode = resolve_bed_update_schedule({
            "inference_plugin_enabled": True,
            "bed_update_mode": "snapshots",
            "bed_snap_cols": [0, 3, 6],
            "km": 0,
        })
        self.assertTrue(active)
        self.assertTrue(is_snapshot)
        self.assertEqual(mode, "snapshots")

        active, is_snapshot, _ = resolve_bed_update_schedule({
            "inference_plugin_enabled": True,
            "bed_update_mode": "snapshots",
            "bed_snap_cols": [0, 3, 6],
            "km": 1,
            "_bed_update_active": True,  # stale private state is ignored
        })
        self.assertFalse(active)
        self.assertFalse(is_snapshot)

    def test_analysis_controls_use_requested_observation_time(self):
        cycle_time, model_time = resolve_analysis_cycle_time(
            {
                "t": np.array([0.0, 1.963636, 2.1818]),
                "obs_t": np.array([2.0]),
            },
            k=1,
            km=0,
        )
        self.assertEqual(cycle_time, 2.0)
        self.assertEqual(model_time, 1.963636)

    def test_analysis_cycle_time_falls_back_for_legacy_schedule(self):
        cycle_time, model_time = resolve_analysis_cycle_time(
            {"t": np.array([0.0, 1.5])}, k=1, km=0
        )
        self.assertEqual(cycle_time, 1.5)
        self.assertEqual(model_time, 1.5)

    def test_first_regularized_snapshot_uses_forecast_as_increment_reference(self):
        analysis = np.full((3, 2), 10.0)
        result = apply_bed_regularized_correction(
            analysis.copy(),
            ["bed"],
            3,
            {
                "physics_bed_inference": True,
                "mesh_coordinates": np.array(
                    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]
                ),
                "_bed_forecast_reference": np.zeros_like(analysis),
                "bed_spatial_regularization": 0.0,
                "bed_update_blend_factor": 1.0,
                "bed_enforce_below_surface": False,
                "bed_inference_start_time": 0.0,
            },
            model_time=2.0,
        )
        np.testing.assert_allclose(result, analysis)

    def test_grounded_only_gate_preserves_each_members_floating_bed(self):
        hdim = 3
        forecast = np.zeros((3 * hdim, 2))
        forecast[0:hdim, :] = [
            [100.0, 50.0],
            [20.0, 30.0],
            [80.0, 90.0],
        ]
        forecast[2 * hdim:3 * hdim, :] = [
            [-50.0, -60.0],
            [-40.0, -50.0],
            [-20.0, -30.0],
        ]
        analysis = forecast.copy()
        analysis[2 * hdim:3 * hdim, :] += 25.0

        result = apply_bed_domain_gate_global(
            analysis,
            forecast,
            ["Thickness", "Surface", "bed"],
            hdim,
            {
                "bed_update_domain": "grounded_only",
                "bed_snap_cols": [2],
                "km": 2,
                "di": 1.0,
            },
        )

        expected_bed = forecast[2 * hdim:3 * hdim, :].copy()
        grounded = (
            forecast[0:hdim, :]
            + forecast[2 * hdim:3 * hdim, :]
            > 0.0
        )
        expected_bed[grounded] += 25.0
        np.testing.assert_allclose(
            result[2 * hdim:3 * hdim, :], expected_bed
        )

    def test_gate_is_inactive_outside_bed_snapshots(self):
        forecast = np.zeros((4, 2))
        analysis = np.ones((4, 2))
        result = apply_bed_domain_gate_global(
            analysis.copy(),
            forecast,
            ["Thickness", "bed"],
            2,
            {
                "bed_update_domain": "grounded_only",
                "bed_snap_cols": [3],
                "km": 2,
            },
        )
        np.testing.assert_array_equal(result, analysis)

    def test_observed_only_gate_restores_every_unobserved_node(self):
        hdim = 4
        forecast = np.zeros((2 * hdim, 3))
        analysis = forecast.copy()
        analysis[hdim:, :] = 25.0
        support = np.zeros((hdim, 2), dtype=bool)
        support[[0, 2], 1] = True

        result = apply_bed_domain_gate_global(
            analysis,
            forecast,
            ["Thickness", "bed"],
            hdim,
            {
                "bed_update_domain": "observed_only",
                "bed_snap_cols": [1],
                "bed_mask_map_cols": {"bed": support},
                "km": 1,
            },
        )

        expected = np.zeros((hdim, 3))
        expected[[0, 2], :] = 25.0
        np.testing.assert_array_equal(result[hdim:, :], expected)

    def test_direct_bed_anchor_reduces_only_observed_innovations(self):
        hdim = 3
        analysis = np.zeros((2 * hdim, 2))
        observations = np.full((2 * hdim, 1), np.nan)
        observations[hdim + 1, 0] = 100.0

        result = apply_bed_observation_anchor_global(
            analysis,
            ["Thickness", "bed"],
            hdim,
            {
                "inference_plugin_enabled": True,
                "bed_update_mode": "snapshots",
                "bed_snap_cols": [0],
                "km": 0,
                "hu_obs_loaded": observations,
                "bed_observation_nudge_factor": 0.8,
            },
            stage="pre",
        )

        expected = np.zeros_like(analysis)
        expected[hdim + 1, :] = 80.0
        np.testing.assert_allclose(result, expected)

    def test_compact_bed_anchor_matches_dense_anchor(self):
        hdim = 4
        vec_inputs = ["Thickness", "Surface", "bed", "Friction"]
        analysis = np.arange(4 * hdim * 2, dtype=float).reshape(4 * hdim, 2)
        observations = np.full((4 * hdim, 2), np.nan)
        observations[2 * hdim + 1, 1] = -310.0
        observations[2 * hdim + 3, 1] = -725.0
        common = {
            "inference_plugin_enabled": True,
            "bed_update_mode": "snapshots",
            "bed_snap_cols": [1],
            "km": 1,
            "bed_observation_nudge_factor": 0.65,
        }
        dense = apply_bed_observation_anchor_global(
            analysis.copy(), vec_inputs, hdim,
            {**common, "hu_obs_loaded": observations},
            stage="pre",
        )

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "synthetic_obs.h5")
            rows = np.array([0, 2 * hdim + 1, 2 * hdim + 2, 2 * hdim + 3])
            values = np.zeros((rows.size, 2))
            active = np.zeros((rows.size, 2), dtype=np.uint8)
            values[1, 1] = -310.0
            values[2, 1] = 9999.0  # Inactive storage must never be anchored.
            values[3, 1] = -725.0
            active[[1, 3], 1] = 1
            with h5py.File(path, "w") as handle:
                handle.create_dataset("obs_indices", data=rows)
                handle.create_dataset("hu_obs_compact", data=values)
                handle.create_dataset("obs_active", data=active)
            compact = apply_bed_observation_anchor_global(
                analysis.copy(), vec_inputs, hdim,
                {**common, "synthetic_obs_file": path},
                stage="pre",
            )

        np.testing.assert_allclose(compact, dense)


if __name__ == "__main__":
    unittest.main()
