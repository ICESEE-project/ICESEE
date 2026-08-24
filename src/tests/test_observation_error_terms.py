import unittest

import numpy as np

from ICESEE.src.utils.localization import (
    active_observation_std,
    stochastic_observation_terms,
)


class ObservationErrorTermsTests(unittest.TestCase):
    def test_terms_are_reproducible_and_share_the_same_eta(self):
        rng = np.random.default_rng(3)
        ha = rng.normal(size=(4, 40))
        observations = np.arange(4.0)
        sigma = np.array([1.0, 2.0, 3.0, 4.0])

        hp1, eta1, dp1 = stochastic_observation_terms(
            ha, observations, sigma, 19
        )
        hp2, eta2, dp2 = stochastic_observation_terms(
            ha, observations, sigma, 19
        )

        np.testing.assert_allclose(hp1, hp2)
        np.testing.assert_allclose(eta1, eta2)
        np.testing.assert_allclose(dp1, dp2)
        np.testing.assert_allclose(np.mean(eta1, axis=1), 0.0, atol=1.0e-14)
        np.testing.assert_allclose(dp1, observations[:, None] + eta1 - ha)

    def test_active_std_accepts_both_file_orientations(self):
        rows = np.array([0, 2, 4])
        time_first = np.arange(1.0, 11.0).reshape(2, 5)
        expected = time_first[1, rows]

        np.testing.assert_allclose(
            active_observation_std({"error_R": time_first}, 1, rows),
            expected,
        )
        np.testing.assert_allclose(
            active_observation_std({"error_R": time_first.T}, 1, rows),
            expected,
        )


if __name__ == "__main__":
    unittest.main()
