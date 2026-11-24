"""RTPS moments, ownership, fallback and shared analysis-control integration."""
import numpy as np
import pytest

from ICESEE.src.utils.ensemble_spread import relax_analysis_spread, spread_relaxation_settings


@pytest.mark.parametrize("alpha", [0., .25, .5, 1.])
def test_rtps_has_the_requested_sample_sd_and_preserves_analysis_mean(alpha):
    rng = np.random.default_rng(17)
    prior = rng.normal(size=(7, 32)) * np.arange(1, 8)[:, None]
    analysis = 100 + rng.normal(size=prior.shape)*.003
    original = analysis.copy()
    result = relax_analysis_spread(analysis, prior, np.arange(7), {"h": (0, 7)},
                                   {"analysis_spread_relaxation": alpha})
    expected = (1-alpha)*analysis.std(axis=1, ddof=1) + alpha*prior.std(axis=1, ddof=1)
    np.testing.assert_allclose(result.mean(axis=1), analysis.mean(axis=1), atol=3e-14)
    np.testing.assert_allclose(result.std(axis=1, ddof=1), expected, atol=2e-14)
    np.testing.assert_array_equal(analysis, original)
    if alpha == 0:
        assert result is analysis


def test_zero_analysis_spread_uses_prior_direction_without_random_draws():
    prior = np.array([[-1., 0., 1.], [5., 5., 5.]])
    analysis = np.array([[10., 10., 10.], [2., 2., 2.]])
    result = relax_analysis_spread(analysis, prior, np.arange(2), {"field": (0, 2)},
                                   {"analysis_spread_relaxation": .5})
    np.testing.assert_allclose(result, [[9.5, 10, 10.5], [2, 2, 2]])


def test_heterogeneous_fields_member_permutations_and_chunk_boundaries():
    rng = np.random.default_rng(4)
    prior = rng.normal(size=(9, 8)); analysis = 10+rng.normal(size=(9, 8))*.01
    rows = np.array([0, 2, 3, 5, 8, 9, 10, 11, 13])
    fields = {"h": (0, 4), "melt": (4, 10), "frozen": (10, 14)}
    options = {"analysis_spread_relaxation_factors": {"h": .3, "melt": .7}}
    whole = relax_analysis_spread(analysis, prior, rows, fields, options)
    chunks = np.concatenate([relax_analysis_spread(analysis[i:i+2], prior[i:i+2], rows[i:i+2], fields, options)
                             for i in range(0, len(rows), 2)])
    np.testing.assert_array_equal(whole, chunks)
    members = rng.permutation(8)
    permuted = relax_analysis_spread(analysis[:, members], prior[:, members], rows, fields, options)
    np.testing.assert_allclose(permuted, whole[:, members], atol=1e-14)
    np.testing.assert_array_equal(whole[rows >= 10], analysis[rows >= 10])


@pytest.mark.parametrize("value", [-.1, 1.1, np.nan, np.inf, True, "typo"])
def test_invalid_coefficients_are_refused(value):
    with pytest.raises(ValueError, match="spread-relaxation"):
        spread_relaxation_settings({"analysis_spread_relaxation": value}, ["h"])


def test_misspelled_fields_and_bad_layouts_are_refused():
    with pytest.raises(ValueError, match="unknown"):
        spread_relaxation_settings({"analysis_spread_relaxation_factors": {"thickness": .5}}, ["h"])
    with pytest.raises(ValueError, match="map"):
        spread_relaxation_settings({"analysis_spread_relaxation_factors": [1]}, ["h"])
    a = np.ones((2, 3))
    with pytest.raises(ValueError, match="cover"):
        relax_analysis_spread(a, a, np.array([0, 5]), {"h": (0, 2)}, {"analysis_spread_relaxation": .5})
    with pytest.raises(ValueError, match="two members"):
        relax_analysis_spread(a[:, :1], a[:, :1], np.arange(2), {"h": (0, 2)}, {"analysis_spread_relaxation": .5})


def test_shared_controls_keep_bounds_and_frozen_variables_authoritative():
    from ICESEE.src.parallelization._mpi_analysis_functions import apply_analysis_controls_local
    prior = np.array([[0., 10., 20.], [2., 3., 4.]])
    analysis = np.array([[7., 7., 7.], [90., 100., 110.]])
    options = dict(nd=2, vec_inputs=["melt", "frozen"], num_state_vars=1,
                   analysis_spread_relaxation=1., inflation_factor=1.,
                   analysis_bounds={"melt": [0, 10]}, frozen_analysis_vars=["frozen"])
    result = apply_analysis_controls_local(analysis, prior, np.arange(2), options)
    np.testing.assert_allclose(result[0], [0, 7, 10])
    np.testing.assert_array_equal(result[1], prior[1])


def test_weak_observation_limit_does_not_manufacture_uncertainty():
    rng = np.random.default_rng(15)
    prior = rng.normal(size=(3, 32))
    result = relax_analysis_spread(prior, prior, np.arange(3), {"field": (0, 3)},
                                   {"analysis_spread_relaxation": .8})
    np.testing.assert_allclose(result, prior, atol=1e-15)
