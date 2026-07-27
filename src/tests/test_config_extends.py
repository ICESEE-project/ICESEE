from ICESEE.config.config_loader import load_yaml_to_dict
from pathlib import Path
import yaml


def test_yaml_extends_deep_merges_sections(tmp_path):
    base = tmp_path / "base.yaml"
    base.write_text("section:\n  keep: 1\n  replace: 2\n")
    child = tmp_path / "child.yaml"
    child.write_text("extends: base.yaml\nsection:\n  replace: 3\n")

    assert load_yaml_to_dict(child) == {
        "section": {"keep": 1, "replace": 3}
    }


def test_reviewer_controls_inherit_hybrid_observation_design():
    root = (
        "applications/issm_model/examples/ISMIP_Choi/"
        "reviewer_experiments"
    )
    hybrid = load_yaml_to_dict(f"{root}/friction_inversion_hybrid.yaml")[
        "enkf-parameters"
    ]
    enkf_only = load_yaml_to_dict(f"{root}/friction_enkf_only.yaml")[
        "enkf-parameters"
    ]
    fixed = load_yaml_to_dict(f"{root}/wrong_friction_fixed.yaml")[
        "enkf-parameters"
    ]

    for control in (enkf_only, fixed):
        assert control["bed_obs_snapshot"] == hybrid["bed_obs_snapshot"]
        assert control["bed_obs_stride"] == hybrid["bed_obs_stride"]
        assert control["bed_update_domain"] == "grounded_only"
        assert control["inversion_flag"] == 0
    assert enkf_only["frozen_analysis_vars"] == []
    assert fixed["frozen_analysis_vars"] == ["coefficient"]


def test_tuned_low_prior_hybrid_delays_inversion_and_keeps_bed_sparse():
    path = (
        "applications/issm_model/examples/ISMIP_Choi/"
        "reviewer_experiments/"
        "friction_inversion_hybrid_low_prior_tuned.yaml"
    )
    tuned = load_yaml_to_dict(path)["enkf-parameters"]

    assert tuned["freq_obs"] == 1
    assert tuned["obs_max_time"] == 30
    assert tuned["bed_obs_snapshot"] == [2, 8, 14, 20, 24]
    assert tuned["inversion_flag"] == 1
    assert tuned["inversion_start_time"] == 6.0
    assert tuned["frozen_analysis_vars"] == ["coefficient"]
    assert tuned["min_friction"] == 1500
    assert tuned["max_friction"] == 4500
    assert tuned["bed_max_update_per_cycle"] == 20.0
    assert tuned["analysis_increment_limits"]["bed"] == 20.0
    assert tuned["initial_thickness_scale"] == 0.85
    assert tuned["initial_bed_offset_m"] == -80.0


def test_heterogeneous_hybrid_changes_prior_not_assimilation_design():
    root = (
        "applications/issm_model/examples/ISMIP_Choi/"
        "reviewer_experiments/"
    )
    tuned = load_yaml_to_dict(
        f"{root}friction_inversion_hybrid_low_prior_tuned.yaml"
    )["enkf-parameters"]
    mixed = load_yaml_to_dict(
        f"{root}friction_inversion_hybrid_heterogeneous.yaml"
    )["enkf-parameters"]

    # The robustness run must remain an apples-to-apples filter comparison.
    for key in (
        "freq_obs",
        "obs_max_time",
        "bed_obs_snapshot",
        "inversion_start_time",
        "min_friction",
        "max_friction",
        "bed_max_update_per_cycle",
    ):
        assert mixed[key] == tuned[key]

    assert mixed["initial_thickness_scale"] == 1.0
    assert mixed["initial_thickness_anomaly_fraction"] == 0.0
    assert mixed["initial_thickness_anomaly_m"] == 120.0
    assert mixed["initial_thickness_delta_min_m"] == -180.0
    assert mixed["initial_thickness_delta_max_m"] == 180.0
    assert mixed["initial_floating_thickness_anomaly_factor"] == 0.25
    assert mixed["initial_bed_gl_buffer_m"] == 25000.0
    assert mixed["initial_bed_offset_m"] == -80.0
    assert mixed["initial_bed_anomaly_m"] == 120.0
    assert mixed["initial_bed_delta_min_m"] == -250.0
    assert mixed["initial_bed_delta_max_m"] == 150.0
    assert mixed["initial_prior_length_x_m"] == 120000.0
    assert mixed["initial_prior_length_y_m"] == 40000.0
    assert mixed["initial_bed_background_domain"] == "grounded_only"

    check = load_yaml_to_dict(f"{root}heterogeneous_ic_check.yaml")[
        "enkf-parameters"
    ]
    assert check["initial_state_only"] is True
    assert check["generate_true_wrong_state_only"] is False
    assert check["data_path"] == "_reviewer_heterogeneous_ic_check"


def test_synchronized_ibf_wbf_ebf_profiles_differ_only_in_friction_method():
    root = (
        "applications/issm_model/examples/ISMIP_Choi/"
        "rebutal_experiments/"
    )
    common_path = Path(f"{root}param.yaml")
    common_raw = yaml.safe_load(common_path.read_text())
    assert "extends" not in common_raw

    for name in ("ibf", "wbf", "ebf"):
        child_raw = yaml.safe_load(Path(f"{root}param_{name}.yaml").read_text())
        assert child_raw["extends"] == "param.yaml"

    preflight = load_yaml_to_dict(f"{root}param_ic_check.yaml")["enkf-parameters"]
    assert preflight["initial_state_only"] is True
    assert preflight["generate_true_wrong_state_only"] is False
    assert (
        preflight["data_path"]
        == "_modelrun_datasets_rebuttal_ic_stretched_seed_bed_v6"
    )
    assert (
        preflight["initial_bed_background_domain"]
        == "grounded_plus_tapered_floating"
    )
    assert preflight["initial_bed_gl_buffer_m"] == 40000.0

    profiles = {
        name: load_yaml_to_dict(f"{root}param_{name}.yaml")
        for name in ("ibf", "wbf", "ebf")
    }

    common_keys = (
        "freq_obs",
        "obs_start_time",
        "obs_max_time",
        "bed_obs_snapshot",
        "observed_vars",
        "initial_thickness_anomaly_m",
        "initial_bed_offset_m",
        "initial_bed_background_domain",
        "initial_bed_gl_buffer_m",
        "initial_floating_bed_anomaly_factor",
        "initial_floating_bed_max_error_m",
        "initial_floating_bed_transition_m",
        "initial_floating_bed_flotation_margin_m",
        "initial_bed_smoothing_iterations",
        "initial_bed_smoothing_strength",
        "initial_bed_seed_max_x_m",
        "initial_bed_downstream_anomaly_factor",
        "seed",
    )
    ibf = profiles["ibf"]["enkf-parameters"]
    for profile in profiles.values():
        assert profile["modeling-parameters"]["num_years"] == 100
        enkf = profile["enkf-parameters"]
        for key in common_keys:
            assert enkf[key] == ibf[key]
        assert enkf["obs_max_time"] == 55
        assert enkf["bed_obs_snapshot"] == [
            2, 8, 14, 20, 24, 30, 36, 40, 46, 50
        ]
        assert (
            enkf["initial_bed_background_domain"]
            == "grounded_plus_tapered_floating"
        )
        assert enkf["initial_bed_gl_buffer_m"] == 40000.0
        assert enkf["initial_floating_bed_anomaly_factor"] == 0.95
        assert enkf["initial_floating_bed_max_error_m"] == 100.0
        assert enkf["initial_floating_bed_transition_m"] == 10000.0
        assert enkf["initial_floating_bed_flotation_margin_m"] == 5.0
        assert enkf["initial_bed_smoothing_iterations"] == 35
        assert enkf["initial_bed_smoothing_strength"] == 0.65
        assert enkf["initial_bed_seed_max_x_m"] == 300000.0
        assert enkf["initial_bed_downstream_anomaly_factor"] == 0.98

    assert ibf["data_path"] == "_modelrun_datasets_ibf_1"
    assert ibf["inversion_flag"] == 1
    assert ibf["inversion_start_time"] == 2.0
    assert ibf["frozen_analysis_vars"] == ["coefficient"]

    wbf = profiles["wbf"]["enkf-parameters"]
    assert wbf["data_path"] == "_modelrun_datasets_wbf_1"
    assert wbf["inversion_flag"] == 0
    assert wbf["frozen_analysis_vars"] == ["coefficient"]
    assert wbf["localized_vars"] == ["bed"]

    ebf = profiles["ebf"]["enkf-parameters"]
    assert ebf["data_path"] == "_modelrun_datasets_ebf_1"
    assert ebf["inversion_flag"] == 0
    assert ebf["frozen_analysis_vars"] == []
    assert ebf["localized_vars"] == ["bed", "coefficient"]
    assert ebf["analysis_increment_limits"]["coefficient"] == 300.0


def test_seaward_gl_prior_changes_only_the_initial_grounding_zone():
    root = (
        "applications/issm_model/examples/ISMIP_Choi/"
        "rebutal_experiments/seaward_gl_prior/"
    )
    raw = yaml.safe_load(Path(f"{root}param.yaml").read_text())
    assert raw["extends"] == "../param.yaml"

    profile = load_yaml_to_dict(f"{root}param.yaml")["enkf-parameters"]
    baseline = load_yaml_to_dict(
        "applications/issm_model/examples/ISMIP_Choi/"
        "rebutal_experiments/param.yaml"
    )["enkf-parameters"]
    assert profile["initial_gl_seaward_thickness_m"] == 500.0
    assert profile["initial_gl_seaward_width_m"] == 40000.0
    for key in (
        "initial_thickness_scale",
        "initial_bed_offset_m",
        "initial_bed_anomaly_m",
        "bed_obs_snapshot",
        "localization_radius",
        "analysis_increment_limits",
        "min_friction",
        "max_friction",
    ):
        assert profile[key] == baseline[key]

    check = load_yaml_to_dict(f"{root}param_ic_check.yaml")["enkf-parameters"]
    assert check["initial_state_only"] is True
    assert check["data_path"] == "_modelrun_datasets_rebuttal_ic_seaward_gl_v1"

    ibf = load_yaml_to_dict(f"{root}param_ibf.yaml")["enkf-parameters"]
    assert ibf["inversion_flag"] == 1
    assert ibf["frozen_analysis_vars"] == ["coefficient"]
