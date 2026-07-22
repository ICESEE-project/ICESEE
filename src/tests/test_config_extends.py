from ICESEE.config.config_loader import load_yaml_to_dict


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
