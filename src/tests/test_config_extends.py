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
