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
