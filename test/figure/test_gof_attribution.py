import pytest
import pandas as pd

from seismic_waveform_factory.figure.generate_figure_synthetics import (
    build_src_lookup_for_plots,
)


class MockWaveformPlot:
    """Mock WaveformFigureGenerator for testing"""

    def __init__(self, plt_id, synthetics, enabled=True):
        self.plt_id = plt_id
        self.plt_cfg = {"synthetics": synthetics}
        self.enabled = enabled


@pytest.mark.parametrize(
    "synthetics_config,plot_synthetics,expected_order",
    [
        # Test case 1: Single synthetic (axitra_1)
        (
            [
                {"name": "axitra_1", "source_files": ["file_a1.h5"]},
                {"name": "seissol_1", "outputs": ["output_s1"]},
            ],
            ["axitra_1"],
            [("axitra_1", "file_a1.h5")],
        ),
        # Test case 2: Two synthetics in plot order [seissol_1, axitra_1]
        (
            [
                {"name": "axitra_1", "source_files": ["file_a1.h5", "file_a2.h5"]},
                {"name": "seissol_1", "outputs": ["output_s1"]},
            ],
            ["axitra_1", "seissol_1"],
            [
                ("axitra_1", "file_a1.h5"),
                ("axitra_1", "file_a2.h5"),
                ("seissol_1", "output_s1"),
            ],
        ),
        # Test case 3: Two synthetics in REVERSE of config order [axitra_1, seissol_1]
        (
            [
                {"name": "axitra_1", "source_files": ["file_a1.h5"]},
                {"name": "seissol_1", "outputs": ["output_s1"]},
            ],
            ["seissol_1", "axitra_1"],
            [("seissol_1", "output_s1"), ("axitra_1", "file_a1.h5")],
        ),
        # Test case 4: Multiple files per synthetic
        (
            [
                {"name": "seissol_1", "outputs": ["out_s1", "out_s2"]},
                {"name": "axitra_1", "source_files": ["file_a1.h5"]},
            ],
            ["seissol_1", "axitra_1"],
            [
                ("seissol_1", "out_s1"),
                ("seissol_1", "out_s2"),
                ("axitra_1", "file_a1.h5"),
            ],
        ),
    ],
)
def test_gof_attribution_order(synthetics_config, plot_synthetics, expected_order):
    """Test that src_loop_up is built in correct plot order, not config order"""

    # Setup
    cfg = {"synthetics": synthetics_config}
    wf_plot = MockWaveformPlot(plt_id=0, synthetics=plot_synthetics, enabled=True)
    wf_plots = [wf_plot]

    # Execute
    src_loop_up = build_src_lookup_for_plots(wf_plots, cfg)

    # Assert: Check that the order matches plot config, not synthetics config
    assert src_loop_up["0"] == expected_order, (
        f"Expected {expected_order}, but got {src_loop_up['0']}. "
        f"Order should match plot synthetics {plot_synthetics}, not config order."
    )


def test_gof_column_attribution():
    """Test that GOF column indices correctly map to synthetic files"""

    # Setup: 2 synthetics, axitra has 1 file, seissol has 1 output
    cfg = {
        "synthetics": [
            {"name": "seissol_1", "outputs": ["output_s1"]},
            {"name": "axitra_1", "source_files": ["file_a1.h5"]},
        ]
    }
    wf_plot = MockWaveformPlot(
        plt_id=0, synthetics=["seissol_1", "axitra_1"], enabled=True
    )
    wf_plots = [wf_plot]

    # Build lookup
    src_loop_up = build_src_lookup_for_plots(wf_plots, cfg)

    # Simulate GOF dataframe with column names like: g0_E0, g0_E1, g0_Z0, g0_Z1
    # Model index 0 = seissol_1, Model index 1 = axitra_1
    df_station_average = pd.DataFrame(
        {
            "gofa_name": ["g0_E0", "g0_E1", "g0_Z0", "g0_Z1"],
            "gofa": [0.5, 0.6, 0.7, 0.8],
        }
    )

    # Extract plot_id and file_id from column names (as done in main)
    plot_id = (
        df_station_average["gofa_name"].str.extract(r"(\d+)", expand=False).astype(int)
    )
    file_id = (
        df_station_average["gofa_name"].str.extract(r"(\d+)$", expand=False).astype(int)
    )

    # Map to actual synthetics
    point_srcs = [
        src_loop_up[f"{p_id}"][f_id] for (p_id, f_id) in zip(plot_id, file_id)
    ]

    # Assert correct attribution
    expected = [
        ("seissol_1", "output_s1"),  # g0_E0: model 0
        ("axitra_1", "file_a1.h5"),  # g0_E1: model 1
        ("seissol_1", "output_s1"),  # g0_Z0: model 0
        ("axitra_1", "file_a1.h5"),  # g0_Z1: model 1
    ]
    assert point_srcs == expected, (
        f"GOF columns not attributed correctly. "
        f"Expected {expected}, got {point_srcs}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
