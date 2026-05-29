"""
Tests for per-fly tortuosity metrics (scripts/tortuosity.py).

Covers:
  * test_straight_line             -- T==1, S==1, theta==0 for a straight path.
  * test_right_angle_L             -- known T for an L-shape (2 == sqrt(2)/1 ?
                                       no -- L of equal arms: T = 2/sqrt(2) = sqrt(2)).
  * test_closed_loop               -- net_displacement == 0 -> straightness == 0,
                                       tortuosity == NaN (closed bout convention).
  * test_degenerate_short_track    -- <2 points -> all metrics NaN; coincident
                                       points -> path_length == 0, both ratios NaN.
  * test_compute_tortuosity_table  -- end-to-end: synthetic df_filtered with a
                                       'particle' column plus a synthetic df_fng
                                       produce one row per (vial, event, particle).
  * test_cohort_mode_no_particles  -- without a 'particle' column the function
                                       returns an empty, correctly-typed table.

No video decoding required; everything operates on tiny in-memory DataFrames.
"""

import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import tortuosity as tort  # noqa: E402


# ---------------------------------------------------------------------------
# bout_metrics: pure-function unit tests
# ---------------------------------------------------------------------------
def test_straight_line():
    xs = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    ys = np.zeros_like(xs)
    m = tort.bout_metrics(xs, ys)
    assert m['n_points'] == 5
    assert m['path_length'] == pytest.approx(4.0)
    assert m['net_displacement'] == pytest.approx(4.0)
    assert m['tortuosity'] == pytest.approx(1.0)
    assert m['straightness'] == pytest.approx(1.0)
    assert m['mean_turning_angle_rad'] == pytest.approx(0.0, abs=1e-12)


def test_right_angle_L():
    # Two equal-length arms at 90 degrees: total length = 2, net = sqrt(2).
    xs = np.array([0.0, 1.0, 1.0])
    ys = np.array([0.0, 0.0, 1.0])
    m = tort.bout_metrics(xs, ys)
    assert m['path_length'] == pytest.approx(2.0)
    assert m['net_displacement'] == pytest.approx(math.sqrt(2.0))
    assert m['tortuosity'] == pytest.approx(2.0 / math.sqrt(2.0))
    assert m['straightness'] == pytest.approx(math.sqrt(2.0) / 2.0)
    # Single 90-degree turn -> pi/2.
    assert m['mean_turning_angle_rad'] == pytest.approx(math.pi / 2.0)


def test_closed_loop():
    # Triangle returning to its origin: net displacement == 0.
    xs = np.array([0.0, 1.0, 1.0, 0.0])
    ys = np.array([0.0, 0.0, 1.0, 0.0])
    m = tort.bout_metrics(xs, ys)
    assert m['path_length'] == pytest.approx(2.0 + math.sqrt(2.0))
    assert m['net_displacement'] == pytest.approx(0.0)
    # Convention: tortuosity is NaN for D == 0 (would otherwise be infinite),
    # straightness is exactly 0.
    assert math.isnan(m['tortuosity'])
    assert m['straightness'] == pytest.approx(0.0)
    # Turning angles are still well defined.
    assert not math.isnan(m['mean_turning_angle_rad'])


def test_single_point_track():
    m = tort.bout_metrics([5.0], [5.0])
    assert m['n_points'] == 1
    assert math.isnan(m['path_length'])
    assert math.isnan(m['net_displacement'])
    assert math.isnan(m['tortuosity'])
    assert math.isnan(m['straightness'])
    assert math.isnan(m['mean_turning_angle_rad'])


def test_coincident_points():
    # All points at the same coordinate -> path_length == 0.
    xs = np.array([2.0, 2.0, 2.0, 2.0])
    ys = np.array([3.0, 3.0, 3.0, 3.0])
    m = tort.bout_metrics(xs, ys)
    assert m['n_points'] == 4
    assert m['path_length'] == pytest.approx(0.0)
    assert m['net_displacement'] == pytest.approx(0.0)
    assert math.isnan(m['tortuosity'])
    assert math.isnan(m['straightness'])
    # No non-zero step vectors -> turning angle is undefined.
    assert math.isnan(m['mean_turning_angle_rad'])
    # Vertical efficiency: clamped at 0 when path_length == 0.
    assert m['vertical_efficiency'] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# vertical_efficiency: primary climbing metric (Fix #1)
# ---------------------------------------------------------------------------
def test_vertical_efficiency_perfect_climb():
    """Pure vertical climb -> vertical_efficiency = 1.0."""
    ys = np.linspace(0.0, 10.0, 50)
    xs = np.zeros_like(ys)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(1.0, abs=1e-3)


def test_vertical_efficiency_horizontal_walk():
    """Pure horizontal walk -> vertical_efficiency = 0.0."""
    xs = np.linspace(0.0, 10.0, 50)
    ys = np.zeros_like(xs)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(0.0, abs=1e-3)


def test_vertical_efficiency_descending():
    """Net descent -> vertical_efficiency clamped at 0.0."""
    xs = np.zeros(50)
    ys = np.linspace(10.0, 0.0, 50)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(0.0, abs=1e-3)


def test_vertical_efficiency_diagonal():
    """45-degree diagonal climb -> vertical_efficiency = sin(45 deg) ~ 0.707."""
    xs = np.linspace(0.0, 10.0, 50)
    ys = np.linspace(0.0, 10.0, 50)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(math.sin(math.pi / 4.0), abs=1e-3)


# ---------------------------------------------------------------------------
# compute_tortuosity_table: integration unit tests on synthetic DataFrames
# ---------------------------------------------------------------------------
def _make_df_filtered():
    """Synthetic linked detections: two flies in vial 1 over 20 frames.
    Fly 100 climbs straight; fly 101 takes a zig-zag of identical net motion.
    """
    n = 20
    frames = np.arange(n)

    # Fly 100: straight downward y (climbing in image coords).
    fly_a = pd.DataFrame({
        'frame': frames, 't': frames / 25.0, 'vial': 1,
        'x': np.full(n, 100.0),
        'y': np.linspace(0.0, 19.0, n),
        'particle': 100,
    })

    # Fly 101: alternating +/-1 x perturbation, same start/end y.
    x_zig = np.full(n, 110.0) + np.where(frames % 2 == 0, 0.0, 1.0)
    fly_b = pd.DataFrame({
        'frame': frames, 't': frames / 25.0, 'vial': 1,
        'x': x_zig,
        'y': np.linspace(0.0, 19.0, n),
        'particle': 101,
    })

    return pd.concat([fly_a, fly_b], ignore_index=True)


def _make_df_fng():
    """Two FNG events in vial 1 to drive bout windows."""
    return pd.DataFrame([
        {'vial': 1, 'event_idx': 1, 'frame_peak': 9, 'frame_fall_end': 12},
        {'vial': 1, 'event_idx': 2, 'frame_peak': 19, 'frame_fall_end': 19},
    ])


def test_compute_tortuosity_table_per_event_per_particle():
    df = _make_df_filtered()
    df_fng = _make_df_fng()
    out = tort.compute_tortuosity_table(df, df_fng)

    # 2 events x 2 particles = 4 rows.
    assert len(out) == 4
    assert set(out['event_idx']) == {1, 2}
    assert set(out['particle']) == {100, 101}

    # Fly 100 took a straight path in every bout -> tortuosity == 1.
    straight = out[out['particle'] == 100]
    assert np.allclose(straight['tortuosity'].to_numpy(), 1.0)
    assert np.allclose(straight['straightness'].to_numpy(), 1.0)

    # Fly 101's zig-zag has longer path than net displacement -> T > 1.
    zigzag = out[out['particle'] == 101]
    assert (zigzag['tortuosity'] > 1.0).all()
    assert (zigzag['straightness'] < 1.0).all()


def test_compute_tortuosity_table_no_particle_column():
    """Cohort-mode df_filtered (no 'particle' column) returns an empty table."""
    df = _make_df_filtered().drop(columns=['particle'])
    out = tort.compute_tortuosity_table(df, _make_df_fng())
    assert out.empty
    assert list(out.columns) == tort.TORTUOSITY_BOUT_COLUMNS


def test_compute_tortuosity_table_no_fng_events():
    df = _make_df_filtered()
    out = tort.compute_tortuosity_table(df, pd.DataFrame(columns=['vial', 'event_idx', 'frame_peak', 'frame_fall_end']))
    assert out.empty
    assert list(out.columns) == tort.TORTUOSITY_BOUT_COLUMNS


def test_compute_tortuosity_table_missing_required_column_raises():
    df = _make_df_filtered().drop(columns=['x'])
    with pytest.raises(ValueError):
        tort.compute_tortuosity_table(df, _make_df_fng())


# ---------------------------------------------------------------------------
# Integration with the detector hook -- writes <video>.tortuosity.csv
# ---------------------------------------------------------------------------
import types  # noqa: E402

sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402


def _bind(det, *names):
    for name in names:
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def test_detector_compute_tortuosity_writes_csv(tmp_path):
    """detector.compute_tortuosity() writes BOTH *.tortuosity_bouts.csv and
    *.tortuosity_particle.csv next to the other outputs."""
    df = _make_df_filtered()
    df_fng = _make_df_fng()

    det = types.SimpleNamespace(
        debug=False,
        df_filtered=df,
        df_fng=df_fng,
        name_nosuffix=os.path.join(str(tmp_path), "synthetic"),
    )
    _bind(det, "compute_tortuosity")
    det.compute_tortuosity()

    bouts_path = os.path.join(str(tmp_path), "synthetic.tortuosity_bouts.csv")
    particle_path = os.path.join(str(tmp_path), "synthetic.tortuosity_particle.csv")
    assert os.path.exists(bouts_path)
    assert os.path.exists(particle_path)

    bouts = pd.read_csv(bouts_path)
    assert list(bouts.columns) == tort.TORTUOSITY_BOUT_COLUMNS
    assert len(bouts) == 4
    assert 'vertical_efficiency' in bouts.columns

    particle = pd.read_csv(particle_path)
    assert list(particle.columns) == tort.TORTUOSITY_PARTICLE_COLUMNS
    assert 'median_vertical_efficiency' in particle.columns
    assert len(particle) == 2  # two particles -> two rows


def test_detector_compute_tortuosity_cohort_mode_writes_empty_csv(tmp_path):
    """Cohort mode (no 'particle' column) still emits header-only csvs so
    downstream tooling doesn't have to special-case missing files."""
    df = _make_df_filtered().drop(columns=['particle'])
    det = types.SimpleNamespace(
        debug=False,
        df_filtered=df,
        df_fng=_make_df_fng(),
        name_nosuffix=os.path.join(str(tmp_path), "cohort"),
    )
    _bind(det, "compute_tortuosity")
    det.compute_tortuosity()

    bouts_path = os.path.join(str(tmp_path), "cohort.tortuosity_bouts.csv")
    particle_path = os.path.join(str(tmp_path), "cohort.tortuosity_particle.csv")
    assert os.path.exists(bouts_path)
    assert os.path.exists(particle_path)

    bouts = pd.read_csv(bouts_path)
    assert bouts.empty
    assert list(bouts.columns) == tort.TORTUOSITY_BOUT_COLUMNS
    particle = pd.read_csv(particle_path)
    assert particle.empty
    assert list(particle.columns) == tort.TORTUOSITY_PARTICLE_COLUMNS
