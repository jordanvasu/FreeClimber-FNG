"""
Tests for the failure-to-climb (FTC) measure (detector.compute_ftc) and the
FNG noise floor that keeps non-climbing vials from being scored as falls.

FTC is judged per fly over the whole clip: a fly with no climbing bout and no
fall failed to climb. No height line or time limit is involved.

Covers:
  * test_ftc_needs_individual_mode   -- cohort mode writes nothing.
  * test_ftc_particle_outcomes       -- climber / fng / ftc / unscored per fly,
    including a partial climb (a climber) and a partial climb that falls (fng).
  * test_ftc_detected_counts         -- per-vial counts from the tracks.
  * test_ftc_expected_counts_and_warnings -- with flies_per_vial, flies never
    seen moving count as failures; shortfalls and split tracks are flagged.
  * test_ftc_detections_dropped      -- flies that vanish are flagged.
  * test_ftc_disabled                -- ftc_enabled=False writes nothing.
  * test_fng_ignores_non_climbing_vial -- a vial of motionless, jittery flies
    produces no FNG events (it used to produce many).

Synthetic data only; no ffmpeg or video decoding is required.
"""

import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402

PX_PER_CM = 20.0
FPS = 10.0
N_FRAMES = 100


def _bind(det):
    for name in ("compute_ftc", "_climbing_bouts", "_flies_expected",
                 "_detect_fng_series", "_vial_label", "_relabel_vial_col"):
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def _make_det(df, tmp_path, **cfg):
    det = types.SimpleNamespace(
        debug=False, vials=int(df.vial.max()), vial_labels=None,
        pixel_to_cm=PX_PER_CM, frame_rate=FPS, n_frames=N_FRAMES,
        crop_0=0, crop_n=N_FRAMES,
        fng_smooth_window=5, fng_climb_thresh=0.10, fng_fall_thresh=0.10,
        fng_min_gap=5, ftc_eval_frames=5, ftc_min_coverage=0.8,
        df_filtered=df, file_details={"geno": "w1118"},
        name_nosuffix=str(tmp_path / "clip"),
    )
    for key, val in cfg.items():
        setattr(det, key, val)
    return _bind(det)


def _track(vial, particle, heights_cm, frames=None, jitter=0.0, seed=0):
    frames = np.arange(len(heights_cm)) if frames is None else np.asarray(frames)
    rng = np.random.default_rng(seed)
    y = np.asarray(heights_cm, dtype=float) * PX_PER_CM + rng.normal(0, jitter, len(frames))
    return pd.DataFrame({"frame": frames, "vial": vial, "particle": particle,
                         "x": 10.0 * vial, "y": y})


def _climb(to_cm, over=30, n=N_FRAMES):
    """Rise linearly from 0.1 cm to to_cm over 'over' frames, then hold."""
    return np.concatenate([np.linspace(0.1, to_cm, over), np.full(n - over, to_cm)])


def _scenario():
    """vial 1: 3 climbers; vial 2: 3 motionless flies; vial 3: a fly that
    climbs high then falls, a short climber (1 cm), a brief low fragment, and a
    short climb that ends in a fall."""
    parts = []
    for p in range(3):
        parts.append(_track(1, 10 + p, _climb(5.0), jitter=0.3, seed=p))
    for p in range(3):
        parts.append(_track(2, 20 + p, np.full(N_FRAMES, 0.2), jitter=0.5, seed=10 + p))
    up_then_down = np.concatenate([np.linspace(0.1, 5.0, 30), np.full(20, 5.0),
                                   np.linspace(5.0, 0.3, 5), np.full(45, 0.3)])
    partial_fall = np.concatenate([np.linspace(0.1, 1.2, 25), np.full(10, 1.2),
                                   np.linspace(1.2, 0.2, 5), np.full(60, 0.2)])
    parts.append(_track(3, 30, up_then_down, jitter=0.2, seed=20))
    parts.append(_track(3, 31, _climb(1.0), jitter=0.2, seed=21))
    parts.append(_track(3, 32, np.full(10, 0.2), frames=np.arange(10)))
    parts.append(_track(3, 33, partial_fall, jitter=0.2, seed=22))
    return pd.concat(parts, ignore_index=True)


def test_ftc_needs_individual_mode(tmp_path):
    det = _make_det(_scenario().drop(columns="particle"), tmp_path)
    det.compute_ftc()
    assert det.df_ftc.empty and det.df_ftc_particle.empty
    assert not os.path.exists(str(tmp_path / "clip.ftc.csv"))


def test_ftc_particle_outcomes(tmp_path):
    det = _make_det(_scenario(), tmp_path)
    det.compute_ftc()
    parts = det.df_ftc_particle.set_index("particle")

    assert (parts.loc[[10, 11, 12], "outcome"] == "climber").all()
    assert (parts.loc[[20, 21, 22], "outcome"] == "ftc").all()
    assert (parts.loc[[20, 21, 22], "n_climbing_bouts"] == 0).all()
    assert parts.loc[30, "outcome"] == "fng"
    # distance does not matter: a 1 cm climb is still a climb
    assert parts.loc[31, "outcome"] == "climber"
    assert parts.loc[32, "outcome"] == "unscored"
    # a short climb that ends in a fall is a fall
    assert parts.loc[33, "outcome"] == "fng" and parts.loc[33, "n_falls"] >= 1

    assert parts.loc[10, "latency_sec"] == pytest.approx(0.0, abs=0.3)
    assert parts.loc[[20, 32], "latency_sec"].isna().all()
    assert parts.loc[10, "max_rise_cm"] == pytest.approx(4.9, abs=0.3)
    assert parts.loc[20, "max_rise_cm"] < 0.2
    assert os.path.exists(str(tmp_path / "clip.ftc_particle.csv"))


def test_ftc_detected_counts(tmp_path):
    det = _make_det(_scenario(), tmp_path)
    det.compute_ftc()
    out = det.df_ftc.set_index("vial")
    assert list(out.method.unique()) == ["detected"]
    assert out.loc[1, ["n_tracks_climber", "ftc_count", "ftc_fraction"]].tolist() == [3, 0, 0.0]
    assert out.loc[2, ["n_tracks_ftc", "ftc_count", "ftc_fraction"]].tolist() == [3, 3, 1.0]
    assert out.loc[3, ["n_tracks_fng", "n_tracks_climber", "n_tracks_unscored",
                       "ftc_count"]].tolist() == [2, 1, 1, 0]
    assert out.clip_sec.iloc[0] == 10.0

    written = pd.read_csv(str(tmp_path / "clip.ftc.csv"))
    assert written.columns[0] == "geno"  # naming-convention details prepended


def test_ftc_expected_counts_and_warnings(tmp_path):
    det = _make_det(_scenario(), tmp_path, flies_per_vial=[3, 5, 2])
    det.compute_ftc()
    out = det.df_ftc.set_index("vial")
    assert list(out.method.unique()) == ["expected"]
    assert out.loc[1, "ftc_count"] == 0 and out.loc[1, "count_warning"] == ""
    # 5 loaded, none ever seen moving: all 5 failed, and the shortfall is flagged
    assert out.loc[2, ["ftc_count", "ftc_fraction"]].tolist() == [5, 1.0]
    assert out.loc[2, "count_warning"] == "fewer_detected_than_expected"
    # 3 moving tracks for 2 flies: a fly was split into several tracks
    assert "more_tracks_than_expected" in out.loc[3, "count_warning"]
    assert out.loc[3, "ftc_count"] == 0


def test_ftc_detections_dropped(tmp_path):
    df = _scenario()
    df = df[~((df.vial == 1) & (df.frame >= 50))]  # vial 1 flies vanish at frame 50
    det = _make_det(df, tmp_path)
    det.compute_ftc()
    row = det.df_ftc.set_index("vial").loc[1]
    assert row.n_tracks_climber == 3  # they climbed before vanishing
    assert row.n_detected_end == 0 and row.count_warning == "detections_dropped"


def test_ftc_disabled(tmp_path):
    det = _make_det(_scenario(), tmp_path, ftc_enabled=False)
    det.compute_ftc()
    assert det.df_ftc.empty
    assert not os.path.exists(str(tmp_path / "clip.ftc.csv"))


def test_fng_ignores_non_climbing_vial():
    """Motionless flies (sub-millimetre jitter) must not be scored as falls:
    the 0..1 normalisation used to stretch the jitter into fake climb/falls."""
    rng = np.random.default_rng(0)
    jitter = pd.Series(5 + rng.normal(0, 1.5, 300), index=range(300))
    det = types.SimpleNamespace(pixel_to_cm=26.0, frame_rate=30.0, n_frames=300)
    detect = types.MethodType(dfng.detector._detect_fng_series, det)

    assert len(detect(jitter)) == 0
    det.fng_min_range_cm = 0.0  # the old behaviour, for comparison
    assert len(detect(jitter)) > 5
