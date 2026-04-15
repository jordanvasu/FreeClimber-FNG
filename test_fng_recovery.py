"""
Self-contained FNG recovery-speed test.
Lifts _height_traces, _detect_fng_series, and compute_fng verbatim from
scripts/detector.py and runs them against the pre-existing
climbing_1.filtered.csv -- no ffmpeg / trackpy needed.
"""

import types
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

# --------------------------------------------------------------------------
# Minimal stand-in for the detector object
# --------------------------------------------------------------------------
det = types.SimpleNamespace(
    debug               = False,
    frame_rate          = 29.0,
    pixel_to_cm         = 46.0,
    vials               = 6,
    name_nosuffix       = "example_other/ex_1/climbing_1",
    fng_enabled         = True,
    fng_smooth_window   = 5,
    fng_climb_thresh    = 0.10,
    fng_fall_thresh     = 0.10,
    fng_min_gap         = 5,
    fng_recovery_thresh = 0.05,
)

# --------------------------------------------------------------------------
# Copy of _height_traces
# --------------------------------------------------------------------------
def _height_traces(self):
    df = getattr(self, "df_filtered", None)
    if df is None or df.empty:
        return pd.DataFrame()
    H = (df.groupby(["frame", "vial"]).y.mean()
           .unstack("vial")
           .sort_index())
    return H

# --------------------------------------------------------------------------
# Copy of _detect_fng_series (with recovery speed)
# --------------------------------------------------------------------------
def _detect_fng_series(self, series,
                       smooth_window=None, climb_thresh=None,
                       fall_thresh=None, min_gap=None, recovery_thresh=None):
    smooth_window   = smooth_window   if smooth_window   is not None else 1
    climb_thresh    = climb_thresh    if climb_thresh    is not None else 0.0
    fall_thresh     = fall_thresh     if fall_thresh     is not None else 0.0
    min_gap         = min_gap         if min_gap         is not None else 1
    recovery_thresh = recovery_thresh if recovery_thresh is not None else 0.05

    frame_rate = float(getattr(self, "frame_rate", 1))
    if frame_rate <= 0:
        frame_rate = 1.0

    s = pd.Series(series).rolling(window=max(1, int(smooth_window)),
                                  center=True).mean().bfill().ffill()
    rng = s.max() - s.min()
    if not pd.notnull(rng) or rng == 0:
        return []
    nv = ((s - s.min()) / rng).values
    frame_index = s.index.to_numpy()
    peaks, _ = find_peaks(nv)

    events = []
    last_event_frame = -(10**9)
    w = max(1, int(smooth_window)) * 2

    px_to_cm = getattr(self, "pixel_to_cm", None)
    if px_to_cm is not None:
        try:
            px_to_cm = float(px_to_cm)
            if px_to_cm <= 0:
                px_to_cm = None
        except Exception:
            px_to_cm = None

    for p in peaks:
        if p - last_event_frame < min_gap:
            continue
        left_start  = max(0, p - w);  left_end  = p + 1
        right_start = p;              right_end = min(len(nv), p + w)
        left  = nv[left_start:left_end]
        right = nv[right_start:right_end]
        if len(left) < 2 or len(right) < 2:
            continue
        right_rel_idx = right.argmin()
        right_idx     = right_start + right_rel_idx
        rise      = nv[p] - left.min()
        drop_norm = nv[p] - right.min()

        if rise >= climb_thresh and drop_norm >= fall_thresh:
            drop_px = float(s.iloc[p]) - float(s.iloc[right_idx])
            drop_cm = drop_px / px_to_cm if px_to_cm is not None else np.nan

            fall_duration_frames = int(frame_index[right_idx]) - int(frame_index[p])
            fall_duration_sec    = fall_duration_frames / frame_rate

            next_after     = peaks[peaks > p]
            recovery_bound = int(next_after[0]) if len(next_after) > 0 else len(nv)
            rec_val                  = nv[right_idx] + recovery_thresh
            frame_recovery           = np.nan
            recovery_duration_frames = np.nan
            recovery_duration_sec    = np.nan

            for i in range(right_idx + 1, recovery_bound):
                if nv[i] >= rec_val:
                    frame_recovery           = int(frame_index[i])
                    recovery_duration_frames = frame_recovery - int(frame_index[right_idx])
                    recovery_duration_sec    = recovery_duration_frames / frame_rate
                    break

            events.append({
                "frame_peak"              : int(frame_index[p]),
                "frame_fall_start"        : int(frame_index[p]),
                "frame_fall_end"          : int(frame_index[right_idx]),
                "rise_norm"               : float(rise),
                "drop_norm"               : float(drop_norm),
                "drop_px"                 : float(drop_px),
                "drop_cm"                 : float(drop_cm),
                "fall_duration_frames"    : fall_duration_frames,
                "fall_duration_sec"       : fall_duration_sec,
                "frame_recovery"          : frame_recovery,
                "recovery_duration_frames": recovery_duration_frames,
                "recovery_duration_sec"   : recovery_duration_sec,
            })
            last_event_frame = p
    return events

# --------------------------------------------------------------------------
# Copy of compute_fng (with recovery speed)
# --------------------------------------------------------------------------
def compute_fng(self):
    _all_columns = [
        "vial", "event_idx",
        "frame_peak", "frame_fall_start", "frame_fall_end",
        "rise_norm", "drop_norm", "drop_px", "drop_cm",
        "fall_duration_frames", "fall_duration_sec",
        "frame_recovery", "recovery_duration_frames", "recovery_duration_sec",
    ]
    if not getattr(self, "fng_enabled", True):
        self.df_fng = pd.DataFrame(columns=_all_columns)
        return

    H        = _height_traces(self)
    path_fng = self.name_nosuffix + ".fng.csv"

    if H.empty:
        self.df_fng = pd.DataFrame(columns=_all_columns)
        self.df_fng.to_csv(path_fng, index=False)
        return

    records = []
    for vial in H.columns:
        evs = _detect_fng_series(
            self, H[vial],
            smooth_window   = getattr(self, "fng_smooth_window",    5),
            climb_thresh    = getattr(self, "fng_climb_thresh",   0.10),
            fall_thresh     = getattr(self, "fng_fall_thresh",    0.10),
            min_gap         = getattr(self, "fng_min_gap",           5),
            recovery_thresh = getattr(self, "fng_recovery_thresh", 0.05),
        )
        for idx, ev in enumerate(evs, start=1):
            rec_sec = ev["recovery_duration_sec"]
            records.append({
                "vial"                    : int(vial),
                "event_idx"               : idx,
                "frame_peak"              : int(ev["frame_peak"]),
                "frame_fall_start"        : int(ev["frame_fall_start"]),
                "frame_fall_end"          : int(ev["frame_fall_end"]),
                "rise_norm"               : round(ev["rise_norm"], 4),
                "drop_norm"               : round(ev["drop_norm"], 4),
                "drop_px"                 : round(ev["drop_px"], 4),
                "drop_cm"                 : round(ev["drop_cm"], 4),
                "fall_duration_frames"    : ev["fall_duration_frames"],
                "fall_duration_sec"       : round(ev["fall_duration_sec"], 4),
                "frame_recovery"          : ev["frame_recovery"],
                "recovery_duration_frames": ev["recovery_duration_frames"],
                "recovery_duration_sec"   : (round(rec_sec, 4)
                                             if not (isinstance(rec_sec, float)
                                                     and np.isnan(rec_sec))
                                             else np.nan),
            })

    self.df_fng = pd.DataFrame.from_records(records, columns=_all_columns)
    counts = (self.df_fng.groupby("vial").size()
                .rename("fng_count")
                .reindex(range(1, self.vials + 1), fill_value=0)
                .reset_index())
    self.df_fng_counts = counts
    self.df_fng.to_csv(path_fng, index=False)
    print(f"  --> Saved: {path_fng}")

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
print("Loading climbing_1.filtered.csv ...")
det.df_filtered = pd.read_csv("example_other/ex_1/climbing_1.filtered.csv")
print(f"  shape: {det.df_filtered.shape}")

print("\nRunning compute_fng() ...")
compute_fng(det)

print("\n=== df_fng columns ===")
print(det.df_fng.columns.tolist())

new_cols = ["fall_duration_frames", "fall_duration_sec",
            "frame_recovery", "recovery_duration_frames", "recovery_duration_sec"]
missing  = [c for c in new_cols if c not in det.df_fng.columns]
print("\nNew columns present:", ("ALL OK" if not missing else f"MISSING {missing}"))

print("\n=== df_fng (all rows) ===")
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 240)
pd.set_option("display.float_format", "{:.4f}".format)
print(det.df_fng.to_string(index=False))

print("\n=== Per-vial FNG counts ===")
print(det.df_fng_counts.to_string(index=False))

print("\n=== Sanity checks ===")
if not det.df_fng.empty:
    fd = det.df_fng["fall_duration_frames"]
    print(f"fall_duration_frames   : min={fd.min()}, max={fd.max()}, any_negative={(fd < 0).any()}")

    rd = det.df_fng["recovery_duration_frames"].dropna()
    if len(rd):
        print(f"recovery_duration_frames (non-NaN, n={len(rd)}): "
              f"min={rd.min()}, max={rd.max()}, any_negative={(rd < 0).any()}")
    else:
        print("recovery_duration_frames: all NaN (no recoveries found within bound)")

    frac = det.df_fng["frame_recovery"].notna().mean()
    print(f"Fraction of events with recovery found: {frac:.0%}")

    # frame_fall_end + recovery_duration == frame_recovery
    has_rec   = det.df_fng["frame_recovery"].notna()
    if has_rec.any():
        check = (det.df_fng.loc[has_rec, "frame_fall_end"]
                 + det.df_fng.loc[has_rec, "recovery_duration_frames"]
                 == det.df_fng.loc[has_rec, "frame_recovery"])
        print(f"frame_fall_end + recovery_duration_frames == frame_recovery: {check.all()}")

    # frame_fall_end - frame_fall_start == fall_duration_frames
    dur_check = (det.df_fng["frame_fall_end"] - det.df_fng["frame_fall_start"]
                 == det.df_fng["fall_duration_frames"])
    print(f"frame_fall_end - frame_fall_start == fall_duration_frames:      {dur_check.all()}")
else:
    print("df_fng is empty -- no events detected")
