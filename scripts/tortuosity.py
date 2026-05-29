#!/usr/bin/env python
# -*- coding: utf-8 -*-

## File name : tortuosity.py
## Purpose   : Per-fly trajectory tortuosity metrics for FreeClimber-FNG.
##             Computed from linked tracks (analysis_mode == 'individual').

"""
Per-fly tortuosity metrics for FreeClimber-FNG.

Four metrics are computed per climbing bout, per fly:

  1. tortuosity (T)            classical tortuosity = path_length / net_displacement
                               T >= 1 for any non-zero net displacement; larger
                               values mean a more meandering path.
  2. straightness (S)          straightness index   = net_displacement / path_length
                               S in [0, 1]; S = 1 is a perfectly straight path.
                               (S == 1 / T; we record both because each is
                                conventional in different sub-fields and the
                                inverse is undefined when path_length == 0.)
  3. vertical_efficiency       max(0, y_end - y_start) / path_length, in [0, 1].
                               1.0 == every unit of path translates into one unit
                               of upward displacement; 0.0 == no net climb (or a
                               net descent). The PRIMARY climbing metric.
  4. mean_turning_angle (rad)  AUXILIARY path-shape metric: mean absolute
                               turning angle between consecutive step vectors,
                               in radians. 0 == straight, pi/2 == random walk,
                               pi == reversals every step. Kept for closed-loop
                               bouts (where T and S degenerate) and as a
                               directional-persistence indicator independent of
                               T / S / vertical_efficiency.

Y-CONVENTION
  detector.step_5 inverts y via detector.invert_y (detector_fng.py:1492), so in
  self.df_filtered "climbing" means INCREASING y. Confirmed empirically against
  the cohort climbing-slope sign reported by detector.local_linear_regression
  on the synthetic FNG validation fixture: slope is positive (+4.7428 px/frame)
  during the climb, with mean y rising from ~7.8 (early frames) to ~155
  (late frames). vertical_efficiency relies on this convention.

UNITS
  vertical_efficiency is unit-free, but is defined for path_length and
  (y_end - y_start) measured in the SAME linear units. compute_tortuosity in
  detector_fng.py converts pixel coordinates to MILLIMETRES (using
  self.pixel_to_cm) before invoking these helpers, so reported path lengths
  and net displacements are in mm.
"""

import numpy as np
import pandas as pd

TORTUOSITY_BOUT_COLUMNS = [
    'vial', 'event_idx', 'particle',
    'frame_start', 'frame_end', 'n_points',
    'path_length', 'net_displacement',
    'tortuosity', 'straightness', 'vertical_efficiency',
    'mean_turning_angle_rad',
]

# Backward-compat alias for any external callers / tests that imported
# TORTUOSITY_COLUMNS from earlier versions of this module.
TORTUOSITY_COLUMNS = TORTUOSITY_BOUT_COLUMNS

TORTUOSITY_PARTICLE_COLUMNS = [
    'vial', 'particle', 'n_bouts', 'median_vertical_efficiency',
]


def path_length(xs, ys):
    """Cumulative Euclidean distance along the (xs, ys) polyline.

    Returns 0.0 for sequences with fewer than 2 points.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2:
        return 0.0
    dx = np.diff(xs)
    dy = np.diff(ys)
    return float(np.sum(np.sqrt(dx * dx + dy * dy)))


def net_displacement(xs, ys):
    """Straight-line distance from first to last point of the (xs, ys) polyline.

    Returns 0.0 for sequences with fewer than 2 points.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2:
        return 0.0
    return float(np.hypot(xs[-1] - xs[0], ys[-1] - ys[0]))


def vertical_efficiency(xs, ys):
    """Fraction of (smoothed) path length that maps to net upward climb.

    Definition (per spec): vertical_efficiency = max(0, y_end - y_start) / L
    where L is the smoothed path length over (xs, ys) and (xs, ys) are in
    millimetres. Pure function: callers in compute_tortuosity smooth the
    coordinates (Savitzky-Golay, polyorder=2) and convert to mm before
    invoking this helper.

    Y-convention: in self.df_filtered, "climbing" means INCREASING y -- see
    detector.invert_y at detector_fng.py:745-757 and the call site at
    detector_fng.py:1492. Empirically confirmed against the sign of the
    cohort climbing slope from detector.local_linear_regression on the
    synthetic FNG validation fixture (positive slope during climb).

    Range: [0, 1]. Returns 0.0 when path_length == 0 or when y_end <= y_start.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2:
        return 0.0
    L = path_length(xs, ys)
    if L == 0.0:
        return 0.0
    rise = float(ys[-1] - ys[0])
    if rise <= 0.0:
        return 0.0
    eff = rise / L
    # Numerical guard: rise can never exceed total path length, but clip
    # to keep the contract [0, 1] watertight in the face of FP noise.
    if eff > 1.0:
        return 1.0
    return eff


def mean_absolute_turning_angle(xs, ys):
    """Mean absolute turning angle (radians) between consecutive step vectors.

    For a polyline with n points there are n-1 steps and n-2 turning angles.
    Steps with zero length are skipped (they have no defined direction).
    Returns NaN if fewer than 2 well-defined turning angles can be computed.

    Auxiliary path-shape metric: not the primary climbing metric (see
    vertical_efficiency) but kept for closed-loop bouts where T and S
    degenerate and as a directional-persistence indicator.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 3:
        return float('nan')

    vx = np.diff(xs)
    vy = np.diff(ys)
    norms = np.sqrt(vx * vx + vy * vy)

    # Drop zero-length steps; a zero-length step has no defined direction.
    valid = norms > 0
    if valid.sum() < 2:
        return float('nan')
    vx = vx[valid]
    vy = vy[valid]
    norms = norms[valid]

    # Cosine of turning angle between consecutive step vectors.
    cos_theta = (vx[:-1] * vx[1:] + vy[:-1] * vy[1:]) / (norms[:-1] * norms[1:])
    # Floating-point can push |cos| just past 1, which would NaN arccos.
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    angles = np.arccos(cos_theta)
    return float(np.mean(np.abs(angles)))


def bout_metrics(xs, ys):
    """Compute the per-bout tortuosity metric bundle for a single polyline.

    Inputs:
      xs, ys : 1-D arrays of equal length, ordered by frame (ascending). For
               spec conformance compute_tortuosity passes smoothed (mm)
               coordinates; this pure function does not require it.

    Returns:
      dict with keys: n_points, path_length, net_displacement, tortuosity,
                     straightness, vertical_efficiency, mean_turning_angle_rad.

    Degenerate handling:
      * fewer than 2 points -> all metrics NaN (vertical_efficiency = 0).
      * net_displacement == 0 (closed loop) -> tortuosity is NaN
        (path_length / 0 has no useful value); straightness is 0.
      * path_length == 0 (all points coincide) -> straightness is NaN
        and tortuosity is NaN; vertical_efficiency is 0.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    n = int(xs.size)

    if n < 2:
        return {
            'n_points': n,
            'path_length': float('nan'),
            'net_displacement': float('nan'),
            'tortuosity': float('nan'),
            'straightness': float('nan'),
            'vertical_efficiency': 0.0,
            'mean_turning_angle_rad': float('nan'),
        }

    L = path_length(xs, ys)
    D = net_displacement(xs, ys)

    if L == 0.0:
        T = float('nan')
        S = float('nan')
    elif D == 0.0:
        # Closed loop: classical tortuosity diverges; straightness is exactly 0.
        T = float('nan')
        S = 0.0
    else:
        T = L / D
        S = D / L

    theta = mean_absolute_turning_angle(xs, ys)
    ve = vertical_efficiency(xs, ys)

    return {
        'n_points': n,
        'path_length': L,
        'net_displacement': D,
        'tortuosity': T,
        'straightness': S,
        'vertical_efficiency': ve,
        'mean_turning_angle_rad': theta,
    }


def _bout_windows_for_vial(df_fng_vial, frame_max):
    """Yield (event_idx, frame_start, frame_end) for each climbing bout in a vial.

    The bout window is [previous fall_end + 1, frame_peak]. The first bout
    starts at frame 0. Events are processed in ascending frame_peak order.
    """
    df = df_fng_vial.sort_values('frame_peak').reset_index(drop=True)
    prev_fall_end = -1
    for _, row in df.iterrows():
        start = max(0, prev_fall_end + 1)
        end = int(row['frame_peak'])
        if end >= start:
            yield int(row['event_idx']), int(start), int(end)
        prev_fall_end = int(row['frame_fall_end'])


def compute_tortuosity_table(df_filtered, df_fng):
    """Build the per-(vial, event, particle) tortuosity-bout table.

    Inputs:
      df_filtered : DataFrame with columns frame, vial, x, y, particle (the
                    linked output written by detector.link_trajectories).
                    A 'particle' column is required; if absent, an empty
                    table is returned (cohort mode has no per-fly tracks).
      df_fng      : DataFrame with columns vial, event_idx, frame_peak,
                    frame_fall_end -- the event table written by compute_fng().

    Returns:
      DataFrame with columns TORTUOSITY_BOUT_COLUMNS. Empty (correctly-typed
      headers only) when there is nothing to compute.
    """
    if df_filtered is None or df_filtered.empty or 'particle' not in df_filtered.columns:
        return pd.DataFrame(columns=TORTUOSITY_BOUT_COLUMNS)
    if df_fng is None or df_fng.empty:
        return pd.DataFrame(columns=TORTUOSITY_BOUT_COLUMNS)

    needed = {'frame', 'vial', 'x', 'y', 'particle'}
    missing = needed - set(df_filtered.columns)
    if missing:
        raise ValueError('df_filtered missing required column(s): %s' % sorted(missing))

    frame_max = int(df_filtered['frame'].max())
    records = []

    for vial, df_vial in df_filtered.groupby('vial'):
        ev_vial = df_fng[df_fng['vial'] == vial]
        if ev_vial.empty:
            continue

        for event_idx, fstart, fend in _bout_windows_for_vial(ev_vial, frame_max):
            window = df_vial[(df_vial['frame'] >= fstart) & (df_vial['frame'] <= fend)]
            if window.empty:
                continue

            for particle, dfp in window.groupby('particle'):
                dfp = dfp.sort_values('frame')
                xs = dfp['x'].to_numpy()
                ys = dfp['y'].to_numpy()
                m = bout_metrics(xs, ys)
                records.append({
                    'vial': int(vial),
                    'event_idx': int(event_idx),
                    'particle': int(particle),
                    'frame_start': int(dfp['frame'].iloc[0]),
                    'frame_end': int(dfp['frame'].iloc[-1]),
                    'n_points': m['n_points'],
                    'path_length': m['path_length'],
                    'net_displacement': m['net_displacement'],
                    'tortuosity': m['tortuosity'],
                    'straightness': m['straightness'],
                    'vertical_efficiency': m['vertical_efficiency'],
                    'mean_turning_angle_rad': m['mean_turning_angle_rad'],
                })

    if not records:
        return pd.DataFrame(columns=TORTUOSITY_BOUT_COLUMNS)
    return pd.DataFrame.from_records(records, columns=TORTUOSITY_BOUT_COLUMNS)


def compute_particle_table(bouts):
    """Build the per-particle aggregate table from a per-bout table.

    Inputs:
      bouts : DataFrame with the TORTUOSITY_BOUT_COLUMNS schema.

    Returns:
      DataFrame with columns TORTUOSITY_PARTICLE_COLUMNS:
        vial, particle, n_bouts, median_vertical_efficiency.

      The spec for the particle CSV currently names median_vertical_efficiency
      explicitly; the file's primary purpose is per-fly aggregation of the
      primary climbing metric.
    """
    if bouts is None or bouts.empty:
        return pd.DataFrame(columns=TORTUOSITY_PARTICLE_COLUMNS)

    grouped = (bouts
               .groupby(['vial', 'particle'], as_index=False)
               .agg(n_bouts=('vial', 'size'),
                    median_vertical_efficiency=('vertical_efficiency', 'median')))
    return grouped[TORTUOSITY_PARTICLE_COLUMNS]
