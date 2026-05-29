#!/usr/bin/env python
# -*- coding: utf-8 -*-

## File name : tortuosity.py
## Purpose   : Per-fly trajectory tortuosity metrics for FreeClimber-FNG.
##             Computed from linked tracks (analysis_mode == 'individual')
##             once per FNG climbing bout.

"""
Per-fly tortuosity metrics for FreeClimber-FNG.

Three metrics are computed per climbing bout, per fly:

  1. tortuosity (T)            classical tortuosity = path_length / net_displacement
                               T >= 1 for any non-zero net displacement; larger
                               values mean a more meandering path.
  2. straightness (S)          straightness index   = net_displacement / path_length
                               S in [0, 1]; S = 1 is a perfectly straight path.
                               (S == 1 / T; we record both because each is
                                conventional in different sub-fields and the
                                inverse is undefined when path_length == 0.)
  3. mean_turning_angle (rad)  mean absolute turning angle between consecutive
                               step vectors, in radians. 0 == straight, pi/2 ==
                               random walk, pi == reversals every step. This is
                               a path-shape metric independent of T / S and is
                               not destroyed by zero net displacement (closed
                               loops).

  ---
  The task brief described two metrics in full (T and S) and was truncated
  before the third was named; mean absolute turning angle is the conservative
  choice. It is the standard third-axis metric in movement ecology, well
  defined for closed bouts (where T and S degenerate), and adds genuinely new
  information rather than restating T or S.

A "climbing bout" is the contiguous frame window from the previous FNG event's
fall-end (or the start of the recording, for the first event) up to and
including the current event's peak. Bouts are defined per vial -- the same
vial-level FNG events that drive compute_fng() also drive bout windows here --
and metrics are then computed per linked particle (fly) inside each window.
"""

import numpy as np
import pandas as pd

TORTUOSITY_COLUMNS = [
    'vial', 'event_idx', 'particle',
    'frame_start', 'frame_end', 'n_points',
    'path_length', 'net_displacement',
    'tortuosity', 'straightness', 'mean_turning_angle_rad',
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


def mean_absolute_turning_angle(xs, ys):
    """Mean absolute turning angle (radians) between consecutive step vectors.

    For a polyline with n points there are n-1 steps and n-2 turning angles.
    Steps with zero length are skipped (they have no defined direction).
    Returns NaN if fewer than 2 well-defined turning angles can be computed.
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
    """Compute the three tortuosity metrics for a single bout polyline.

    Inputs:
      xs, ys : 1-D arrays of equal length, ordered by frame (ascending).

    Returns:
      dict with keys: n_points, path_length, net_displacement, tortuosity,
                     straightness, mean_turning_angle_rad.

    Degenerate handling:
      * fewer than 2 points -> all metrics NaN.
      * net_displacement == 0 (closed loop) -> tortuosity is NaN
        (path_length / 0 has no useful value); straightness is 0.
      * path_length == 0 (all points coincide) -> straightness is NaN
        and tortuosity is NaN.
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

    return {
        'n_points': n,
        'path_length': L,
        'net_displacement': D,
        'tortuosity': T,
        'straightness': S,
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
    """Build the per-(vial, event, particle) tortuosity table.

    Inputs:
      df_filtered : DataFrame with columns frame, vial, x, y, particle (the
                    linked output written by detector.link_trajectories).
                    A 'particle' column is required; if absent, an empty
                    table is returned (cohort mode has no per-fly tracks).
      df_fng      : DataFrame with columns vial, event_idx, frame_peak,
                    frame_fall_end -- the event table written by compute_fng().

    Returns:
      DataFrame with columns TORTUOSITY_COLUMNS. Empty (correctly-typed
      headers only) when there is nothing to compute.
    """
    if df_filtered is None or df_filtered.empty or 'particle' not in df_filtered.columns:
        return pd.DataFrame(columns=TORTUOSITY_COLUMNS)
    if df_fng is None or df_fng.empty:
        return pd.DataFrame(columns=TORTUOSITY_COLUMNS)

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
                    'mean_turning_angle_rad': m['mean_turning_angle_rad'],
                })

    if not records:
        return pd.DataFrame(columns=TORTUOSITY_COLUMNS)
    return pd.DataFrame.from_records(records, columns=TORTUOSITY_COLUMNS)
