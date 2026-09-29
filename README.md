# FreeClimber-FNG

FreeClimber-FNG is a Python tool for video analysis of *Drosophila melanogaster* negative geotaxis (climbing) assays. Beyond the climbing velocity of the original [FreeClimber](https://github.com/adamspierer/FreeClimber) platform (Spierer et al., 2020), it measures three separate climbing outcomes:

- **Failed negative geotaxis (FNG)** — a fly climbs, then falls: fall events, fall distance, fall duration and recovery time
- **Failure to climb (FTC)** — a fly never ascends past a line within a time limit
- **Climbing quality** — per-fly tracks and how directly each fly climbs (tortuosity, straightness, vertical efficiency)

> **Note:** This repository is a fork of [adamspierer/FreeClimber](https://github.com/adamspierer/FreeClimber). The core particle detection and climbing velocity pipeline is the work of Adam N. Spierer and colleagues. FreeClimber-FNG adds the analyses above on top of that foundation. Development happens only here, at [jordanvasu/FreeClimber-FNG](https://github.com/jordanvasu/FreeClimber-FNG).

**Contents:** [What's new](#whats-new-in-freeclimber-fng) · [Installation](#installation) · [Usage](#usage) · [Output files](#output-files) · [Individual-fly tracking](#individual-fly-tracking-mode) · [FNG parameters](#fng-detection-parameters) · [Failure to climb](#failure-to-climb-ftc) · [Upgrade notes](#upgrade-notes-september-2026) · [Validation and tests](#synthetic-validation) · [Citing](#citing-this-work)

---

## What's New in FreeClimber-FNG

Added by Jordan Vasu (2025–2026):

- **FNG detection** — automated identification of climb-to-fall transitions (failed negative geotaxis events) in each vial
- **Fall measurements** — for each fall: the distance fallen (`drop_cm`, from the top of the climb to the bottom of the fall), the fall duration, and the recovery time until the fly starts climbing again (`recovery_duration_sec`)
- **Individual-fly tracking mode** — an optional mode that links per-frame detections into per-fly trajectories using [TrackPy](http://soft-matter.github.io/trackpy/) (including predictive linking), in addition to the default cohort (mean-position) analysis
- **Per-fly tortuosity / meandering metrics** — in individual mode, tortuosity, straightness, vertical efficiency and mean turning angle per fly per climbing bout, quantifying how directly (or erratically) each fly climbs
- **Failure to climb (FTC)** *(new, September 2026)* — a separate outcome from FNG: a fly that never passes a line (by default the top of the drawn ROI box) within a time limit, and never falls. Reported per vial in both modes, and per fly in individual mode with its time to reach the line. See [Failure to climb](#failure-to-climb-ftc).
- **Per-folder `vials.txt`** — override the vial count, name the physical vials, and give the number of flies loaded per vial
- **Reliability fixes** *(September 2026)* — step 6 no longer skips videos that have an empty vial or frame; FNG no longer counts jitter in non-climbing vials as falls; heights are measured from the vial floor; `.mov`/`.mp4` names are no longer truncated; the GUI reads the frame rate from the video and reports errors instead of closing. Some outputs change: see [Upgrade notes](#upgrade-notes-september-2026).

All analyses are implemented in `scripts/detector_fng.py` (plus `scripts/tortuosity.py`) and use the standard FreeClimber configuration file and batch workflow.

---

## Installation

FreeClimber-FNG uses the same environment and dependencies as the base FreeClimber platform. We recommend running in an Anaconda virtual environment.

**1. Create and activate a Python 3.8 environment:**

```bash
conda create -n freeclimber python=3.8
conda activate freeclimber
```

**2. Install dependencies:**

```bash
pip install numpy pandas scipy matplotlib trackpy ffmpeg-python wxPython
```

FFmpeg itself must also be installed and on your `PATH` (for example
`conda install -c conda-forge ffmpeg`). The current code is tested with Python
3.8, pandas 2.0, NumPy 1.24, SciPy 1.10, matplotlib 3.7, trackpy 0.7 and
wxPython 4.2.

**3. Clone this repository:**

```bash
git clone https://github.com/jordanvasu/FreeClimber-FNG.git
cd FreeClimber-FNG
```

For more on FFmpeg and wxPython setup, see the original [FreeClimber installation guide](https://github.com/adamspierer/FreeClimber#installing).

---

## Usage

FreeClimber-FNG is run in two stages: set up the analysis on one video in the GUI, which saves a configuration (`.cfg`) file, then batch-process every video with that file from the command line. See [TUTORIAL.md](TUTORIAL.md) for a step-by-step walkthrough; FNG detection is described in the [paper](paper.md).

**1. GUI — calibrate on one video and save a `.cfg`:**

```bash
pythonw ./scripts/FreeClimber_gui.py --video_file ./example/w1118_m_2_1.mov
```

Draw the region of interest (ROI) over the vials with its bottom edge at the
vial floor and its top edge at the finish line, set the detection parameters,
press **Test parameters** to check them, and press **Save configuration**. The
GUI steps are: 1 video and units (the frame rate is read from the video), 2 ROI,
3 spot detection, 4 frames and vials, 5 naming, 6 individual tracking and
tortuosity, 7 failure to climb.

**2. Command line — process every video in the project folder:**

```bash
python ./scripts/FreeClimber_main.py --config_file ./example/example.cfg --process_all
```

Other options: `--process_undone` (only videos without a `.slopes.csv`),
`--process_custom <file>.prc` (a list of video paths), `--no_concat` (skip the
project-level result files), `--optimization_plots` and `--debug` (stop at the
first error with a full traceback).

The batch runner walks `path_project` recursively and processes every file
ending in `file_suffix` (matched case-insensitively, so `mov` finds both `.mov`
and `.MOV`). A video that fails is logged to `<path_project>/log/skipped.log`
with the reason printed to the console, and the batch continues. Unless
`--no_concat` is given, the per-video results are then combined in
`path_project`:

| Project file | Combined from | Notes |
|---|---|---|
| `results.csv` | `*.slopes.csv` | climbing velocity per vial |
| `fng_results.csv` | `*.fng.csv` | one row per FNG event, plus a `video` column |
| `ftc_results.csv` | `*.ftc.csv` | failure to climb per vial, plus a `video` column |

Paths in the `.cfg` (`path_project`, `background_image`) are read verbatim, so
Windows paths with backslashes work; the GUI writes them with forward slashes.

### Per-folder vial count (`vials.txt`)

A single shared `.cfg` can be reused across folders whose videos have different
vial counts. Drop a plain-text file named `vials.txt` into a folder containing
just the integer count (e.g. `3`), and it overrides the `vials` value in the
`.cfg` for every video in that folder. Accepted contents: a non-blank,
non-`#`-comment line as a bare integer (`3`) or `vials=3` / `vials: 3`. A
missing or unparseable file leaves the `.cfg` value unchanged, so batch runs
never break.

#### Naming the surviving vials (`id =`)

When only some vials remain in a video — say the set `b31-35` was reduced to
just vials 33, 34 and 35 — add an `id =` line listing which physical vials are
present, left to right:

```
id = 33, 34, 35
```

These IDs label the vials in every per-vial output (the `vial_ID` in
`*.slopes.csv`/`results.csv`, and the `vial` column of `*.fng.csv`,
`*.tracks.csv` and the tortuosity CSVs) instead of the default positional
`1, 2, 3`. The vial count is taken from how many IDs are listed, so the `id =`
line alone is enough — a separate count line is optional, and if it disagrees
with the number of IDs the ID count wins. Accepted separators are `=` or `:`,
the key may be `id` or `ids`, and non-numeric labels are allowed too. A bare
comma-separated list with no key (`b26, b27, b29`) is also read as an ID line
(it used to be silently ignored). The vials are still detected by their
left-to-right position in the frame, so list the IDs in that same left-to-right
order.

#### Flies loaded per vial (`n =`)

An optional `n =` line (or `flies =`) gives how many flies were loaded into
each vial, left to right — `n = 10, 10, 9`, or a single `n = 10` for every
vial. It overrides the `.cfg` key `flies_per_vial` for that folder and is used by
the [failure-to-climb](#failure-to-climb-ftc) measure. The GUI also applies a
`vials.txt` beside the video when you press *Test parameters*, so the preview
bins vials the same way the batch will.

---

## Output files

Each video's outputs are written next to it, named `<video>.<suffix>`:

| File | Written when | Contents |
|---|---|---|
| `.raw.csv` | always | every detected spot, with filter results and vial |
| `.filtered.csv` | always | kept spots; `y` = height above the vial floor, in pixels |
| `.fng.csv` | `fng_enabled` (default) | one row per FNG event: frames, `drop_cm`, fall duration, `recovery_duration_sec` |
| `.ftc.csv` | `ftc_enabled` (default) | failure to climb per vial |
| `.slopes.csv` | always | climbing velocity per vial (local linear regression) |
| `.diagnostic.png` | always | frames, spots and mean height over time |
| `.tracks.csv` | individual mode | per-fly tracks (`particle` column) |
| `.ftc_particle.csv` | individual mode | one row per fly: outcome, max height, time to the line, falls |
| `.tortuosity_bouts.csv`, `.tortuosity_particle.csv` | individual mode + `tortuosity_enabled` | per-bout and per-fly path metrics |
| `.ROI.png`, `.spot_check.png`, `.processed.png` | `--optimization_plots` or GUI | detection check plots |

The project-level `results.csv`, `fng_results.csv` and `ftc_results.csv` are
described under [Usage](#usage).

---

## Individual-fly tracking mode

By default FreeClimber-FNG runs in **cohort mode**: it analyzes the mean
position of all flies in a vial, exactly as the base FreeClimber platform does.
An optional **individual mode** additionally links per-frame detections into
per-fly trajectories using [TrackPy](http://soft-matter.github.io/trackpy/),
including predictive linking (`trackpy.predict.NearestVelocityPredict`).

Individual mode is opt-in and fully backward compatible — when `analysis_mode`
is unset or `'cohort'`, no linking is done and no per-fly files are written. To enable
it, set `analysis_mode='individual'` in the configuration (`.cfg`) file. When
enabled, a `<video>.tracks.csv` file is written alongside the other outputs
containing the columns `particle, frame, t, vial, x, y` plus the naming-convention
fields. Each vial is linked independently so particle IDs never swap across vials.

The five new configuration keys (defaults shown) are:

| Key | Default | Description |
|---|---|---|
| `analysis_mode` | `'cohort'` | `'cohort'` (mean-position analysis) or `'individual'` (per-fly linking) |
| `link_search_range` | `15` | Maximum inter-frame displacement, in pixels |
| `link_memory` | `3` | Frames a particle may be lost before it is treated as a new track |
| `link_predictor` | `'nearest_velocity'` | `'nearest_velocity'` (predictive) or `'none'` (plain nearest-neighbor) |
| `link_min_track_length` | `5` | Minimum track length, in frames; shorter tracks are dropped |

All five keys are optional. Existing `.cfg` files that omit them keep working
unchanged. See `example/example.cfg` for the keys as commented-out defaults.

In the GUI, individual mode and tortuosity are the two *Step 6: Trajectory
analysis* checkboxes (ticking tortuosity also turns on individual mode).

### Per-fly tortuosity metrics

When individual mode is enabled (and `tortuosity_enabled` is not `False`),
four metrics are computed per fly, per climbing bout:

| Metric | Definition | Range |
|---|---|---|
| `tortuosity` | path length / net displacement | `>= 1` (NaN for closed bouts) |
| `straightness` | net displacement / path length | `[0, 1]` |
| `vertical_efficiency` | upward displacement / path length — the primary climbing metric | `[0, 1]` |
| `mean_turning_angle_rad` | mean absolute turning angle between consecutive step vectors | `[0, pi]` |

A *climbing bout* is a run of consecutive frames in one fly's track whose
vertical velocity exceeds `tortuosity_velocity_threshold` (mm/s), after
Savitzky–Golay smoothing of the track (`tortuosity_smoothing_window`). Bouts
shorter than `tortuosity_bout_min_frames` frames or climbing less than
`tortuosity_bout_min_displacement` mm are dropped. Bout detection is
independent of FNG events. Two files are written:

| File | One row per | Columns |
|---|---|---|
| `<video>.tortuosity_bouts.csv` | `(vial, particle, bout_idx)` | frames, duration, path length and displacements (mm), and the four metrics |
| `<video>.tortuosity_particle.csv` | `(vial, particle)` | `n_bouts`, `median_vertical_efficiency` |

See `scripts/tortuosity.py` for the formal definitions and degenerate-case
handling, and `dashboard/` for a small exploratory plot script.

---

## FNG detection parameters

FNG events are detected on each vial's mean height trace, smoothed and rescaled
to 0–1 so the thresholds are comparable across rigs.

| Key | Default | Description |
|---|---|---|
| `fng_enabled` | `True` | Write `<video>.fng.csv` |
| `fng_smooth_window` | `5` | Rolling-mean window, in frames |
| `fng_climb_thresh` | `0.10` | Minimum rise before a fall, as a fraction of the trace range |
| `fng_fall_thresh` | `0.10` | Minimum drop, as a fraction of the trace range |
| `fng_min_gap` | `5` | Minimum frames between events |
| `fng_min_range_cm` | `2.0` | Smallest range (cm) the trace is rescaled by — see below |
| `fng_recovery_thresh` | — | Accepted for compatibility but **not used**: recovery is the first run of 3 rising frames after the fall |

**`fng_min_range_cm` (new).** Rescaling each vial to its own 0–1 range meant
that in a vial where nobody climbs, a few pixels of jitter were stretched to
the full range and scored as climb-then-fall events. In 30 of the project's real
videos, 33 of 86 FNG events (38%) came from vials whose mean height moved less
than 2 cm in total. The trace is now rescaled by at least `fng_min_range_cm`,
so a fall must be at least `0.10 × 2 cm = 2 mm` on the mean trace. Vials whose
height range exceeds 2 cm are unaffected, and the synthetic validation output is
unchanged. Set `fng_min_range_cm=0` to reproduce the old behavior.

---

## Failure to climb (FTC)

A fly that **fails to climb** never ascends past a height line within a time
limit. This is a different outcome from FNG: an FNG fly climbed and then fell,
an FTC fly never ascended. The measure follows the classic "percentage of flies
past the line" climbing index, and is written for every video in both modes.

Each fly is given exactly one outcome:

| Outcome | Rule |
|---|---|
| `fng` | climbed and fell at least once — **whether or not it reached the line first** (a partial climb that ends in a fall is a fall) |
| `climber` | reached the line, no fall |
| `ftc` | never reached the line and never fell |
| `unscored` | never reached the line or fell, but tracked for less than `ftc_min_coverage` of the window |

**The line** is `ftc_height_cm` above the vial floor when that key is set.
Otherwise it defaults to **the top of the drawn ROI box**, less one spot
`diameter`, because TrackPy cannot detect a fly centred right on the image
edge. Draw the ROI so its top edge is where you want the finish line.

**`<video>.ftc.csv`** — one row per vial, with the naming-convention fields
first. Counts are medians over runs of `ftc_eval_frames` frames:

| Column | Meaning |
|---|---|
| `window_start_frame`, `window_end_frame`, `window_sec` | assessment window (frames relative to `crop_0`) |
| `n_expected` | flies loaded (`flies_per_vial` or `vials.txt` `n =`), else blank |
| `n_detected_start`, `n_detected_end`, `n_detected_max` | flies detected at the start, at the end, and at most at once |
| `n_reached_line` | most flies seen above the line at once during the window |
| `n_above_line_end` | flies above the line at the end of the window (classic index) |
| `ftc_count`, `ftc_fraction` | flies that never reached the line |
| `method` | `expected`: `n_expected − n_reached_line` over `n_expected`; `detected`: `n_detected_max − n_reached_line` over `n_detected_max` |
| `count_warning` | `fewer_detected_than_expected`, `more_above_than_expected`, `no_flies_detected`, `detections_dropped` |
| `n_tracks_climber/_fng/_ftc/_unscored` | per-fly outcome counts (individual mode only) |

"Reached the line" uses the peak number of flies above the line, not the number
above it at the end, so a fly that reached the line and later fell (FNG) or
stopped and was lost from view is not counted as a failure. These per-vial
counts do not use tracking, so in cohort mode a fly that climbed partway and
then fell *is* included in `ftc_count`. Use individual mode, where
`n_tracks_fng` / `n_tracks_ftc` separate those flies.

**`<video>.ftc_particle.csv`** (individual mode) — one row per linked fly:
`first_frame`, `last_frame`, `coverage` (fraction of the window tracked),
`start_height_cm`, `max_height_cm`, `reached_line`, `latency_sec` (time from
the window start to first crossing the line), `n_falls`, `outcome`.
`latency_sec` is blank when the line was not reached, or when the track began
too late to time the climb. Flies that never reach the line are
right-censored at `window_sec`, so latency is suited to survival analysis
(Kaplan–Meier / Cox) as well as to the yes/no FTC fraction.

| Key | Default | Description |
|---|---|---|
| `ftc_enabled` | `True` | Write the FTC files |
| `ftc_height_cm` | `None` | Height of the line above the vial floor, cm; `None` = top of the ROI box (less one `diameter`) |
| `ftc_window_sec` | `None` | Time limit in seconds; `None` = to the end of the cropped video |
| `ftc_start_frame` | `crop_0` | Absolute video frame the time limit starts from (e.g. the tap) |
| `ftc_eval_frames` | `5` | Frames each count is a median over |
| `ftc_min_coverage` | `0.8` | Fraction of the window a non-climbing track must cover to be scored `ftc` |
| `flies_per_vial` | `None` | Flies loaded: an int, or a left-to-right list such as `[10, 10, 9]` |
| `floor_y` | `None` | Vial floor, in pixels from the top of the ROI; `None` = bottom edge of the ROI |
| `background_image` | `None` | Image (or video) of the empty vials to use as the background — see below |

The GUI's *Step 7* row sets the line height (blank = top of the ROI box), time limit, flies per vial and
the vial floor (entered as an image y-coordinate; it is converted to `floor_y`
on save). *Test parameters* then draws the floor (cyan) and the line (magenta
dashes) on the check-frame panel so both can be checked by eye.

### Getting a trustworthy FTC count

1. **Set the floor.** Heights are measured from `floor_y`, or from the bottom
   edge of the ROI. Flies resting on the floor should read about 0 cm. If they
   sit well above 0 in `ftc_particle.csv`, the ROI extends below the floor, so
   set `floor_y` (or the GUI's *Vial floor y*).
2. **Give the fly count.** Flies that never move can be subtracted into the
   background, because the background is the median of the `blank_0`–`blank_n`
   frames. Those are exactly the flies that fail to climb. With `flies_per_vial`
   (or `n =` in `vials.txt`) they are still counted as failures, because the
   count is `n_expected − n_reached_line`.
3. **Use an empty-vial background if you can.** `background_image` takes a
   picture of the same vials, in the same position, without flies. The image
   must have the video's resolution. The background then contains no flies, so
   motionless flies stay visible.
4. **Read `count_warning`.** `detections_dropped` means fewer than half the
   flies seen at the peak are still detected at the end, which usually means
   flies are being lost to background subtraction.

---

## Upgrade notes (September 2026)

The failure-to-climb update also fixed several bugs. These fixes change some
existing outputs, so results from before and after the update should not be
pooled without checking:

- **Heights are measured from the vial floor.** `y` in `*.filtered.csv` and
  `*.tracks.csv` is now height above `floor_y` (default: the ROI bottom). It
  used to be height above the lowest detection in each video, which moved from
  video to video. Differences in `y`, slopes, FNG events and tortuosity are
  unchanged; the `intercept` in `*.slopes.csv` shifts by a constant.
- **File names with 3-letter extensions** (`.mov`, `.mp4`) are no longer
  truncated. Previously the last naming-convention field came out blank or cut
  short (e.g. `vial_ID` `climbing__1`; `rep` blank). `vial_ID` now includes
  that field when `vial_id_vars` covers it.
- **FNG in non-climbing vials** — see `fng_min_range_cm` above.
- **Step 6 no longer aborts the video** when the plotted frame has no
  detections or a vial has no usable detections. In the project's batch logs
  this was the cause of most missing `*.slopes.csv` / `*.diagnostic.png`
  files. Unfittable vials now get a row of blanks in `*.slopes.csv`.
- `window` is kept an integer when the cropped video is shorter than it (it
  used to become a float and crash), `check_frame` is kept inside the cropped
  video, and `blank_0`/`blank_n` are applied relative to `crop_0`.

---

## Repository Structure

| File/Folder | Description |
|---|---|
| `scripts/detector_fng.py` | Core detection, FNG and failure-to-climb logic |
| `scripts/tortuosity.py` | Per-fly tortuosity metrics |
| `scripts/FreeClimber_gui.py` | GUI for calibrating on one video and saving a `.cfg` |
| `scripts/FreeClimber_main.py` | Command-line batch processing |
| `scripts/gather_files.py` | Builds a `.prc` list of videos for `--process_custom` |
| `example/`, `example_other/` | Example videos and configuration files |
| `tests/` | Regression tests and synthetic validation data |
| `dashboard/` | Exploratory plotting scripts for the outputs |
| `paper.md` | JOSS manuscript |
| `TUTORIAL.md` | Usage walkthrough |

---

## Synthetic validation

A synthetic single-fly negative geotaxis video is included in the repository for regression testing against known ground truth. The video was generated with 5 scripted fall events at peak frames 150, 280, 410, 540, and 660 in a 750-frame, 25 fps recording (30 s). Event parameters (rise magnitude, fall distance, recovery duration) are fully specified in the ground-truth CSV.

All synthetic validation assets live under `tests/fixtures/synthetic_validation/`:

| File | Description |
|---|---|
| `freeclimber_fng_validation_video.mp4` | Synthetic test video (750 frames, 25 fps) |
| `freeclimber_fng_validation_video_ground_truth.csv` | Ground-truth fall events (5 events, peak frames 150/280/410/540/660) |
| `freeclimber_fng_validation_video.raw.csv` | Raw TrackPy output used by the regression test |
| `freeclimber_fng_validation_video.fng.csv` | Expected FNG output from the fixed pipeline |
| `freeclimber_fng_test_validation_README.txt` | Full description of video parameters and detection notes |

The regression test (`tests/test_fng_bounds_and_detection.py`) loads the raw CSV, runs the FNG detection pipeline, and asserts that exactly 5 events are detected with peak frames within ±5 frames of ground truth and no impossible (out-of-bounds) frame indices.

**Run the tests:**

```bash
pip install pytest
pytest tests -v
```

No video decoding or FFmpeg is required — the tests use pre-computed CSVs and
synthetic tracks. Besides the FNG regression test, the suite (60 tests) covers
linking, tortuosity, the `vials.txt` sidecar, failure to climb
(`tests/test_ftc.py`) and the September 2026 bug fixes
(`tests/test_bugfixes.py`).

---

## Citing This Work

If you use FreeClimber-FNG in your research, please cite both this tool and the original FreeClimber platform:

**FreeClimber-FNG:**
> Vasu, J. (2026). FreeClimber-FNG (v1.0.4). Zenodo. https://doi.org/10.5281/zenodo.18090957

**Original FreeClimber:**
> Spierer, A. N., Zhuo, L., Zhu, C. T., & Rand, D. M. (2020). FreeClimber: Automated quantification of climbing performance in *Drosophila*. *Journal of Experimental Biology*, 223, jeb229377. https://doi.org/10.1242/jeb.229377

---

## License

This project is licensed under the MIT License, consistent with the original FreeClimber license. Modifications by Jordan Vasu are open-source under the same terms, with attribution to the original authors.

---

## Authors

**FreeClimber-FNG modifications:** Jordan Vasu

**Original FreeClimber:** Adam N. Spierer, Lei Zhuo, and colleagues — Brown University Computational Biology Core
