# FreeClimber-FNG

FreeClimber-FNG is a Python-based tool for automated detection of **Failed Negative Geotaxis (FNG)** events and calculation of fall distance in *Drosophila melanogaster* locomotor assays. It extends the [FreeClimber](https://github.com/adamspierer/FreeClimber) platform (Spierer et al., 2020) with new functionality for identifying climb-to-fall transitions and quantifying recovery distance.

> **Note:** This repository is a fork of [adamspierer/FreeClimber](https://github.com/adamspierer/FreeClimber). The core particle detection and climbing velocity pipeline is the work of Adam N. Spierer and colleagues. FreeClimber-FNG adds FNG-specific detection on top of that foundation.

---

## What's New in FreeClimber-FNG

The following capabilities were added by Jordan Vasu (2025):

- **FNG detection** — automated identification of climb-to-fall transitions (failed negative geotaxis events) within climbing assay videos
- **Fall distance measurement** — calculates the vertical distance between the position at fall initiation and the position at fall recovery; this metric captures the interval during which a fly falls and subsequently recovers, and can be interpreted as a proxy for fall severity
- **Individual-fly tracking mode** — an optional mode that links per-frame detections into per-fly trajectories using [TrackPy](http://soft-matter.github.io/trackpy/) (including predictive linking), in addition to the default cohort (mean-position) analysis
- **Per-fly tortuosity / meandering metrics** — when individual mode is enabled, computes path tortuosity, straightness, and mean turning angle per fly per climbing bout, quantifying how directly (or erratically) each fly climbs

These additions are implemented in `detector_fng.py` and are designed to integrate with the existing FreeClimber parameter configuration and batch processing workflow.

---

## Installation

FreeClimber-FNG uses the same environment and dependencies as the base FreeClimber platform. We recommend running in an Anaconda virtual environment.

**1. Create and activate a Python 3.6 environment:**

```bash
conda create -n freeclimber python=3.6 anaconda
conda activate freeclimber
```

**2. Install dependencies:**

```bash
pip install FreeClimber
```

**3. Clone this repository:**

```bash
git clone https://github.com/jordanvasu/FreeClimber-FNG.git
cd FreeClimber-FNG
```

For full dependency details (FFmpeg, wxPython, trackpy, etc.), see the [FreeClimber installation guide](https://github.com/adamspierer/FreeClimber#installing).

---

## Usage

FreeClimber-FNG can be run via GUI or command line, following the same conventions as the base platform. See [TUTORIAL.md](TUTORIAL.md) for step-by-step instructions on parameter configuration and batch processing (covers base FreeClimber workflow). FNG-specific detection is handled by `detector_fng.py` and is described in the [paper](paper.md).

**GUI:**

```bash
pythonw ./scripts/FreeClimber_gui.py --video_file ./example/<video_file.suffix>
```

**Command line:**

```bash
python FreeClimber_main.py --config_file ./example/example.cfg
```

### Per-folder vial count (`vials.txt`)

A single shared `.cfg` can be reused across folders whose videos have different
vial counts. Drop a plain-text file named `vials.txt` into a folder containing
just the integer count (e.g. `3`), and it overrides the `vials` value in the
`.cfg` for every video in that folder. Accepted contents: the first non-blank,
non-`#`-comment line as a bare integer (`3`) or `vials=3` / `vials: 3`. A
missing or unparseable file leaves the `.cfg` value unchanged, so batch runs
never break.

---

## Individual-fly tracking mode

By default FreeClimber-FNG runs in **cohort mode**: it analyzes the mean
position of all flies in a vial, exactly as the base FreeClimber platform does.
An optional **individual mode** additionally links per-frame detections into
per-fly trajectories using [TrackPy](http://soft-matter.github.io/trackpy/),
including predictive linking (`trackpy.predict.NearestVelocityPredict`).

Individual mode is opt-in and fully backward compatible — when `analysis_mode`
is unset or `'cohort'`, output is byte-identical to previous behavior. To enable
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

> **Note:** Individual mode is currently command-line only. GUI exposure is
> deferred until linking is validated on real multi-fly data.

### Per-fly tortuosity metrics

When individual mode is enabled, three tortuosity metrics are additionally
computed per fly, per climbing bout, and written to `<video>.tortuosity.csv`:

| Metric | Definition | Range |
|---|---|---|
| `tortuosity` | path length / net displacement | `>= 1` (NaN for closed bouts) |
| `straightness` | net displacement / path length | `[0, 1]` |
| `mean_turning_angle_rad` | mean absolute turning angle between consecutive step vectors | `[0, pi]` |

A *climbing bout* is the contiguous frame window from the previous FNG event's
fall-end (or the start of the recording, for the first event) up to and
including the current event's peak frame. Bouts are defined per vial using
the same FNG events emitted to `<video>.fng.csv`, and metrics are then
computed per linked particle inside each window. One row is emitted per
`(vial, event_idx, particle)` tuple. See `scripts/tortuosity.py` for the
formal definitions and degenerate-case handling, and `dashboard/` for a
small exploratory plot script.

---

## Repository Structure

| File/Folder | Description |
|---|---|
| `detector_fng.py` | Core FNG detection logic (FreeClimber-FNG additions) |
| `scripts/` | GUI and command line interface wrappers |
| `example/` | Example video and configuration files |
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

**Run the regression test:**

```bash
pip install pytest numpy pandas scipy
pytest tests/test_fng_bounds_and_detection.py -v
```

No video decoding or FFmpeg is required — the test operates on the pre-computed raw CSV.

---

## Citing This Work

If you use FreeClimber-FNG in your research, please cite both this tool and the original FreeClimber platform:

**FreeClimber-FNG:**
> Vasu, J. (2025). FreeClimber-FNG: Automated detection of failed negative geotaxis and fall distance in *Drosophila* climbing assays. *Journal of Open Source Software* (under review).

**Original FreeClimber:**
> Spierer, A. N., Zhuo, L., Zhu, C. T., & Rand, D. M. (2020). FreeClimber: Automated quantification of climbing performance in *Drosophila*. *Journal of Experimental Biology*, 223, jeb229377. https://doi.org/10.1242/jeb.229377

---

## License

This project is licensed under the MIT License, consistent with the original FreeClimber license. Modifications by Jordan Vasu are open-source under the same terms, with attribution to the original authors.

---

## Authors

**FreeClimber-FNG modifications:** Jordan Vasu

**Original FreeClimber:** Adam N. Spierer, Lei Zhuo, and colleagues — Brown University Computational Biology Core
