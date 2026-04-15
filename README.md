# FreeClimber-FNG

FreeClimber-FNG is a Python-based tool for automated detection of **Failed Negative Geotaxis (FNG)** events and calculation of fall distance in *Drosophila melanogaster* locomotor assays. It extends the [FreeClimber](https://github.com/adamspierer/FreeClimber) platform (Spierer et al., 2020) with new functionality for identifying climb-to-fall transitions and quantifying recovery distance.

> **Note:** This repository is a fork of [adamspierer/FreeClimber](https://github.com/adamspierer/FreeClimber). The core particle detection and climbing velocity pipeline is the work of Adam N. Spierer and colleagues. FreeClimber-FNG adds FNG-specific detection on top of that foundation.

---

## What's New in FreeClimber-FNG

The following capabilities were added by Jordan Vasu (2025):

- **FNG detection** — automated identification of climb-to-fall transitions (failed negative geotaxis events) within climbing assay videos
- **Fall distance measurement** — calculates the vertical distance between the position at fall initiation and the position at fall recovery; this metric captures the interval during which a fly falls and subsequently recovers, and can be interpreted as a proxy for fall severity

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
