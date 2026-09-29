# FreeClimber–FNG Adaptation

This repository contains modifications to the [FreeClimber](https://github.com/adamspierer/FreeClimber) software (Spierer et al., 2020) to enable automated detection of failed negative geotaxis (FNG) and calculation of fall distance in *Drosophila melanogaster*.

## Original Authors
- Adam N. Spierer, Lei Zhuo, and colleagues  
- Citation: Spierer, A. N. et al. (2020). *Journal of Experimental Biology*, 223, jeb229377.

## Modifications by
- Jordan Vasu (2025–2026)  
- Repository: https://github.com/jordanvasu/FreeClimber-FNG
- Added functions for:
  - FNG detection (climb → fall events)
  - Fall measurements -> for each fall, the distance fallen (`drop_cm`, from the top of the climb to the bottom of the fall), the fall duration, and the recovery time until the _Drosophila_ starts to climb again (`recovery_duration_sec`).
  - Individual-fly tracking and per-fly tortuosity metrics.
  - Failure to climb (FTC) -> a fly that makes no climbing movement during the clip, judged per fly from its own track. This is a separate outcome from FNG: it is a failure to ascend, not a fall; a fly that climbs any distance and then falls counts as a fall. See the "Failure to climb (FTC)" section of README.md.
  

## License
The original FreeClimber is released under the MIT License.  
All modifications herein remain open-source under the same license, with attribution to the original authors.
