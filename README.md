# Human Presence Monitoring Safety System using Phased Array Radar

**Authors:** Eldrin Ereqi & Justin Lee  
**Acknowledgments:** Dr. Arno Thielens, Dr. Kevin Gu, and the IEEE Aerospace and Electronic Systems Society  

## Overview
As dangerous robotics become increasingly capable, they are frequently deployed in cluttered, dynamic environments alongside human workers. This project aims to safely integrate these robotics into shared workspaces by developing a **Human Presence Monitoring Safety System** utilizing the Analog Devices ADALM PHASER kit (CN0566). 

The system continuously monitors a 6-meter radius. When a target is detected, it steers the radar beam to collect micro-Doppler signatures, classifies whether the target is a human or a non-human object using Machine Learning, and triggers a 3-Zone Safety Alert System (Warn, Slow Down, Stop) based on proximity.

## System Pipeline
Our real-time processing pipeline consists of several key stages:

1. **Background Calibration:** The system initializes by capturing a range-angle heatmap of the empty environment to establish a baseline.
2. **Target Detection (Search Mode):** The radar continuously sweeps between -60° and +60° in 5° increments. By subtracting the live heatmap from the background baseline, new targets are easily identified.
3. **Beamsteering (Track Mode):** Once a significant target is detected, the radar steers its beam directly at the target to focus on its specific movements.
4. **Signal Processing:** 
   * **MTI Filter:** A Moving Target Indicator filter is applied to isolate moving targets from stationary clutter.
   * **CFAR Filter:** A Constant False Alarm Rate filter is used to distinguish the moving target from surrounding noise.
   * **DBSCAN:** Applied to cluster detections and eliminate radar distortion and artifacts (like radar spearing).
5. **Classification:** Features are extracted from the micro-Doppler signatures and fed into a Support Vector Machine (SVM). The SVM was trained on a custom dataset of 1,600 micro-Doppler heatmaps (covering humans, metal carts with boxes, and ambient noise across various lab/hallway environments and angles).
6. **Safety Zones:** If the target is classified as human, their distance dictates the safety protocol: Warn (Zone 1), Slow Down (Zone 2), or Stop (Zone 3).

## Project Structure
This repository combines our custom DSP and safety logic with the foundational hardware control scripts provided by the Analog Devices (ADI) Phaser tutorial.

* **`main_script.m`** - The core execution loop combining the search, track, process, and UI state machines.
* **`Custom_Functions/`** - Our custom modularized logic for this project:
  * `capture_data.m`, `process_signals.m` - Hardware interfacing and signal processing (MTI, CFAR).
  * `extract_targets.m` - Target isolation using DBSCAN.
  * `steer.m`, `precompute_phases.m` - Beamsteering logic.
  * `setup_ui.m`, `update_dashboard.m`, `update_state.m` - 3-Zone Safety UI and state machine logic.
* **`HelperFunctions/`** - Original hardware configuration, calibration, and control scripts provided by the ADI MATLAB tutorial (e.g., `setupLabRadar.m`, `setupPluto.m`, calibration weights).

## Hardware Setup & Requirements
This project requires the **Analog Devices ADALM PHASER kit (CN0566)**, configured as an FMCW radar. 

**Required Connections:**
* USB-C Power cable inserted into the antenna board.
* USB to Micro-USB cable connecting the computer to the ADI Pluto.
* USB to Ethernet cable connecting the computer to the Raspberry Pi.

*Note: For detailed environment setup, refer to the [MATLAB Phaser Setup Guide](https://wiki.analog.com/resources/eval/user-guides/phaser).*

## Known Limitations & Future Work
While the individual components of this pipeline are successful, full seamless integration requires the ability to dynamically change the number of samples captured by the device during runtime. Currently, hardware buffer constraints make it difficult to switch rapidly between a fast, low-sample sweep for the range-angle heatmap and a high-sample collection for accurate micro-Doppler tracking. 

## Credits & License
* **Base Codebase:** A significant portion of the hardware interfacing and calibration code (`HelperFunctions/`) was created by Analog Devices and provided via their Phaser tutorial repository. 
* **Custom Logic:** The safety pipeline, DBSCAN integration, SVM classification logic, and modular architecture were developed by our team for the IEEE Challenge.
