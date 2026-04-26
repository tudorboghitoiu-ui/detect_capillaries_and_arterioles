# detect_capillaries_and_arterioles

A Fiji/ImageJ macro for automated, compartment-specific quantification of amyloid-β in brain vasculature from confocal Z-stacks.

The macro segments four mutually exclusive compartments — capillary basement membrane, large vessel wall, large vessel smooth-muscle layer (SMA), and non-vascular tissue — and reports per-compartment amyloid load, intensity statistics, and morphometric measurements. It was developed for studies of intramural periarterial drainage (IPAD) and cerebral amyloid angiopathy in mouse hippocampus, but the pipeline applies to other regions and species without modification provided that vessel and nucleus dimensions are comparable to the validation conditions.

**Status:** Version 1.0. Validated on 12 confocal Z-stacks from 5 wild-type mice. A methodology manuscript describing the pipeline, validation, and statistical analysis is in preparation.

---

## Table of contents

- [Quick start](#quick-start)
- [What it does](#what-it-does)
- [Installation](#installation)
- [Input requirements](#input-requirements)
- [Running the macro](#running-the-macro)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Calibration for non-default acquisitions](#calibration-for-non-default-acquisitions)
- [Troubleshooting](#troubleshooting)
- [Validation summary](#validation-summary)
- [Versioning](#versioning)
- [Citation](#citation)
- [License](#license)
- [Authorship and acknowledgements](#authorship-and-acknowledgements)
- [Contact](#contact)

---

## Quick start

If your input is a 4-channel Leica `.lif` Z-stack with channel order C1=SMA, C2=Aβ, C3=ColIV, C4=DAPI, acquired at 20× (NA 0.75, 0.568 µm/pixel), no configuration changes are needed. Otherwise see [Configuration](#configuration) before running.

1. Install Fiji (see [Installation](#installation)).
2. Open Fiji and load your Z-stack via **File > Import > Bio-Formats**.
3. Open `detect_capillaries_and_arterioles.ijm` in the Fiji script editor (**Plugins > New > Macro**, then paste; or **File > Open**).
4. Click **Run**. The macro will prompt once per Fiji session for an output directory (e.g. `~/Desktop/QC/`), then process the active image.
5. Inspect the QC overlay image and the QC summary printed to the Log window.

Per-image processing takes 10–20 seconds on a recent desktop (see [Validation summary](#validation-summary) for indicative timings).

---

## What it does

The macro takes a single confocal Z-stack as input and produces:

- A **QC overlay image** showing the maximum-intensity projection of the original data with the four segmented compartments outlined in colour. This is the primary visual check that the segmentation worked.
- A **per-image QC summary** printed to the Fiji Log, reporting for each compartment: object count, mean / standard deviation / maximum area, and mean / standard deviation pixel intensity in the relevant channel.
- **Per-ROI Results tables** for each compartment, suitable for export to CSV and downstream statistical analysis. Each row corresponds to one segmented object and includes its area, total amyloid pixel area within the object, mean amyloid intensity, particle count, and mean particle size.
- A **non-vascular-tissue amyloid measurement** representing background parenchymal Aβ outside any segmented vessel.
- Optional **binary mask exports** for spatial validation against manual annotation, enabled by setting `CFG_VALIDATION_MODE = 1` in the configuration block.

The pipeline runs in five phases that mirror the figure flowchart in the methodology manuscript: image preparation, channel-by-channel feature extraction, vascular segmentation, compartment classification, and quantification.

---

## Installation

### Required software

- **Fiji** (a distribution of ImageJ bundled with scientific plugins). Download from <https://fiji.sc/>. Tested on Fiji versions distributed in 2025–2026 (ImageJ 1.54 family). Earlier versions may work but are not tested.
- **MorphoLibJ** plugin. Required for morphological filters used in the SMA closing step. Install via **Help > Update... > Manage Update Sites**, then enable the **IJPB-plugins** site and apply changes.
- **Bio-Formats**. Bundled with Fiji and required for loading vendor file formats (Leica `.lif`, Zeiss `.czi`, etc.).

### Macro installation

The macro is a single `.ijm` file. Three installation options:

1. **Run from the script editor.** Open `detect_capillaries_and_arterioles.ijm` in Fiji's script editor and click **Run**. This is the simplest option for occasional use.
2. **Install permanently.** Place the file in `Fiji.app/macros/` and restart Fiji. The macro will appear under **Plugins > Macros**.
3. **Install as a startup macro.** Place the file in `Fiji.app/macros/StartupMacros.fiji.ijm` (concatenated with any existing content). Use this only if you want the macro available without opening the script editor.

No compilation step is required.

---

## Input requirements

| Property | Required value |
|---|---|
| File format | Any format readable by Fiji's Bio-Formats importer (Leica `.lif`, Zeiss `.czi`, OME-TIFF, etc.) |
| Image type | Confocal Z-stack |
| Channels | 3 or 4. Default channel order: C1=SMA (smooth muscle actin), C2=Aβ (amyloid-beta), C3=ColIV (collagen IV, basement membrane), C4=DAPI (nuclei). DAPI is optional. |
| Bit depth | 8-bit or 16-bit grayscale per channel |
| Validation conditions | 1024×1024 pixels, 0.568 µm/pixel, 100–200 Z-slices, mouse hippocampal sections |

The macro processes one image at a time. For batch processing, use Fiji's **Process > Batch > Macro...** or wrap the macro in your own loop over file paths.

If your acquisition differs substantially from the validation conditions — for example a different magnification, a different tissue, or human tissue — see [Calibration for non-default acquisitions](#calibration-for-non-default-acquisitions).

---

## Running the macro

### Single image, interactive

1. Open the Z-stack in Fiji (**File > Import > Bio-Formats**).
2. Open the macro in the script editor and click **Run**.
3. On first run in a Fiji session, a folder picker will appear. Choose where QC outputs should be saved. This choice persists for the rest of the session.
4. Wait 10–20 seconds. Progress is logged to the Log window with timestamps.
5. When the run completes, inspect:
   - The QC overlay image (saved to the chosen output directory and shown on screen).
   - The QC summary in the Log window.
   - The Results tables (one per compartment) shown in the Fiji UI.

### Batch processing

To process every Z-stack in a folder, open Fiji's **Process > Batch > Macro...** dialog, point it at your input folder, and paste the macro into the script box. Set the output directory in the macro configuration (`CFG_QC_DIR`) so it doesn't prompt repeatedly.

For unattended runs that should not prompt for paths, hardcode `CFG_QC_DIR` to an absolute path in Section 1 of the configuration block.

### Validation mode

Setting `CFG_VALIDATION_MODE = 1` in the configuration block enables export of binary masks for each compartment alongside the QC image. These can be loaded into a separate annotation tool to compute spatial overlap statistics (precision, recall, IoU) against manual ground truth. Validation mode roughly doubles processing time and disk usage; leave it disabled for routine analysis.

---

## Configuration

The macro's behaviour is controlled by approximately 110 `CFG_*` variables at the top of the file, organised into ten numbered sections. Most users will only ever change Section 2 (channel assignment).

### Section 2 — channel assignment (the only section most users will touch)

If your acquisition is not a 4-channel Leica `.lif` with the default channel order, edit these five variables:

```javascript
var CFG_CHANNELS_N = 4;          // 3 or 4
var CFG_CHANNEL_SMA = 1;         // input channel for smooth muscle actin
var CFG_CHANNEL_AMYLOID = 2;     // input channel for amyloid-β
var CFG_CHANNEL_COLIV = 3;       // input channel for collagen IV
var CFG_CHANNEL_DAPI = 4;        // input channel for DAPI; ignored if N=3
```

The macro internally remaps these to canonical channel names (C1=SMA, C2=Aβ, C3=ColIV, C4=DAPI) so the rest of the code is independent of how your microscope ordered the channels.

If you have only 3 channels (no DAPI), set `CFG_CHANNELS_N = 3`. The DAPI-bleed-through correction will auto-disable because it requires DAPI.

### When to look at other sections

| Section | When to consider editing |
|---|---|
| 1. FILE I/O | If you want unattended batch processing, hardcode `CFG_QC_DIR` and `CFG_VALIDATION_DIR`. |
| 3. FEATURE TOGGLES | If you want to disable adaptive thresholding, DAPI bleed-through correction, or enable validation export. |
| 4. EMPTY SPACE DETECTION | If your sections have unusually large empty regions (e.g. tissue tears) and the empty-space mask is over- or under-detecting them. |
| 5. SHARED MORPHOLOGICAL CONSTANTS | Almost never. These are tightly tuned defaults. |
| 6. COLLAGEN IV DETECTION | If basement membranes in your data look very different (much thicker, much more diffuse, or substantially different staining intensity). |
| 7. DAPI NUCLEI DETECTION | If nuclei in your data are notably larger or smaller than mouse hippocampal nuclei. |
| 8. SMA DETECTION | If your SMA channel is unusually dim or unusually bright relative to the validation data. The adaptive thresholding handles most of this automatically; manual override is rarely needed. |
| 9. LARGE VESSEL / CAPILLARY CLASSIFICATION | If your tissue contains vessels in size ranges very different from the validation data (e.g. predominantly large arterioles, or predominantly capillaries with no large vessels). |
| 10. AMYLOID DETECTION | If your amyloid staining intensity is substantially different and the threshold of 26 over-segments or under-segments. |

Each section in the macro is preceded by a header comment that explains what the parameters do, what units they are in, and which other parameters they are related to. Read the header before editing any value, and read the warning in Section 6 of the configuration block: changing one parameter without co-adjusting related ones can produce worse results than leaving defaults alone.

---

## Outputs

For each processed image, the macro produces:

### On disk (in `CFG_QC_DIR`)

- `QC_<imagename>.tif` — RGB overlay of the maximum projection with all four compartments outlined.
- If `CFG_VALIDATION_MODE = 1`: `<imagename>_capillaries_mask.tif`, `<imagename>_sma_mask.tif`, `<imagename>_large_vessels_mask.tif`, `<imagename>_amyloid_mask.tif`, `<imagename>_nonvascular_mask.tif`.

### In the Fiji Log window

A summary block per compartment, e.g.:

```
QC SUMMARY: Martor.lif - C10.1.tif
--- CAPILLARIES ---
  Object Count: 157
  Mean Area: 206.36 px2
  StdDev Area: 162.41 px2
  Max Area: 1245.00 px2
  Mean Intensity (C3): 66.58
  StdDev Intensity (C3): 51.51
--- SMA POSITIVE ---
  Object Count: 40
  ...
--- LARGE VESSELS ---
  Object Count: 7
  ...
--- AMYLOID ---
  Object Count: 99
  ...
```

The Log can be saved manually (**File > Save As...** in the Log window) or programmatically with `selectWindow("Log"); saveAs("Text", "/path/to/log.txt");`.

### In the Fiji Results tables

One Results table per compartment with one row per segmented ROI. Columns include the image name, ROI index, ROI area, total amyloid pixel area within the ROI, mean amyloid intensity, particle count, mean particle size, and integrated density. Save each table to CSV via **File > Save As... > CSV** in the Results window.

### Suggested downstream workflow

For statistical analysis, save the per-compartment Results tables to CSV and load them in R or Python. The hierarchical structure of the data (ROIs nested within images nested within mice) requires linear mixed-effects models if you want population-level inferences; per-image summary statistics from the QC log can be used for exploratory plots but are not appropriate as the unit of analysis for hypothesis tests.

---

## Calibration for non-default acquisitions

The macro defaults are tuned for 0.568 µm/pixel imaging of mouse hippocampus. If your acquisition is substantially different, the size-based thresholds will be wrong. Two situations are common.

### Different pixel size

Most size thresholds in the macro are expressed in **pixels squared** (`CFG_COLLAGEN_MIN_SIZE`, `CFG_SMA_MIN_SIZE`, `CFG_LARGE_VESSEL_MIN_SIZE`, `CFG_NUCLEI_MIN_SIZE`, the bled-nuclei size range, the size bin breakpoints in collagen filtering, etc.). If your pixel size differs from the validation, scale these by the ratio of areas:

```
scaling_factor = (validation_pixel_size / your_pixel_size)^2
                = (0.568 µm/px / your_pixel_size)^2
```

**Worked example.** Suppose your acquisition is 0.4 µm/pixel (a 28% finer sampling than validation). Then:

```
scaling_factor = (0.568 / 0.4)^2 = 2.02
```

Multiply every pixel-area threshold by 2.02. For example, `CFG_LARGE_VESSEL_MIN_SIZE = 150` becomes `150 × 2.02 ≈ 303`. The collagen size bins (`<=150`, `<=300`, `<=500`, `<=1000`, `>1000`) become approximately (`<=303`, `<=606`, `<=1010`, `<=2020`, `>2020`). The Phansalkar radii in Section 5 are also pixel-based and should scale by the linear ratio (`0.568 / your_pixel_size`), not the squared ratio.

Do not rescale Yen-factor multipliers, percentile thresholds, or the SMA-collagen overlap level cutoffs — these are dimensionless and independent of pixel size.

### Different tissue or species

If you are imaging human tissue, a different brain region, or a non-mouse species, vessel and nucleus sizes can differ substantially from mouse hippocampus. In this case, pixel-size scaling alone is not sufficient. The recommended workflow is:

1. Manually annotate compartments in 3–5 representative images using ImageJ's ROI Manager or an external annotation tool.
2. Run the macro on those images with the default parameters and compare visually using the QC overlay.
3. If specific compartments are over- or under-detected, adjust the relevant section's parameters guided by what you observe (e.g. if many small capillaries are missed, lower `CFG_COLLAGEN_MIN_SIZE`; if SMA is over-detected on dim vessels, raise the adaptive threshold floor in Section 8).
4. Re-run and re-evaluate. Iterate until the QC overlay matches your annotations on the calibration set.
5. Apply the calibrated parameters to the full dataset.

This calibration process typically takes a few hours and produces parameter changes much smaller than the worked example above. Most parameters are robust across moderate biological variation.

---

## Troubleshooting

### The macro prompts for an output directory every time I run it

Either you are starting a fresh Fiji session each time (the prompt is once per session by design), or you have left `CFG_QC_DIR = ""` and want to suppress the prompt entirely. Set `CFG_QC_DIR` to an absolute path in Section 1 of the configuration block.

### "Channel count mismatch" error

The macro checks that the input image has the number of channels declared in `CFG_CHANNELS_N`. If you set `CFG_CHANNELS_N = 4` but your image is 3-channel, the macro stops with this error. Check your input and update `CFG_CHANNELS_N` accordingly.

### QC overlay shows compartments in the wrong colours

The macro assigns colours by canonical channel role (SMA, Aβ, ColIV, DAPI), not by input channel index. If colours look swapped, your channel assignment in Section 2 is probably wrong. Verify by inspecting your input image's channel order in **Image > Properties** or **Image > Stacks > Tools > Make Substack...**.

### Capillary mask is empty or near-empty

Most often caused by very dim ColIV staining or by the bled-through-nuclei filter removing legitimate vessels. Inspect the Log window for the line `- Detected N bleed-through nuclei regions`. If N is very high (hundreds out of thousands of ROIs), DAPI bleed-through is severe. Either acquire cleaner images or set `CFG_BLEEDTHROUGH_CORRECTION = 0` (which has its own trade-offs — read Section 3 of the configuration block).

### SMA mask captures cellular cytoplasm, not just vessel walls

The SMA antibody is occasionally non-specific. The macro filters SMA detections by spatial overlap with collagen IV (Section 9) and by intensity (Section 8, adaptive). If non-vascular cells are still being detected, raise the L1 / L2 intensity floor in Section 8.

### Run takes much longer than 20 seconds per image

Check the Log timestamps. Empty-space detection (4 median-filtered MAX projections on the full Z-stack) is normally the slowest step at 5–9 seconds. If the run is dominated by SMA classification or the Phansalkar steps, the input image may have many more candidate ROIs than typical (e.g. very dense vasculature or noisy staining). This is usually a property of the data, not a macro bug.

### Java errors in the Fiji editor

Errors of the form `ImageJMacroTokenMaker: could not match input` or `ArrayIndexOutOfBoundsException` in the editor's paint loop are Fiji script-editor UI issues, not macro execution errors. They do not affect the run. Ignore them.

### Reported counts in the QC log differ from the per-ROI Results table

This was a real issue in pre-1.0 versions and was fixed before release. If you observe it in version 1.0, please file a bug report (see [Contact](#contact)) with the input image and the full Log output.

---

## Validation summary

Version 1.0 was validated on **12 confocal Z-stacks from 5 wild-type mice**, acquired on a Leica TCS SP8 confocal microscope at 20× (NA 0.75, 0.568 µm/pixel), 1024×1024 pixels, 116–192 Z-slices per stack. Total processing time across the validation set was approximately 3.5 minutes.

Validation against blinded manual annotation (performed by an independent postdoc on a subset of the dataset) is reported in the methodology manuscript. The validation includes per-compartment precision and recall, intra-rater reliability, and qualitative assessment of edge cases (lumen holes, non-specific SMA staining, regions near tissue tears).

The macro is not currently validated on:

- Human tissue
- Brain regions other than hippocampus
- Magnifications other than 20× / 0.568 µm-per-pixel
- Acquisition modalities other than confocal microscopy (e.g. wide-field, multi-photon, light-sheet)

If you apply it to data outside these conditions, please follow the calibration workflow in [Calibration for non-default acquisitions](#calibration-for-non-default-acquisitions) and report your experience.

---

## Versioning

This macro follows [semantic versioning](https://semver.org/). The version number lives in two places:

- The `Version:` line in the file header of `detect_capillaries_and_arterioles.ijm`.
- The `(v1.0)` tag printed at the start of every run in the Fiji Log.

Each tagged release on GitHub is independently archived on Zenodo with its own DOI (see [Citation](#citation)). The concept DOI always resolves to the latest release.

Patch releases (1.0.x) fix bugs without changing default behaviour. Minor releases (1.x.0) may add new optional features or non-breaking parameter additions. Major releases (x.0.0) may change default behaviour or break existing configurations and will be accompanied by a migration note.

---

## Citation

If you use this macro in published work, please cite both the methodology manuscript and the archived code release.

**Methodology manuscript** (in preparation; placeholder — update when published):

> Boghițoiu T-G, et al. *<TODO_FILL_IN_MANUSCRIPT_TITLE>*. <TODO_FILL_IN_JOURNAL>, 2026. DOI: <TODO_FILL_IN_MANUSCRIPT_DOI>

**Code archive**:

> Boghițoiu T-G. *detect_capillaries_and_arterioles: Compartment-specific quantification of amyloid-β in brain vasculature from confocal Z-stacks* (Version 1.0.0). Zenodo, 2026. DOI: <TODO_FILL_IN_ZENODO_CONCEPT_DOI>

A `CITATION.cff` file is included in the repository; GitHub renders a "Cite this repository" button in the sidebar that uses this metadata.

---

## License

MIT License. See the [LICENSE](LICENSE) file in this repository for the full text.

The MIT license permits commercial and non-commercial use, modification, and redistribution, subject to inclusion of the original copyright notice and disclaimer of warranty.

---

## Authorship and acknowledgements

**Code author**

- Tudor-Gabriel Boghițoiu, Doctoral School of Medicine and Pharmacy, George Emil Palade University of Medicine, Pharmacy, Science, and Technology of Târgu Mureș, Romania

**Collaborators**

The methodology manuscript reports the validation work and full author list, including collaborators who contributed scientific direction, manual annotation, and experimental data. The macro itself was written by the code author listed above; collaborators are acknowledged here for their contributions to the broader project but are not authors of the source code.

**Funding**

Development of this macro was supported by the Romanian National Recovery and Resilience Plan (PNRR), Contract Number CF 63/14.11.2022.

**Data**

The validation dataset consists of confocal Z-stacks of mouse hippocampal sections imaged on a Leica TCS SP8. Acquisition and tissue preparation details are reported in the methodology manuscript.

---

## Contact

For bug reports, feature requests, or questions: please open an issue on the GitHub repository at <https://github.com/tudorboghitoiu-ui/detect_capillaries_and_arterioles>.

For correspondence regarding the methodology, the validation, or scientific applications: contact details are in the methodology manuscript.
