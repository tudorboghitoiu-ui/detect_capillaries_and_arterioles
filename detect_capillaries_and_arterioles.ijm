//============================================================================
//  detect_capillaries_and_arterioles.ijm
//
//  ImageJ/Fiji macro for automated compartment-specific quantification of
//  amyloid-beta in brain vasculature from confocal Z-stacks. Segments four
//  mutually exclusive vascular compartments (capillary basement membrane,
//  large vessel wall, large vessel SMA, non-vascular tissue) and reports
//  per-compartment amyloid load, intensity statistics, and morphometric
//  measurements.
//
//  Validated on mouse hippocampal sections imaged at 20× (NA 0.75,
//  0.568 µm/pixel). The pipeline is applicable to other brain regions
//  without modification, and to other species or magnifications provided
//  vessel and nucleus dimensions are comparable; substantially different
//  acquisition parameters will require recalibration of the configuration
//  block (see README.md and Section 6 of the configuration block).
//
//  Version: 1.0
//  Author:  Tudor-Gabriel Boghițoiu
//           Doctoral School of Medicine and Pharmacy,
//           George Emil Palade University of Medicine, Pharmacy, Science,
//           and Technology of Târgu Mureș, Târgu Mureș, Romania
//
//  License: MIT (see LICENSE file in the repository)
//
//  Copyright (c) 2026 Tudor-Gabriel Boghițoiu
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction. The full
//  license text is provided in the LICENSE file accompanying this macro.
//
//----------------------------------------------------------------------------
//  INPUT REQUIREMENTS
//----------------------------------------------------------------------------
//  - Confocal Z-stack in a Fiji-readable format (Leica .lif tested;
//    Bio-Formats import recommended for other vendors).
//  - 3 or 4 channels, with one channel each for: alpha-smooth muscle actin
//    (SMA), amyloid-beta (Aβ), collagen IV (ColIV), and optionally DAPI.
//  - Channel-to-marker assignment is configured via CFG_CHANNEL_* variables
//    in Section 2 of the configuration block.
//  - Required Fiji plugins: MorphoLibJ (for morphological filters) and
//    Bio-Formats (for vendor file import). Both are bundled with recent
//    Fiji distributions.
//
//----------------------------------------------------------------------------
//  QUICK START
//----------------------------------------------------------------------------
//  1. Open the macro in Fiji's script editor (Plugins > New > Macro,
//     then paste, or File > Open).
//  2. If your input is not a 4-channel Leica acquisition with the default
//     channel order (C1=SMA, C2=Aβ, C3=ColIV, C4=DAPI), edit CFG_CHANNELS_N
//     and the four CFG_CHANNEL_* variables in Section 2 of the
//     configuration block to match your data.
//  3. Click Run. The macro will prompt for an output directory the first
//     time it runs in a given Fiji session, then process the active image
//     and save QC outputs and per-ROI measurements to that directory.
//
//  For full documentation, parameter reference, calibration guide, and
//  troubleshooting, see README.md in the repository.
//
//----------------------------------------------------------------------------
//  CITATION
//----------------------------------------------------------------------------
//  If you use this macro in published work, please cite:
//    Boghițoiu T-G et al. (manuscript in preparation). [Citation details
//    will be updated upon publication.]
//============================================================================

setBatchMode(true);
//============================================================================
//=================== CONFIGURATION BLOCK ====================================
//============================================================================
// All tunable parameters are grouped below by processing stage rather than by
// parameter type, so that all settings relevant to (e.g.) collagen IV
// detection live in a single section. Edit values in place; do not move
// declarations outside this block without understanding dependencies.
// Runtime state variables (computed during processing) are declared in a
// separate section AFTER the CFG block.
//
// RECALIBRATION GUIDANCE FOR NEW DATASETS
// ----------------------------------------------------------------------------
// The only parameter that typically requires recalibration for a new dataset
// is CFG_AMYLOID_THRESHOLD (Section 10). It depends directly on laser power,
// detector gain, and antibody titer for the amyloid channel, all of which
// change between acquisitions.
//
// All other parameters were tuned iteratively on the validation dataset and
// are interdependent. Changing one value (e.g., a Phansalkar k1) without
// co-adjusting related values (the corresponding radius, the relevant
// min-size threshold, the Yen factors downstream) can produce worse results
// than leaving defaults alone. For adaptation to acquisitions with
// substantially different optics, staining protocols, or tissue preparation,
// compare macro output against manual annotation on a small set of
// representative images rather than tuning individual values in isolation.
// Consult the calibration section of the README for worked examples.
//
// Section index:
//   1.  FILE I/O
//   2.  INPUT IMAGE
//   3.  FEATURE TOGGLES
//   4.  EMPTY SPACE DETECTION
//   5.  SHARED MORPHOLOGICAL CONSTANTS
//   6.  COLLAGEN IV DETECTION
//   7.  DAPI NUCLEI DETECTION
//   8.  SMA DETECTION
//   9.  LARGE VESSEL / CAPILLARY CLASSIFICATION
//  10.  AMYLOID DETECTION

//----------------------------------------------------------------------------
// 1. FILE I/O
//----------------------------------------------------------------------------
// Controls where the macro writes its QC images and (optionally) validation
// export masks. Both paths support three modes of specification:
//   (a) leave blank "" to have the macro prompt once per Fiji session,
//   (b) hardcode an absolute path here for unattended batch processing,
//   (c) rely on the prompt cached from a previous call in the same session.
// Choice (a) is recommended for interactive use; (b) for scripted batch runs.
// Session caching uses JVM system properties, which are cleared automatically
// when Fiji is closed — choices do NOT persist across Fiji launches.

var CFG_QC_OUTPUT_DIR = "";

// VALIDATION MODE produces two additional outputs beyond the routine QC image:
// (1) binary compartment masks (capillaries / SMA / large vessels / amyloid)
//     as TIFFs, intended for downstream spatial analysis in R;
// (2) a calibrated diagnostic composite (DIAG file) with masks overlaid on
//     the raw MAX projection, intended for loading external .roi files
//     (e.g., manual expert annotations) for pixel-accurate validation.
// Leave false for routine analysis. Enable only when validating the macro
// against manual ground truth or when downstream spatial statistics are
// needed. Extra I/O adds a small per-image overhead.
var CFG_VALIDATION_MODE = false;
var CFG_VALIDATION_EXPORT_DIR = "";

//----------------------------------------------------------------------------
// 2. INPUT IMAGE
//----------------------------------------------------------------------------
// Describes the structure of the input confocal stack: how many channels,
// and which input channel corresponds to which biological marker.

// Number of input channels in the image. Supported values: 3 or 4.
//   - When 4: the macro operates as originally designed, with DAPI-based
//     bleed-through correction (the standard Leica SP8 four-channel protocol).
//   - When 3: DAPI processing is skipped entirely and CFG_BLEEDTHROUGH_CORRECTION
//     is force-disabled (no DAPI channel available to calibrate it). The
//     macro falls back to Yen intensity thresholding for XS-round collagen
//     validation — a path that, in the validation dataset, produced results
//     very similar to the bleed-through-corrected path.
//   - The channel count of each input image is verified strictly at macro
//     start; mismatch with CFG_CHANNELS_N causes the macro to exit with a
//     clear error rather than silently drop or mis-interpret channels.
var CFG_CHANNELS_N = 4;

// Channel-to-marker assignment. Tells the macro which input channel carries
// which biological marker. The defaults below match the standard Leica SP8
// four-channel protocol used in the validation dataset:
//   C1 = SMA, C2 = Amyloid, C3 = ColIV, C4 = DAPI
// For datasets acquired with a different channel order (e.g., amyloid on
// channel 1), edit ONLY these four values. The macro will rename channels
// internally after splitting so that all downstream code continues to refer
// to "C1 = SMA", "C2 = Amyloid", etc.
// Values must be unique integers in 1..CFG_CHANNELS_N. Duplicate assignments
// or out-of-range values are caught at macro start with a clear error.
var CFG_CHANNEL_SMA     = 1;
var CFG_CHANNEL_AMYLOID = 2;
var CFG_CHANNEL_COLIV   = 3;
var CFG_CHANNEL_DAPI    = 4;  // Ignored when CFG_CHANNELS_N = 3

//----------------------------------------------------------------------------
// 3. FEATURE TOGGLES
//----------------------------------------------------------------------------
// High-level on/off switches for three behaviours that change which code
// paths execute. Unlike numeric parameters elsewhere in the CFG block, these
// flags select between qualitatively different pipeline strategies.

// CFG_BLEEDTHROUGH_CORRECTION
//   - When true (default): the macro segments DAPI-positive nuclei and uses
//     them to identify nuclear-coincident collagen IV ROIs as artefacts. The
//     90th percentile of collagen intensity within these artefacts is used
//     to calibrate a noise-floor threshold for the collagen channel, and the
//     XS round category of the collagen ROI filter uses spatial overlap
//     with the intensity-cleaned mask.
//   - When false: DAPI segmentation is skipped. The collagen mask is built
//     from auto_local_collagen4 via morphological closing only, and the XS
//     round category falls back to Yen intensity threshold validation.
//   - Force-disabled internally whenever CFG_CHANNELS_N = 3 (no DAPI to
//     calibrate against).
//   - Recommended: leave true for datasets with visible nuclear-coincident
//     ColIV artefacts (DAPI emission tail and/or perinuclear lipofuscin
//     autofluorescence). Disabling is supported for 4-channel data but has
//     not been systematically compared against the enabled path on the
//     validation dataset, so the quantitative impact of disabling on
//     otherwise-comparable inputs is not known. The fallback path is the
//     same one that automatically executes on the 2 of 12 validation images
//     with fewer than 30 bled nuclei, so the code is exercised, but it is
//     not a direct substitute.
var CFG_BLEEDTHROUGH_CORRECTION = true;

// CFG_ADAPTIVE_SMA_ENABLED
//   - When true (default): per-overlap-level SMA intensity factors are
//     computed dynamically from the mean C1 intensity of the current image.
//     Dimmer images get more permissive thresholds, brighter images get
//     stricter ones. This corrects for acquisition-to-acquisition brightness
//     variation without requiring per-image parameter tuning.
//   - When false: the fixed CFG_SMA_THRESHOLD_FACTOR (Section 8) is used for
//     all overlap levels, identical to the original macro logic.
//   - Recommended: leave true. The adaptive path was developed specifically
//     to handle the variance observed across the validation dataset.
var CFG_ADAPTIVE_SMA_ENABLED = true;

// CFG_ADAPTIVE_EMPTY_SPACES_ENABLED
//   - When true (default): empty-space detection uses per-channel minimum
//     intensity + tolerance, then AND-combines the per-channel masks. Robust
//     to acquisition differences because each channel is calibrated to its
//     own noise floor.
//   - When false: normalises all channels, sums them, and thresholds the
//     combined intensity against CFG_EMPTY_SPACE_INTENSITY_MAX (static mode).
//   - Recommended: leave true. The static path was the original implementation
//     and is retained for reproducibility of older analyses, not because it
//     is expected to produce better results on new data.
var CFG_ADAPTIVE_EMPTY_SPACES_ENABLED = true;

//----------------------------------------------------------------------------
// 4. EMPTY SPACE DETECTION
//----------------------------------------------------------------------------
// Identifies regions of the image that contain no biological tissue —
// tissue holes, mounting bubbles, ventricular space, and empty coverslip
// area outside the section. Used downstream to (a) exclude empty regions
// from mean-intensity calculations, and (b) optionally define a tissue mask
// for the non-vascular compartment analysis.
//
// Two modes are available, selected by CFG_ADAPTIVE_EMPTY_SPACES_ENABLED
// (Section 3). Only one mode's parameters are read per run.

// Minimum ROI size for empty spaces (px²). Small disconnected dark regions
// are almost always tissue artefacts (lumen, cell boundaries) rather than
// true empty regions; filtering below this size prevents those false
// positives. At the validation calibration (0.568 µm/pixel), 500 px²
// corresponds to roughly 161 µm² of image area.
// This default has not been systematically optimised on the validation
// dataset and may warrant revisiting for datasets where empty space
// detection is used for quantitative analysis (e.g., if the non-vascular
// compartment is a primary outcome). For fields of view that contain small
// intentional empty regions, lower it.
var CFG_EMPTY_SPACE_MIN_SIZE = 500;

// Median filter radius applied to each channel stack before Z-projection.
// Suppresses per-voxel shot noise so that "empty" regions read as genuinely
// near-zero rather than as sparse noise pixels. Larger values give cleaner
// empty-space masks at the cost of slightly dilating tissue-edge detection.
// Default 5 was used across the validation dataset but has not been
// systematically compared against other values.
var CFG_EMPTY_SPACE_MEDIAN_RADIUS = 5;

// --- ADAPTIVE MODE PARAMETERS (used when CFG_ADAPTIVE_EMPTY_SPACES_ENABLED = true) ---
// For each channel, the macro finds the minimum pixel intensity in the
// MAX projection (representing the darkest region of the image), adds the
// per-channel tolerance below, and thresholds all pixels at-or-below this
// value as candidate empty-space. The final empty-space mask is the AND of
// all per-channel candidate masks — a pixel is "empty" only if ALL channels
// read as near-noise at that location.
//
// Higher tolerance = more aggressive empty-space detection (more pixels
// classified as empty). Lower tolerance = more conservative. Defaults of 3
// were used across the validation dataset but have not been systematically
// optimised; raising or lowering them uniformly is a reasonable first
// calibration step for acquisitions with substantially different noise
// floors than the validation hardware.
var CFG_EMPTY_SPACES_TOLERANCE_C1 = 3;   // SMA channel
var CFG_EMPTY_SPACES_TOLERANCE_C2 = 3;   // Amyloid channel
var CFG_EMPTY_SPACES_TOLERANCE_C3 = 3;   // ColIV channel
var CFG_EMPTY_SPACES_TOLERANCE_C4 = 3;   // DAPI channel (only used when CFG_CHANNELS_N = 4)

// --- STATIC MODE PARAMETER (used when CFG_ADAPTIVE_EMPTY_SPACES_ENABLED = false) ---
// Maximum summed intensity across all channels (on 0-255 per-channel scale,
// so 0-1020 for 4 channels) to classify as empty. Only used when the
// adaptive mode is disabled. Retained for reproducibility of analyses run
// with the pre-adaptive macro; the adaptive mode is recommended for new work.
var CFG_EMPTY_SPACE_INTENSITY_MAX = 5;

//----------------------------------------------------------------------------
// 5. SHARED MORPHOLOGICAL CONSTANTS
//----------------------------------------------------------------------------
// Vessel lumen hole detection. After the SMA and capillary masks are built,
// the macro looks for "dark" regions fully enclosed within each vessel ROI
// (the lumen of an arteriole, the cleared centre of a capillary cross-
// section) and SUBTRACTS them from the ROI's area. Lumens are not part of
// the vessel wall, so including them in the ROI area would distort
// downstream measurements: an amyloid load computed as IntDen/Area should
// not divide the integrated amyloid intensity by an area that includes
// empty lumen pixels.
//
// The two parameters below govern this hole detection. They are shared
// between the SMA-vessel and capillary-vessel pipelines because the
// underlying logic is the same in both contexts:
//   - CFG_MIN_HOLE_SIZE: minimum size (px²) for a detected hole to count
//     as a vessel lumen. Holes smaller than this are likely thresholding
//     artefacts (a few stray dark pixels) and are kept inside the ROI area.
//   - CFG_HOLE_DARK_THRESHOLD: a pixel inside the ROI is considered "dark"
//     (and therefore lumen-like) if its intensity is below
//     mean_ROI_intensity * CFG_HOLE_DARK_THRESHOLD.
//     0.5 means a pixel must be at least 50% dimmer than the mean ROI
//     intensity to qualify. Higher values (e.g., 0.7) are more permissive
//     and detect more pixels as lumen; lower values (e.g., 0.3) are stricter.
var CFG_MIN_HOLE_SIZE = 50;
var CFG_HOLE_DARK_THRESHOLD = 0.5;

//----------------------------------------------------------------------------
// 6. COLLAGEN IV DETECTION
//----------------------------------------------------------------------------
// Collagen IV is the marker used to define vessel basement membranes, from
// which both the capillary compartment and the arteriolar/large-vessel
// compartment are ultimately derived. This section contains the largest
// number of parameters in the macro because ColIV staining patterns vary
// widely in intensity across vessel types (capillaries are typically much
// dimmer than arterioles) and because nuclear bleed-through (or perinuclear
// lipofuscin autofluorescence) produces ColIV-like artefacts that need to
// be filtered out.
//
// Pipeline stages within this section:
//   (a) Preprocessing (background subtraction, median filter)
//   (b) Initial mask segmentation via Phansalkar auto local threshold
//   (c) ROI extraction and size/circularity binning
//   (d) Bleed-through correction (optional, see CFG_BLEEDTHROUGH_CORRECTION)
//   (e) Per-bin validation against intensity and/or overlap criteria
//   (f) Morphological closing of surviving ROIs to produce the final mask

// --- Morphological closing radii (Diamond structuring element) ---
// Both radii are used exclusively in the construction of the cleaned
// collagen IV mask (clean_collagen4_mask). They are NOT shared with the
// SMA pipeline, despite the section-5-style naming.
//   - CFG_CLOSING_RADIUS_SMALL (2 px): applied after noise-floor intensity
//     thresholding, when bleed-through correction is enabled and bled
//     nuclei were found. The smaller radius preserves the morphology of
//     thin vessel rims that survive the intensity filter.
//   - CFG_CLOSING_RADIUS_LARGE (3 px): applied in two fallback paths —
//     (i) bleed-through correction enabled but no bled nuclei found, and
//     (ii) bleed-through correction disabled entirely. The larger radius
//     compensates for the absence of intensity filtering by absorbing more
//     small-scale noise into the mask.
var CFG_CLOSING_RADIUS_SMALL = 2;
var CFG_CLOSING_RADIUS_LARGE = 3;

// --- (a) Preprocessing ---
// Rolling-ball background subtraction. The radius is set deliberately small
// because nuclei tend to cluster, and a small rolling-ball radius helps
// separate clustered bleed-through nuclei so that they do not agglutinate
// into a single large blob. Keeping clustered nuclei separated makes them
// easier to distinguish from real capillary cross-sections later in the
// pipeline. (5 px ≈ 2.8 µm at validation calibration.)
var CFG_BG_ROLLING_COLLAGEN = 5;

// Median filter radius (px). Suppresses shot noise before thresholding.
// Small radii preserve thin vessel rims; large radii over-smooth and
// produce gaps in the thresholded mask.
var CFG_MEDIAN_RADIUS_COLLAGEN = 2;

// --- (b) Initial mask segmentation ---
// Phansalkar's auto local threshold. Local thresholding adapts to intensity
// variations across the image, which is essential given the wide dynamic
// range of ColIV signal (bright arterioles vs. dim capillaries in the same
// field). Parameters:
//   - radius: size of the local neighbourhood (px) used to compute the
//     local threshold. The radius for both ColIV and SMA was chosen
//     iteratively to maximize sensitivity without introducing too much
//     noise. The size of the underlying biological structures was a
//     secondary consideration: confocal sections cut vessels at arbitrary
//     angles, so the in-plane footprint of any given structure varies
//     widely from image to image, and a single radius cannot be matched
//     exactly to "the" feature size.
//   - k1: sensitivity parameter. Lower = more sensitive (more pixels pass
//     threshold, including noise). Higher = more conservative.
//     0.5 is Phansalkar's published default and was left unchanged.
//   - k2: dynamic-range parameter. 0.5 is Phansalkar's published default.
// k1, k2 and radius are interdependent — changing one without considering
// the others will likely degrade performance.
var CFG_PHANSALKAR_RADIUS_COLLAGEN = 50;
var CFG_PHANSALKAR_K1_COLLAGEN = 0.5;
var CFG_PHANSALKAR_K2_COLLAGEN = 0.5;

// Minimum ROI size (px²) after initial thresholding. 10 px² corresponds to
// the cross-sectional area of a 2-µm-diameter capillary lumen at validation
// calibration: π × (1 µm)² = 3.14 µm² × (1.76 px/µm)² ≈ 9.7 px² ≈ 10 px².
// This is below the smallest plausible capillary cross-section (the lower
// end of the published lumen-diameter range for hippocampal capillaries,
// 2-8 µm). Pixels or pixel pairs that pass the Phansalkar threshold but
// fall below this size are rejected as thresholding artefacts.
var CFG_COLLAGEN_MIN_SIZE = 10;

// --- (c) ROI size binning ---
// ROIs that survive minimum-size filtering are binned into 5 size classes,
// each of which gets its own downstream validation strategy. Bin boundaries
// (px²) are inclusive at the upper edge of each bin. Values were chosen
// iteratively on the validation dataset to separate capillary-scale
// structures (XS, S) from arteriole/vessel-scale structures (M, L, XL).
//   - XS: <= 150 px² — capillary cross-sections, small wall fragments
//   - S:  <= 300 px² — larger capillary cross-sections, wall segments
//   - M:  <= 500 px² — small arteriole cross-sections
//   - L:  <= 1000 px² — larger arteriole cross-sections
//   - XL: > 1000 px² — large arterioles and arteries
var CFG_COLLAGEN_XS_MAX = 150;
var CFG_COLLAGEN_S_MAX = 300;
var CFG_COLLAGEN_M_MAX = 500;
var CFG_COLLAGEN_L_MAX = 1000;

// --- Circularity thresholds (per size: linear < ELONG, elongated < ROUND, ambiguous, ROUND <= round) ---
// Within each size bin, ROIs are further classified by circularity into
// four shape categories: linear, elongated, ambiguous, round. Round ROIs
// are more likely to be artefacts (bled-through nuclei, punctate noise)
// and must pass a stricter intensity test to survive. Elongated and linear
// ROIs are more likely to represent true vessel cross-sections cut at an
// angle, and are validated with more permissive factors, so that dim
// capillaries are not lost.
//
// Values were obtained by trial and error on the validation dataset. They
// should be treated as starting points for new datasets; the only way to
// know whether they transfer is to compare automated output against manual
// annotation on a small set of images from the target acquisition.
// The round thresholds decrease slightly at larger sizes (0.60 at XS/S →
// 0.55 at L) because larger round structures are less likely to be nucleus
// artefacts than small ones; this small relaxation was chosen empirically.
var CFG_CIRC_XS_LINEAR = 0.25;
var CFG_CIRC_XS_ELONG = 0.40;
var CFG_CIRC_XS_ROUND = 0.60;
var CFG_CIRC_S_LINEAR = 0.25;
var CFG_CIRC_S_ELONG = 0.40;
var CFG_CIRC_S_ROUND = 0.60;
var CFG_CIRC_M_LINEAR = 0.20;
var CFG_CIRC_M_ELONG = 0.40;
var CFG_CIRC_M_ROUND = 0.60;
var CFG_CIRC_L_LINEAR = 0.20;
var CFG_CIRC_L_ELONG = 0.35;
var CFG_CIRC_L_ROUND = 0.55;

// --- (d) Bleed-through correction parameters ---
// Used only when CFG_BLEEDTHROUGH_CORRECTION = true (Section 3). The macro
// identifies ColIV ROIs that co-locate with DAPI nuclei and treats them as
// nuclear-coincident artefacts (DAPI emission tail into the C3 channel,
// and/or perinuclear lipofuscin autofluorescence — see Section 7 intro for
// the full rationale).
//
// CFG_NUCLEI_OVERLAP_THRESH and CFG_COLLAGEN_OVERLAP_THRESH are NOT the
// same operation. They apply to different mask-pair comparisons at
// different pipeline stages:
//
//   - CFG_NUCLEI_OVERLAP_THRESH (used at line ~932):
//       Operates on candidate NUCLEI ROIs versus the DAPI mask. A ROI in
//       the candidate-nucleus set is REJECTED if its overlap with the DAPI
//       mask is below this threshold. Filters the nucleus set itself —
//       the question is "is this candidate really a nucleus?".
//
//   - CFG_COLLAGEN_OVERLAP_THRESH (used at line ~1210, XS-round only):
//       Operates on candidate XS-ROUND COLLAGEN ROIs versus the
//       intensity-cleaned collagen mask. An XS-round collagen ROI is KEPT
//       if its overlap with the cleaned mask is at or above this
//       threshold. Filters the final collagen ROI set — the question is
//       "is this small round ColIV ROI really collagen, or a residual
//       artefact?".
//
// They both happen to default to 70% because both are conservative-overlap
// tests on different mask pairs; the coincidence is empirical rather than
// structural. Lowering either makes the corresponding stage more permissive
// (more ROIs survive); raising it makes the stage stricter.
var CFG_NUCLEI_OVERLAP_THRESH = 70;
var CFG_COLLAGEN_OVERLAP_THRESH = 70;

// Minimum number of bled-through nuclei detected before the overlap-based
// XS-round validation is used. Below this, the macro falls back to Yen
// intensity thresholding for XS-round because the overlap-mask statistics
// become unreliable with too few bleed-through samples. 30 was a pragmatic
// threshold; in the validation dataset, 2 of 12 images fell below it
// (C11.11, C13.6) and used the Yen fallback successfully.
var CFG_COLLAGEN_MIN_BLED_NUCLEI = 30;

// Percentile of the collagen intensity distribution within bled-through
// ROIs used as the noise-floor threshold. 0.90 = 90th percentile. Pixels
// below this intensity are treated as noise in the "cleaned" collagen mask.
// Higher values (e.g., 0.95) are more conservative (keep more real signal
// but also more noise); lower values (e.g., 0.80) are more aggressive.
var CFG_NOISE_PERCENTILE = 0.90;

// --- (e) Per-bin Yen intensity factors ---
// For ROIs not validated by the overlap method (all non-XS-ROUND categories,
// and XS-ROUND in fallback), the Yen auto-threshold of the collagen channel
// is multiplied by a per-bin factor to produce the effective intensity
// requirement. The test compares the ROI's MAXIMUM pixel intensity against
// (mean of Yen mask × factor): a ROI passes if its peak signal is bright
// enough relative to the Yen-derived global threshold.
//
// Because the test uses ROI MAX, larger ROIs have an inherent advantage
// regardless of their typical brightness — more pixels means more chances
// to contain a single pixel that crosses the threshold. The per-bin factors
// compensate for this and encode the macro's confidence in each bin:
//
//   - Round ROIs at XS/S sizes are most often bled-through nuclei (DAPI
//     tail or perinuclear lipofuscin), so they are held to a strict
//     threshold (factor = 1.00) — the ROI's MAX must reach the Yen mean.
//   - Round ROIs at M/L sizes are even more likely to be artefacts: real
//     arterioles are very rarely solid round structures at this scale, and
//     because of the MAX-based test, larger objects have a built-in
//     advantage that needs counteracting. They are therefore held to a
//     stricter factor (1.50) — the ROI's MAX must reach 1.5× the Yen mean.
//   - Linear and elongated ROIs are much more likely to represent real
//     vessel cross-sections cut obliquely; they are held to lenient
//     factors (0.05-0.15) so dim oblique capillaries are not lost.
//   - Ambiguous shapes get an intermediate factor (0.50-0.60).
//   - XL ROIs (>1000 px²) get a moderate factor (0.50) regardless of
//     shape: at this size, the structure is almost always a real large
//     vessel.
//
// All values were obtained by trial and error on the validation dataset
// and should be treated as starting points for new data. The within-bin
// ratios (round >> ambiguous >> elongated > linear) should usually be
// preserved even if absolute magnitudes are recalibrated.
var CFG_YEN_XS_ROUND = 1.00;
var CFG_YEN_XS_AMBIG = 0.60;
var CFG_YEN_XS_ELONG = 0.10;
var CFG_YEN_XS_LINEAR = 0.05;
var CFG_YEN_S_ROUND = 1.00;
var CFG_YEN_S_AMBIG = 0.50;
var CFG_YEN_S_ELONG = 0.10;
var CFG_YEN_S_LINEAR = 0.05;
var CFG_YEN_M_ROUND = 1.50;
var CFG_YEN_M_AMBIG = 0.60;
var CFG_YEN_M_ELONG = 0.15;
var CFG_YEN_M_LINEAR = 0.05;
var CFG_YEN_L_ROUND = 1.50;
var CFG_YEN_L_AMBIG = 0.50;
var CFG_YEN_L_ELONG = 0.10;
var CFG_YEN_L_LINEAR = 0.05;
var CFG_YEN_XL = 0.50;

// --- (f) Morphological closing of surviving ROIs ---
// After per-bin validation, ROIs are grouped into two "runs" and closed
// before being merged into the final collagen4 mask.
//   - POST: closing applied AFTER per-ROI filtering, for the XS+S group.
//     The ordering matters: bleed-through nuclei tend to cluster, and when
//     thresholded they often connect into non-round composite shapes. If
//     closing were applied before filtering, the round-shape exclusion
//     used to identify bleed-through artefacts would be defeated by these
//     pre-fused composites. Filtering first (using the still-discrete
//     round shapes), then closing (to bridge real vessel-rim gaps among
//     the survivors), prevents bleed-through nuclei from contaminating
//     the final ColIV mask.
//   - PRE: closing applied BEFORE per-ROI filtering, for the M+L+XL group.
//     At these sizes, bleed-through nuclei are not the dominant artefact,
//     so closing first allows large arteriole walls with thresholding gaps
//     to be reconstructed before the per-ROI measurement.
var CFG_COLLAGEN_CLOSING_POST = 3;
var CFG_COLLAGEN_CLOSING_PRE = 3;

//----------------------------------------------------------------------------
// 7. DAPI NUCLEI DETECTION
//----------------------------------------------------------------------------
// This section is ONLY executed when CFG_CHANNELS_N = 4 AND
// CFG_BLEEDTHROUGH_CORRECTION = true.
//
// Rationale for DAPI inclusion in the staining protocol: DAPI is a standard
// nuclear counterstain in neuroscience confocal imaging. It is included for
// general tissue orientation, for tissue-quality QC (detecting section
// tears, necrosis, or gross artefacts that would invalidate an image), and
// to provide context for any manual review of the automated output. During
// macro development it became apparent that nuclear-coincident structures
// also appear in the ColIV channel — these may be DAPI emission tail
// crossing into C3 and/or perinuclear lipofuscin autofluorescence (the two
// are not separable from the data we have, and either alone would justify
// the correction). Whichever the source, the artefacts co-localise with
// nuclei, so the DAPI mask provides a natural scaffold for identifying and
// removing them. Nuclei are not a biological compartment of interest in
// this analysis; they are detected solely as a co-localisation reference
// for artefact removal.

// --- Preprocessing ---
// Rolling-ball background subtraction.
var CFG_BG_ROLLING_NUCLEI = 15;

// Median filter radius. Kept small (1 px) because DAPI signal is usually
// bright and relatively clean; heavier filtering would merge adjacent nuclei.
var CFG_MEDIAN_RADIUS_NUCLEI = 1;

// --- Mask segmentation ---
// Phansalkar auto local threshold parameters. As with ColIV and SMA, the
// radius was chosen iteratively to maximise sensitivity without introducing
// excessive noise; matching it exactly to nucleus dimensions is impractical
// because confocal sections cut nuclei at arbitrary angles, so the in-plane
// nucleus footprint varies widely from image to image.
//   - radius: 5 px.
//   - k1, k2: 0.5 each, Phansalkar published defaults, not modified.
// As elsewhere, k1, k2 and radius are interdependent — isolated changes
// are unlikely to improve behaviour.
var CFG_PHANSALKAR_RADIUS_NUCLEI = 5;
var CFG_PHANSALKAR_K1_NUCLEI = 0.5;
var CFG_PHANSALKAR_K2_NUCLEI = 0.5;

// Minimum ROI size (px²) for a thresholded region to count as a nucleus.
// 39 px² ≈ 12.6 µm² at validation calibration, below the cross-section of
// a typical small nucleus. Too-small regions are mostly bright punctate
// noise rather than real nuclear profiles.
var CFG_NUCLEI_MIN_SIZE = 39;

// --- Bled-nuclei size window ---
// After nuclei are detected, the macro finds their intersection with the
// initial collagen IV ROI set. Those collagen ROIs overlapping nuclei are
// candidate "bled-through" artefacts. The macro then filters these
// candidates by size, retaining only candidates in the expected range for
// bleed-through artefacts:
//   - Too small (< MIN): noise, not representative of nucleus-sized
//     bleed-through.
//   - Too large (> MAX): likely a real vessel structure that happens to
//     sit near a nucleus; retaining it as "bleed-through" would remove
//     a real vessel.
// 5-100 px² window was tuned empirically against the validation dataset.
var CFG_BLED_NUCLEI_MIN_SIZE = 5;
var CFG_BLED_NUCLEI_MAX_SIZE = 100;

//----------------------------------------------------------------------------
// 8. SMA DETECTION
//----------------------------------------------------------------------------
// Alpha-smooth muscle actin (SMA) identifies vessels with a contractile
// wall: arterioles, arteries, and (at lower expression) some venules. In
// the IPAD context, SMA-positive vessels are the primary drainage route,
// so accurate SMA detection is essential for the downstream
// arteriole-vs-capillary classification.
//
// SMA signal is more variable across acquisitions than ColIV — antibody
// choice and staining protocol produce wider intensity differences than
// for the ColIV marker. The section therefore includes an adaptive
// thresholding path (CFG_ADAPTIVE_SMA_ENABLED in Section 3) that scales
// with local image intensity.
//
// Pipeline stages within this section:
//   (a) Preprocessing (background subtraction, median filter)
//   (b) Initial mask segmentation via Phansalkar auto local threshold
//   (c) Classification of SMA ROIs by overlap with collagen IV
//   (d) Per-level intensity validation (adaptive or fixed)

// --- (a) Preprocessing ---
// Rolling-ball background subtraction.
var CFG_BG_ROLLING_SMA = 25;

// Median filter radius. 2 px same as ColIV; balances noise suppression
// against preservation of thin circumferential wall staining.
var CFG_MEDIAN_RADIUS_SMA = 2;

// --- (b) Initial mask segmentation ---
// Phansalkar parameters for SMA. As with ColIV and nuclei, the radius was
// chosen iteratively to maximise sensitivity without introducing excessive
// noise; the size of the underlying biological structures was a secondary
// consideration, because confocal sections cut SMA-positive walls at
// arbitrary angles and the in-plane footprint varies widely between images.
//   - radius: 15 px.
//   - k1, k2: 0.5 each, Phansalkar published defaults.
// k1, k2 and radius are interdependent — isolated changes are not
// recommended.
var CFG_PHANSALKAR_RADIUS_SMA = 15;
var CFG_PHANSALKAR_K1_SMA = 0.5;
var CFG_PHANSALKAR_K2_SMA = 0.5;

// Minimum size (px²) for an SMA ROI. 10 px² is the same minimum as ColIV,
// below which objects are almost certainly thresholding noise rather than
// smooth-muscle-actin-bearing tissue.
var CFG_SMA_MIN_SIZE = 10;

// Morphological closing radius applied to the SMA mask before ROI extraction.
// Small value (2 px, Diamond element) bridges gaps in wall staining without
// merging adjacent vessels. Important because SMA immunostaining often
// produces slightly discontinuous wall segments even within a single vessel.
var CFG_SMA_CLOSING_RADIUS = 2;

// --- (c) Classification of SMA ROIs by collagen overlap ---
// SMA ROIs are classified into 5 confidence levels based on their %
// overlap with the auto_local_collagen4 mask. The rationale:
//   - An SMA ROI that overlaps heavily with ColIV is almost certainly a
//     real vessel (both basement membrane AND smooth muscle detected).
//   - An SMA ROI with no ColIV overlap might be a real vessel with weak
//     ColIV staining, or might be non-vascular SMA-positive signal (e.g.,
//     pericytes, staining artefact, bleed-through).
// Higher confidence levels get more permissive thresholds (the macro trusts
// them); lower levels get stricter thresholds (the macro demands strong
// evidence). See the adaptive threshold parameters below for the per-level
// ranges.
//
// Level boundaries (minimum % overlap to qualify for that level):
var CFG_SMA_OL_LEVEL5 = 95;  // Highest confidence: strong ColIV + strong SMA
var CFG_SMA_OL_LEVEL4 = 80;
var CFG_SMA_OL_LEVEL3 = 60;
var CFG_SMA_OL_LEVEL2 = 35;
var CFG_SMA_OL_LEVEL1 = 0;   // Lowest confidence: SMA-only, no ColIV support

// --- (d) Per-level adaptive intensity thresholds ---
// Used when CFG_ADAPTIVE_SMA_ENABLED = true (the default).
// For each overlap level, the effective SMA intensity threshold is a
// function of the mean C1 intensity in the current image:
//   threshold = factor(meanC1)
// where factor(x) linearly interpolates between FACTOR_MIN at INTENSITY_LOW
// and FACTOR_MAX at INTENSITY_HIGH, clamped to [FACTOR_MIN, FACTOR_MAX]
// outside that range. Dimmer images get more permissive factors (closer
// to FACTOR_MIN); brighter images get stricter ones (closer to FACTOR_MAX).
//
// Across all five levels, INTENSITY_LOW = 1 and INTENSITY_HIGH = 30 on the
// 0-255 scale. These anchor points were obtained by trial and error on the
// validation dataset; they bracket the range of mean C1 intensities
// observed across the 12 validation images.
//
// FACTOR ranges differ by level — higher-confidence levels (L5) get more
// permissive ranges (lower minimums), lower-confidence levels (L1) get
// stricter ranges. All factor values were obtained by trial and error.

// Level 5 (>=95% collagen overlap): highest confidence, most permissive threshold
var CFG_ADAPTIVE_L5_INTENSITY_LOW = 1;
var CFG_ADAPTIVE_L5_INTENSITY_HIGH = 30;
var CFG_ADAPTIVE_L5_FACTOR_MIN = 0.25;
var CFG_ADAPTIVE_L5_FACTOR_MAX = 0.65;

// Level 4 (>=80% collagen overlap)
var CFG_ADAPTIVE_L4_INTENSITY_LOW = 1;
var CFG_ADAPTIVE_L4_INTENSITY_HIGH = 30;
var CFG_ADAPTIVE_L4_FACTOR_MIN = 0.30;
var CFG_ADAPTIVE_L4_FACTOR_MAX = 0.65;

// Level 3 (>=60% collagen overlap)
var CFG_ADAPTIVE_L3_INTENSITY_LOW = 1;
var CFG_ADAPTIVE_L3_INTENSITY_HIGH = 30;
var CFG_ADAPTIVE_L3_FACTOR_MIN = 0.35;
var CFG_ADAPTIVE_L3_FACTOR_MAX = 0.70;

// Level 2 (>=35% collagen overlap)
var CFG_ADAPTIVE_L2_INTENSITY_LOW = 1;
var CFG_ADAPTIVE_L2_INTENSITY_HIGH = 30;
var CFG_ADAPTIVE_L2_FACTOR_MIN = 0.35;
var CFG_ADAPTIVE_L2_FACTOR_MAX = 0.85;

// Level 1 (>=0% collagen overlap): lowest confidence, strictest threshold
var CFG_ADAPTIVE_L1_INTENSITY_LOW = 1;
var CFG_ADAPTIVE_L1_INTENSITY_HIGH = 30;
var CFG_ADAPTIVE_L1_FACTOR_MIN = 0.50;
var CFG_ADAPTIVE_L1_FACTOR_MAX = 0.90;

// Fixed (non-adaptive) threshold. Used only when CFG_ADAPTIVE_SMA_ENABLED
// = false. Threshold = Mean of Yen Mask * Factor (applied to ROI Max).
// Retained for reproducibility of analyses run before the adaptive path
// was added. 0.50 is a middle-of-the-road value that worked acceptably
// on the validation dataset when adaptive mode was disabled.
var CFG_SMA_THRESHOLD_FACTOR = 0.50;

//----------------------------------------------------------------------------
// 9. LARGE VESSEL / CAPILLARY CLASSIFICATION
//----------------------------------------------------------------------------
// After ColIV and SMA masks are built, each ColIV-positive ROI is classified
// as either:
//   - a "large vessel" if it overlaps with the SMA mask, with its
//     footprint defined as the ColIV rim minus the SMA signal itself; or
//   - a "capillary" (SMA-negative basement membrane) otherwise.
// "Large vessel" includes arterioles and arteries, and may also include
// SMA-positive venules (some venules express low levels of SMA). Venules
// with little or no SMA signal will fall below CFG_SMA_EXCLUSION_THRESH
// and be classified as capillaries by default. The macro does not attempt
// a separate venule compartment because at the resolution of the
// validation acquisition, SMA-low venules and capillaries are not reliably
// separable, and any rule for splitting them would introduce uncontrolled
// observer subjectivity. This section contains the two parameters that
// govern the large-vessel-vs-capillary classification.

// Minimum size (px²) for a structure to count as a large vessel after the
// SMA footprint is subtracted. Below this, the surviving ring-fragment is
// almost always fragmented wall material rather than a biologically
// meaningful arteriole outline. 150 px² ≈ 48 µm² at validation calibration.
// Lower it only if a higher-magnification acquisition is producing large
// vessel fragments that fall below threshold.
var CFG_LV_MIN_SIZE = 150;

// Minimum % overlap between a ColIV ROI and the SMA mask for the ColIV
// ROI to be excluded from the capillary compartment (i.e., classified as a
// large vessel instead). The default of 0 means "any overlap at all" — any
// ColIV ROI with even a single pixel of SMA overlap is removed from the
// capillary set. This is the most aggressive exclusion threshold and
// minimises false-positive capillaries (capillaries incorrectly including
// SMA-positive vessel fragments).
var CFG_SMA_EXCLUSION_THRESH = 0;

//----------------------------------------------------------------------------
// 10. AMYLOID DETECTION
//----------------------------------------------------------------------------
// Amyloid is the primary quantitative output of the pipeline. Unlike the
// other channels, amyloid detection uses a simple global intensity threshold
// rather than a local auto-threshold, because the amyloid signal distribution
// is strongly bimodal (background vs. positive deposits) and a global cutoff
// is both easier to calibrate and more reproducible across acquisitions.

// --- Preprocessing ---
// Rolling-ball background subtraction (same 25 px radius as SMA).
var CFG_BG_ROLLING_AMYLOID = 25;

// Median filter radius. Kept small (1 px) to preserve the morphology of
// amyloid deposits, which can be small, punctate, and clinically meaningful
// at scales of a few pixels.
var CFG_MEDIAN_RADIUS_AMYLOID = 1;

// --- Detection ---
// CFG_AMYLOID_THRESHOLD is the single most dataset-dependent parameter in
// the entire macro. It sets the lower cutoff for amyloid classification on
// the 0-255 intensity scale; pixels with amyloid-channel intensity >= this
// value are classified as amyloid-positive.
//
// Recommended recalibration procedure for a new dataset: image your control
// group (no fluorescent amyloid injected) with the same laser power and
// detector gain you will use for experimental animals, measure the C2
// autofluorescence distribution in tissue regions, and set this threshold
// just above the autofluorescence ceiling (e.g., max + small buffer, or
// 95th percentile of C2 intensity in control tissue).
var CFG_AMYLOID_THRESHOLD = 26;

//============================================================================
//=================== RUNTIME STATE ==========================================
//============================================================================
// These variables are NOT configuration — they are populated during processing
// and read by downstream code (QC logging, adaptive threshold computation).
// Declared at global scope because ImageJ macro language has no module system.

// Mean channel intensities (measured on MAX projections, excluding empty space)
var meanIntensity_C1 = 0;
var meanIntensity_C2 = 0;
var meanIntensity_C3 = 0;

// Empty-space summary statistics
var emptySpaceArea = 0;
var emptySpaceMeanIntensity = 0;

// Per-level adaptive SMA factors and counts, populated during SMA processing
var sma_factor_L1 = 0; var sma_factor_L2 = 0; var sma_factor_L3 = 0; var sma_factor_L4 = 0; var sma_factor_L5 = 0;
var sma_thresh_L1 = 0; var sma_thresh_L2 = 0; var sma_thresh_L3 = 0; var sma_thresh_L4 = 0; var sma_thresh_L5 = 0;
var sma_kept_L1 = 0; var sma_kept_L2 = 0; var sma_kept_L3 = 0; var sma_kept_L4 = 0; var sma_kept_L5 = 0;
var sma_removed_L1 = 0; var sma_removed_L2 = 0; var sma_removed_L3 = 0; var sma_removed_L4 = 0; var sma_removed_L5 = 0;

//============================================================================
//=================== HELPER FUNCTIONS =======================================
//============================================================================

var macroStartTime = getTime();

function getElapsedTime() {
    var elapsed = (getTime() - macroStartTime) / 1000;
    var minutes = floor(elapsed / 60);
    var seconds = elapsed - (minutes * 60);
    if (minutes > 0) {
        return "[" + d2s(minutes, 0) + "m " + d2s(seconds, 1) + "s]";
    } else {
        return "[" + d2s(seconds, 1) + "s]";
    }
}

function logStep(message) {
    print(getElapsedTime() + " " + message);
}

function safeClose(imageName) {
    if (isOpen(imageName)) { selectImage(imageName); close(); }
}

function computeAdaptiveFactor(meanC1val, intLow, intHigh, factorMin, factorMax) {
    var range = intHigh - intLow;
    if (range > 0) {
        var normalized = (meanC1val - intLow) / range;
        if (normalized < 0) normalized = 0;
        if (normalized > 1) normalized = 1;
        return factorMin + (factorMax - factorMin) * normalized;
    }
    return factorMin;
}

//============================================================================
// v9.0: SESSION-SCOPE OUTPUT DIRECTORY RESOLUTION
//============================================================================
// Resolves output directory paths with priority: CFG > JVM system property > prompt.
// JVM system properties are alive for the duration of the Fiji process; they are
// cleared automatically when Fiji closes. This gives "once per Fiji launch"
// prompt behavior without persisting user choices to disk between sessions.
//
// Returns the resolved path (with trailing slash). Exits the macro if the user
// cancels the dialog.
function resolveSessionDir(cfgVal, propKey, dialogTitle) {
    // 1. CFG value set in code: use it directly, no prompt
    if (cfgVal != "") {
        return cfgVal;
    }
    // 2. Check JVM system property (survives for this Fiji session only)
    var saved = call("java.lang.System.getProperty", propKey, "");
    if (saved != "" && File.exists(saved)) {
        return saved;
    }
    // 3. Prompt the user
    var chosen = getDirectory(dialogTitle);
    if (chosen == "") {
        exit("No directory selected for '" + dialogTitle + "'. Macro cancelled.");
    }
    // Save to JVM system property for the rest of this Fiji session
    call("java.lang.System.setProperty", propKey, chosen);
    return chosen;
}

function resolveOutputDirs() {
    CFG_QC_OUTPUT_DIR = resolveSessionDir(CFG_QC_OUTPUT_DIR, "ipad.qc_dir", "Select QC output folder");
    if (CFG_VALIDATION_MODE) {
        CFG_VALIDATION_EXPORT_DIR = resolveSessionDir(CFG_VALIDATION_EXPORT_DIR, "ipad.val_dir", "Select validation export folder");
    }
}

//============================================================================
// v9.0: CHANNEL-TO-MARKER ASSIGNMENT
//============================================================================
// The macro's internal logic references channels by canonical names:
//   C1-A1 = SMA, C2-A1 = Amyloid, C3-A1 = ColIV, C4-A1 = DAPI
// If the user's acquisition uses a different channel order, the CFG_CHANNEL_*
// variables specify the mapping. After splitting channels, we rename them to
// the canonical internal names so the rest of the macro runs unchanged.

function validateChannelAssignment() {
    // v9.0: Validate CFG_CHANNELS_N first, then validate only the active channels.
    if (CFG_CHANNELS_N != 3 && CFG_CHANNELS_N != 4) {
        exit("Invalid CFG_CHANNELS_N = " + CFG_CHANNELS_N + ". Must be 3 or 4.");
    }
    
    // Collect the active channel assignments (DAPI excluded when N=3)
    var vals;
    var names;
    if (CFG_CHANNELS_N == 4) {
        vals = newArray(CFG_CHANNEL_SMA, CFG_CHANNEL_AMYLOID, CFG_CHANNEL_COLIV, CFG_CHANNEL_DAPI);
        names = newArray("SMA", "Amyloid", "ColIV", "DAPI");
    } else {
        vals = newArray(CFG_CHANNEL_SMA, CFG_CHANNEL_AMYLOID, CFG_CHANNEL_COLIV);
        names = newArray("SMA", "Amyloid", "ColIV");
    }
    var nActive = vals.length;
    
    // Each active value must be in 1..CFG_CHANNELS_N
    for (var i = 0; i < nActive; i++) {
        if (vals[i] < 1 || vals[i] > CFG_CHANNELS_N) {
            exit("Invalid channel assignment: CFG_CHANNEL_" + names[i] + " = " + vals[i] +
                 ". Must be between 1 and " + CFG_CHANNELS_N + ".");
        }
    }
    // All active values must be unique
    for (var i = 0; i < nActive; i++) {
        for (var j = i + 1; j < nActive; j++) {
            if (vals[i] == vals[j]) {
                exit("Duplicate channel assignment: CFG_CHANNEL_" + names[i] +
                     " and CFG_CHANNEL_" + names[j] + " both set to " + vals[i] +
                     ". Each channel must be mapped to a unique marker.");
            }
        }
    }
}

// v9.0: Resolve final channel configuration and reconcile CFG_BLEEDTHROUGH_CORRECTION
// with CFG_CHANNELS_N. Must be called after validateChannelAssignment() and before
// any DAPI-dependent code runs. If the user set CFG_CHANNELS_N=3 but left
// CFG_BLEEDTHROUGH_CORRECTION=true, the latter is force-disabled with a warning,
// since there is no DAPI channel available to calibrate it.
function resolveChannelConfig() {
    if (CFG_CHANNELS_N == 3 && CFG_BLEEDTHROUGH_CORRECTION) {
        logStep("WARNING: CFG_CHANNELS_N = 3 but CFG_BLEEDTHROUGH_CORRECTION = true.");
        logStep("  - No DAPI channel is available to calibrate the correction.");
        logStep("  - Auto-disabling bleed-through correction for this run.");
        CFG_BLEEDTHROUGH_CORRECTION = false;
    }
}

// Renames split channels from their acquisition-order names (C<n>-A1) to the
// canonical internal names (C1-A1 = SMA, etc.). Uses a two-step rename
// (to temporary names first, then to final names) to avoid name collisions
// when the input order already partially matches the canonical order.
// v9.0: Handles 3 or 4 channels per CFG_CHANNELS_N. The actual channel count
// of the input image is verified against CFG_CHANNELS_N at macro start before
// this function is called, so we can rely on exactly CFG_CHANNELS_N channels
// being present here.
function remapChannels() {
    // Step 1: rename each active input channel to a unique temporary name
    selectImage("C" + CFG_CHANNEL_SMA + "-A1");     rename("_TMP_SMA");
    selectImage("C" + CFG_CHANNEL_AMYLOID + "-A1"); rename("_TMP_AMYLOID");
    selectImage("C" + CFG_CHANNEL_COLIV + "-A1");   rename("_TMP_COLIV");
    if (CFG_CHANNELS_N == 4) {
        selectImage("C" + CFG_CHANNEL_DAPI + "-A1"); rename("_TMP_DAPI");
    }
    
    // Step 2: rename to canonical internal names
    selectImage("_TMP_SMA");     rename("C1-A1");
    selectImage("_TMP_AMYLOID"); rename("C2-A1");
    selectImage("_TMP_COLIV");   rename("C3-A1");
    if (CFG_CHANNELS_N == 4) {
        selectImage("_TMP_DAPI"); rename("C4-A1");
    }
}

function logConfig() {
    print("================================================================================");
    print("CONFIGURATION PARAMETERS");
    print("================================================================================");
    print("Size Thresholds (px2):");
    print("  - Collagen min: " + CFG_COLLAGEN_MIN_SIZE);
    print("  - Nuclei min: " + CFG_NUCLEI_MIN_SIZE);
    print("  - Bled nuclei: " + CFG_BLED_NUCLEI_MIN_SIZE + "-" + CFG_BLED_NUCLEI_MAX_SIZE);
    print("  - SMA min: " + CFG_SMA_MIN_SIZE);
    print("  - Large vessel min: " + CFG_LV_MIN_SIZE);
    print("Overlap Thresholds (%):");
    print("  - Nuclei bleed-through: " + CFG_NUCLEI_OVERLAP_THRESH + "%");
    print("Collagen Filtering (5 sizes x 4 shapes):");
    print("  - Size: XS<=" + CFG_COLLAGEN_XS_MAX + " | S<=" + CFG_COLLAGEN_S_MAX + " | M<=" + CFG_COLLAGEN_M_MAX + " | L<=" + CFG_COLLAGEN_L_MAX + " | XL>" + CFG_COLLAGEN_L_MAX);
    print("  - Circ XS: lin<" + CFG_CIRC_XS_LINEAR + " elong<" + CFG_CIRC_XS_ELONG + " ambig round>=" + CFG_CIRC_XS_ROUND);
    print("  - Circ S:  lin<" + CFG_CIRC_S_LINEAR + " elong<" + CFG_CIRC_S_ELONG + " ambig round>=" + CFG_CIRC_S_ROUND);
    print("  - Circ M:  lin<" + CFG_CIRC_M_LINEAR + " elong<" + CFG_CIRC_M_ELONG + " ambig round>=" + CFG_CIRC_M_ROUND);
    print("  - Circ L:  lin<" + CFG_CIRC_L_LINEAR + " elong<" + CFG_CIRC_L_ELONG + " ambig round>=" + CFG_CIRC_L_ROUND);
    print("  - XS round: overlap>=" + CFG_COLLAGEN_OVERLAP_THRESH + "% (Yen*" + CFG_YEN_XS_ROUND + " backup if bled nuclei<" + CFG_COLLAGEN_MIN_BLED_NUCLEI + ")");
    print("  - Yen XS: R=" + CFG_YEN_XS_ROUND + " A=" + CFG_YEN_XS_AMBIG + " E=" + CFG_YEN_XS_ELONG + " L=" + CFG_YEN_XS_LINEAR);
    print("  - Yen S:  R=" + CFG_YEN_S_ROUND + " A=" + CFG_YEN_S_AMBIG + " E=" + CFG_YEN_S_ELONG + " L=" + CFG_YEN_S_LINEAR);
    print("  - Yen M:  R=" + CFG_YEN_M_ROUND + " A=" + CFG_YEN_M_AMBIG + " E=" + CFG_YEN_M_ELONG + " L=" + CFG_YEN_M_LINEAR);
    print("  - Yen L:  R=" + CFG_YEN_L_ROUND + " A=" + CFG_YEN_L_AMBIG + " E=" + CFG_YEN_L_ELONG + " L=" + CFG_YEN_L_LINEAR);
    print("  - Yen XL: " + CFG_YEN_XL);
    print("  - Closing: post(XS+S)=" + CFG_COLLAGEN_CLOSING_POST + " | pre(M+L+XL)=" + CFG_COLLAGEN_CLOSING_PRE);
    print("  - SMA closing: Diamond r=" + CFG_SMA_CLOSING_RADIUS);
    print("  - SMA-collagen overlap levels: L5>=" + CFG_SMA_OL_LEVEL5 + "% L4>=" + CFG_SMA_OL_LEVEL4 + "% L3>=" + CFG_SMA_OL_LEVEL3 + "% L2>=" + CFG_SMA_OL_LEVEL2 + "% L1>=" + CFG_SMA_OL_LEVEL1 + "%");
    print("  - SMA exclusion: " + CFG_SMA_EXCLUSION_THRESH + "%");
    print("Intensity Parameters:");
    print("  - Noise percentile: " + (CFG_NOISE_PERCENTILE * 100) + "%");
    print("  - SMA Threshold Factor (fixed, when adaptive=false): " + CFG_SMA_THRESHOLD_FACTOR);
    print("  - Adaptive SMA Enabled: " + CFG_ADAPTIVE_SMA_ENABLED);
    if (CFG_ADAPTIVE_SMA_ENABLED) {
        print("    - L5: intLow=" + CFG_ADAPTIVE_L5_INTENSITY_LOW + " intHigh=" + CFG_ADAPTIVE_L5_INTENSITY_HIGH + " fMin=" + CFG_ADAPTIVE_L5_FACTOR_MIN + " fMax=" + CFG_ADAPTIVE_L5_FACTOR_MAX);
        print("    - L4: intLow=" + CFG_ADAPTIVE_L4_INTENSITY_LOW + " intHigh=" + CFG_ADAPTIVE_L4_INTENSITY_HIGH + " fMin=" + CFG_ADAPTIVE_L4_FACTOR_MIN + " fMax=" + CFG_ADAPTIVE_L4_FACTOR_MAX);
        print("    - L3: intLow=" + CFG_ADAPTIVE_L3_INTENSITY_LOW + " intHigh=" + CFG_ADAPTIVE_L3_INTENSITY_HIGH + " fMin=" + CFG_ADAPTIVE_L3_FACTOR_MIN + " fMax=" + CFG_ADAPTIVE_L3_FACTOR_MAX);
        print("    - L2: intLow=" + CFG_ADAPTIVE_L2_INTENSITY_LOW + " intHigh=" + CFG_ADAPTIVE_L2_INTENSITY_HIGH + " fMin=" + CFG_ADAPTIVE_L2_FACTOR_MIN + " fMax=" + CFG_ADAPTIVE_L2_FACTOR_MAX);
        print("    - L1: intLow=" + CFG_ADAPTIVE_L1_INTENSITY_LOW + " intHigh=" + CFG_ADAPTIVE_L1_INTENSITY_HIGH + " fMin=" + CFG_ADAPTIVE_L1_FACTOR_MIN + " fMax=" + CFG_ADAPTIVE_L1_FACTOR_MAX);
    }
    print("  - Amyloid Threshold: >= " + CFG_AMYLOID_THRESHOLD);
    print("  - Empty Space Min Size: " + CFG_EMPTY_SPACE_MIN_SIZE + " px2");
    print("  - Empty Space Median Radius: " + CFG_EMPTY_SPACE_MEDIAN_RADIUS);
    print("  - Adaptive Empty Spaces Enabled: " + CFG_ADAPTIVE_EMPTY_SPACES_ENABLED);
    if (CFG_ADAPTIVE_EMPTY_SPACES_ENABLED) {
        if (CFG_CHANNELS_N == 4) {
            print("  - Empty Space Tolerance: C1=" + CFG_EMPTY_SPACES_TOLERANCE_C1 + " C2=" + CFG_EMPTY_SPACES_TOLERANCE_C2 + " C3=" + CFG_EMPTY_SPACES_TOLERANCE_C3 + " C4=" + CFG_EMPTY_SPACES_TOLERANCE_C4);
        } else {
            print("  - Empty Space Tolerance: C1=" + CFG_EMPTY_SPACES_TOLERANCE_C1 + " C2=" + CFG_EMPTY_SPACES_TOLERANCE_C2 + " C3=" + CFG_EMPTY_SPACES_TOLERANCE_C3);
        }
    } else {
        print("  - Empty Space Intensity Max (static): " + CFG_EMPTY_SPACE_INTENSITY_MAX);
    }
    print("Morphological Parameters:");
    print("  - Phansalkar (collagen): r=" + CFG_PHANSALKAR_RADIUS_COLLAGEN + " k1=" + CFG_PHANSALKAR_K1_COLLAGEN + " k2=" + CFG_PHANSALKAR_K2_COLLAGEN);
    print("  - Phansalkar (nuclei): r=" + CFG_PHANSALKAR_RADIUS_NUCLEI + " k1=" + CFG_PHANSALKAR_K1_NUCLEI + " k2=" + CFG_PHANSALKAR_K2_NUCLEI);
    print("  - Phansalkar (SMA): r=" + CFG_PHANSALKAR_RADIUS_SMA + " k1=" + CFG_PHANSALKAR_K1_SMA + " k2=" + CFG_PHANSALKAR_K2_SMA);
    print("  - Closing radius (small): " + CFG_CLOSING_RADIUS_SMALL);
    print("  - Closing radius (large): " + CFG_CLOSING_RADIUS_LARGE);
    print("  - Minimum hole size: " + CFG_MIN_HOLE_SIZE + " px2");
    print("  - Hole dark threshold: " + (CFG_HOLE_DARK_THRESHOLD * 100) + "% of mean intensity");
    print("Background Subtraction (rolling ball):");
    print("  - Collagen: " + CFG_BG_ROLLING_COLLAGEN);
    print("  - Nuclei: " + CFG_BG_ROLLING_NUCLEI);
    print("  - SMA: " + CFG_BG_ROLLING_SMA);
    print("  - Amyloid: " + CFG_BG_ROLLING_AMYLOID);
    print("Median Filter Radii:");
    print("  - Collagen: " + CFG_MEDIAN_RADIUS_COLLAGEN);
    print("  - SMA: " + CFG_MEDIAN_RADIUS_SMA);
    print("  - Amyloid: " + CFG_MEDIAN_RADIUS_AMYLOID);
    print("  - Nuclei: " + CFG_MEDIAN_RADIUS_NUCLEI);
    print("Validation Mode: " + CFG_VALIDATION_MODE);
    print("Output directories:");
    print("  - QC: " + CFG_QC_OUTPUT_DIR);
    if (CFG_VALIDATION_MODE) {
        print("  - Validation Export: " + CFG_VALIDATION_EXPORT_DIR);
    }
    print("Channel Assignment:");
    print("  - N channels = " + CFG_CHANNELS_N);
    print("  - SMA     = input channel " + CFG_CHANNEL_SMA);
    print("  - Amyloid = input channel " + CFG_CHANNEL_AMYLOID);
    print("  - ColIV   = input channel " + CFG_CHANNEL_COLIV);
    if (CFG_CHANNELS_N == 4) {
        print("  - DAPI    = input channel " + CFG_CHANNEL_DAPI);
    } else {
        print("  - DAPI    = N/A (3-channel mode)");
    }
    print("DAPI Bleed-Through Correction: " + CFG_BLEEDTHROUGH_CORRECTION);
    print("================================================================================");
}

//============================================================================
//=================== MAIN MACRO STARTS HERE =================================
//============================================================================
//
// PIPELINE STRUCTURE
// ---------------------------------------------------------------------------
// The pipeline below implements the 20 sequential processing steps grouped
// into 5 phases shown in the manuscript's pipeline flowchart figure. Each
// step is anchored by a "// STEP N: ..." marker at the canonical point
// where its main work begins.
//
// Important: the figure presents the pipeline organised by OPERATION
// (e.g., Step 4: median filter ALL channels → Step 5: rolling-ball ALL
// channels → Step 6: ColIV threshold → ...), while the code below is
// organised by CHANNEL (preprocess one channel completely, then the next).
// Both views describe the same algorithm. Two consequences for readers:
//   (a) Steps 4 and 5 (median filter, rolling-ball) appear once per channel
//       in the code — markers indicate which channel each occurrence covers.
//   (b) Step 9 (empty space detection) executes before Steps 6-8 (channel
//       thresholding) in the code, even though it sits later in the figure's
//       Phase 2; this is intentional because empty space is computed from
//       raw signal and is needed before per-channel processing.

// v9.0: Resolve output directories (prompts once per Fiji session if CFG blank)
resolveOutputDirs();

// v9.0: Validate channel assignment (exits with clear error on misconfiguration)
validateChannelAssignment();

// v9.0: Reconcile CFG_CHANNELS_N with CFG_BLEEDTHROUGH_CORRECTION (auto-disable if 3-channel)
resolveChannelConfig();

//============================================================================
// PHASE 1: PRE-PROCESSING AND PROJECTION (Steps 1–5)
//============================================================================

// STEP 1: Load 4-channel confocal z-stack
// The image is loaded by the caller before the macro runs; the macro
// operates on the current foreground image.

run("Set Scale...", "distance=0 known=0 unit=pixel");
roiManager("reset");
run("Clear Results");
setOption("BlackBackground", true);
var originalName = getTitle();

// Global variables to store SMA statistics for QC
var sma_roi_count_for_qc = 0;
var sma_areas_for_qc = newArray(1000);
var sma_areas_count = 0;

print("================================================================================");
print("MACRO START: " + originalName + " (v1.0)");
print("================================================================================");
logConfig();

logStep("Processing image: " + originalName);

var imgWidth = getWidth();
var imgHeight = getHeight();
var imgSlices = nSlices;
logStep("Image dimensions: " + imgWidth + " x " + imgHeight + " x " + imgSlices + " slices");

// v9.0: STRICT CHANNEL COUNT VALIDATION
// Verify that the actual number of channels in the input image matches
// CFG_CHANNELS_N. Exit immediately on mismatch — silently processing an
// image with a different channel count than configured would either crash
// later (if N > actual) or silently drop channel data (if N < actual).
var _dimW, _dimH, _dimC, _dimZ, _dimT;
getDimensions(_dimW, _dimH, _dimC, _dimZ, _dimT);
if (_dimC != CFG_CHANNELS_N) {
    exit("Channel count mismatch: image '" + originalName + "' has " + _dimC +
         " channel(s), but CFG_CHANNELS_N = " + CFG_CHANNELS_N + ".\n" +
         "Either update CFG_CHANNELS_N to match your input data, or use an input " +
         "image with " + CFG_CHANNELS_N + " channel(s).");
}
logStep("  - Channel count verified: " + _dimC + " channels (matches CFG_CHANNELS_N=" + CFG_CHANNELS_N + ")");

//============================================================================
// DUPLICATE AND SPLIT CHANNELS
//============================================================================
// STEP 2: Split channels into independent stacks
logStep("Duplicating and splitting channels...");
run("Duplicate...", "duplicate");
rename("A1");
run("Split Channels");

// v9.0: Remap split channels from acquisition order to canonical internal names
// (C1-A1 = SMA, C2-A1 = Amyloid, C3-A1 = ColIV, C4-A1 = DAPI)
remapChannels();

logStep("Channels split and remapped:");
logStep("  - Input C" + CFG_CHANNEL_SMA     + " -> C1-A1: SMA (Smooth Muscle Actin)");
logStep("  - Input C" + CFG_CHANNEL_AMYLOID + " -> C2-A1: Amyloid");
logStep("  - Input C" + CFG_CHANNEL_COLIV   + " -> C3-A1: Collagen-4 (Basement membrane)");
if (CFG_CHANNELS_N == 4) {
    logStep("  - Input C" + CFG_CHANNEL_DAPI    + " -> C4-A1: DAPI (Nuclei)");
} else {
    logStep("  - DAPI: not used (3-channel mode)");
}

//============================================================================
// RAW MAX PROJECTIONS (before any preprocessing) for validation intensity analysis
//============================================================================
if (CFG_VALIDATION_MODE) {
    logStep("Creating raw MAX projections for validation...");
    selectImage("C1-A1");
    run("Z Project...", "projection=[Max Intensity]");
    rename("RAW_MAX_C1");
    
    selectImage("C2-A1");
    run("Z Project...", "projection=[Max Intensity]");
    rename("RAW_MAX_C2");
    
    selectImage("C3-A1");
    run("Z Project...", "projection=[Max Intensity]");
    rename("RAW_MAX_C3");
    logStep("  - 3 raw MAX projections created (C1, C2, C3)");
}

//============================================================================
// CREATE MEDIAN-FILTERED MAX PROJECTIONS FOR EMPTY SPACE DETECTION
//============================================================================
// STEP 3: MAX intensity projection (one 2D image per channel)
// STEP 4: Median filtering (applied to stacks before projection, with the
//         empty-space-specific radius CFG_EMPTY_SPACE_MEDIAN_RADIUS rather
//         than the per-channel radii used later in per-channel preprocessing)
// Note: Steps 3 and 4 are also applied (with channel-specific median radii)
// inside each per-channel preprocessing section below. This block produces
// a separate set of MAX projections specifically for empty-space detection.
logStep("Preparing empty space detection (median r=" + CFG_EMPTY_SPACE_MEDIAN_RADIUS + " on stacks)...");

selectImage("C1-A1");
run("Duplicate...", "duplicate title=ES_C1_stack");
run("Median...", "radius=" + CFG_EMPTY_SPACE_MEDIAN_RADIUS + " stack");
run("Z Project...", "projection=[Max Intensity]");
rename("ES_C1");
safeClose("ES_C1_stack");

selectImage("C2-A1");
run("Duplicate...", "duplicate title=ES_C2_stack");
run("Median...", "radius=" + CFG_EMPTY_SPACE_MEDIAN_RADIUS + " stack");
run("Z Project...", "projection=[Max Intensity]");
rename("ES_C2");
safeClose("ES_C2_stack");

selectImage("C3-A1");
run("Duplicate...", "duplicate title=ES_C3_stack");
run("Median...", "radius=" + CFG_EMPTY_SPACE_MEDIAN_RADIUS + " stack");
run("Z Project...", "projection=[Max Intensity]");
rename("ES_C3");
safeClose("ES_C3_stack");

if (CFG_CHANNELS_N == 4) {
    selectImage("C4-A1");
    run("Duplicate...", "duplicate title=ES_C4_stack");
    run("Median...", "radius=" + CFG_EMPTY_SPACE_MEDIAN_RADIUS + " stack");
    run("Z Project...", "projection=[Max Intensity]");
    rename("ES_C4");
    safeClose("ES_C4_stack");
}

logStep("  - " + CFG_CHANNELS_N + " median-filtered MAX projections created");

//============================================================================
// PHASE 2: CHANNEL SEGMENTATION (Steps 6–9)
//============================================================================
// In the figure, Phase 2 contains Steps 6 (ColIV thresholding), 7 (DAPI
// thresholding), 8 (SMA thresholding + closing), and 9 (empty space). In
// the code, Step 9 executes first because it operates on raw stacks before
// any per-channel preprocessing, while Steps 6-8 require the per-channel
// preprocessed images produced in the sections that follow.

//============================================================================
// EMPTY SPACE DETECTION
//============================================================================
// STEP 9: Empty space — per-channel minimum + tolerance, AND across channels
if (CFG_ADAPTIVE_EMPTY_SPACES_ENABLED) {
    // ADAPTIVE MODE: per-channel minimum + tolerance, AND overlap
    logStep("Detecting empty spaces (ADAPTIVE: per-channel min + tolerance, AND overlap)...");
    
    selectImage("ES_C1");
    var esA1, esM1, minC1, esX1;
    getStatistics(esA1, esM1, minC1, esX1);
    var esThreshC1 = minC1 + CFG_EMPTY_SPACES_TOLERANCE_C1;
    run("Duplicate...", "title=ES_mask_C1");
    selectImage("ES_mask_C1");
    setThreshold(0, esThreshC1);
    setOption("BlackBackground", true);
    run("Convert to Mask");
    logStep("  - C1: min=" + minC1 + " tol=" + CFG_EMPTY_SPACES_TOLERANCE_C1 + " thresh=0-" + esThreshC1);
    
    selectImage("ES_C2");
    var esA2, esM2, minC2, esX2;
    getStatistics(esA2, esM2, minC2, esX2);
    var esThreshC2 = minC2 + CFG_EMPTY_SPACES_TOLERANCE_C2;
    run("Duplicate...", "title=ES_mask_C2");
    selectImage("ES_mask_C2");
    setThreshold(0, esThreshC2);
    setOption("BlackBackground", true);
    run("Convert to Mask");
    logStep("  - C2: min=" + minC2 + " tol=" + CFG_EMPTY_SPACES_TOLERANCE_C2 + " thresh=0-" + esThreshC2);
    
    selectImage("ES_C3");
    var esA3, esM3, minC3, esX3;
    getStatistics(esA3, esM3, minC3, esX3);
    var esThreshC3 = minC3 + CFG_EMPTY_SPACES_TOLERANCE_C3;
    run("Duplicate...", "title=ES_mask_C3");
    selectImage("ES_mask_C3");
    setThreshold(0, esThreshC3);
    setOption("BlackBackground", true);
    run("Convert to Mask");
    logStep("  - C3: min=" + minC3 + " tol=" + CFG_EMPTY_SPACES_TOLERANCE_C3 + " thresh=0-" + esThreshC3);
    
    if (CFG_CHANNELS_N == 4) {
        selectImage("ES_C4");
        var esA4, esM4, minC4, esX4;
        getStatistics(esA4, esM4, minC4, esX4);
        var esThreshC4 = minC4 + CFG_EMPTY_SPACES_TOLERANCE_C4;
        run("Duplicate...", "title=ES_mask_C4");
        selectImage("ES_mask_C4");
        setThreshold(0, esThreshC4);
        setOption("BlackBackground", true);
        run("Convert to Mask");
        logStep("  - C4: min=" + minC4 + " tol=" + CFG_EMPTY_SPACES_TOLERANCE_C4 + " thresh=0-" + esThreshC4);
    }
    
    // v9.0: AND chain adapts to 3 or 4 channel mode
    imageCalculator("AND create", "ES_mask_C1", "ES_mask_C2");
    rename("es_and_12");
    if (CFG_CHANNELS_N == 4) {
        imageCalculator("AND create", "es_and_12", "ES_mask_C3");
        rename("es_and_123");
        safeClose("es_and_12");
        imageCalculator("AND create", "es_and_123", "ES_mask_C4");
        rename("empty_space_raw");
        safeClose("es_and_123");
    } else {
        imageCalculator("AND create", "es_and_12", "ES_mask_C3");
        rename("empty_space_raw");
        safeClose("es_and_12");
    }
    safeClose("ES_mask_C1");
    safeClose("ES_mask_C2");
    safeClose("ES_mask_C3");
    safeClose("ES_mask_C4");
    logStep("  - AND-ed all " + CFG_CHANNELS_N + " channel masks");
    
    // Size filter
    selectImage("empty_space_raw");
    run("Analyze Particles...", "size=" + CFG_EMPTY_SPACE_MIN_SIZE + "-Infinity show=Masks");
    if (isOpen("Mask of empty_space_raw")) {
        selectImage("Mask of empty_space_raw");
        rename("empty_space_mask");
    } else {
        newImage("empty_space_mask", "8-bit black", imgWidth, imgHeight, 1);
    }
    safeClose("empty_space_raw");
    
    // Measure QC on un-normalized ES channels
    selectImage("empty_space_mask");
    var emptyStatArea, emptyStatMean;
    getStatistics(emptyStatArea, emptyStatMean);
    emptySpaceArea = emptyStatArea * emptyStatMean / 255;
    
    if (emptySpaceArea > 0) {
        selectImage("empty_space_mask");
        run("Create Selection");
        if (selectionType() != -1) {
            var esSumMean = 0;
            var esTmpA, esTmpM;
            selectImage("ES_C1"); run("Restore Selection"); getStatistics(esTmpA, esTmpM); esSumMean += esTmpM;
            selectImage("ES_C2"); run("Restore Selection"); getStatistics(esTmpA, esTmpM); esSumMean += esTmpM;
            selectImage("ES_C3"); run("Restore Selection"); getStatistics(esTmpA, esTmpM); esSumMean += esTmpM;
            if (CFG_CHANNELS_N == 4) {
                selectImage("ES_C4"); run("Restore Selection"); getStatistics(esTmpA, esTmpM); esSumMean += esTmpM;
            }
            emptySpaceMeanIntensity = esSumMean / CFG_CHANNELS_N;
            run("Select None");
        }
    }
    logStep("  - Empty space: area=" + d2s(emptySpaceArea, 0) + " px, mean intensity: " + d2s(emptySpaceMeanIntensity, 2));
    
} else {
    // STATIC MODE: normalize + sum all channels, threshold combined intensity
    logStep("Detecting empty spaces (STATIC: normalize+sum, threshold 0-" + CFG_EMPTY_SPACE_INTENSITY_MAX + ")...");
    
    selectImage("ES_C1");
    run("Enhance Contrast...", "saturated=0 normalize");
    selectImage("ES_C2");
    run("Enhance Contrast...", "saturated=0 normalize");
    selectImage("ES_C3");
    run("Enhance Contrast...", "saturated=0 normalize");
    if (CFG_CHANNELS_N == 4) {
        selectImage("ES_C4");
        run("Enhance Contrast...", "saturated=0 normalize");
    }
    logStep("  - Normalized " + CFG_CHANNELS_N + " median-filtered projections to 0-255");
    
    imageCalculator("Add create", "ES_C1", "ES_C2");
    rename("es_add_12");
    if (CFG_CHANNELS_N == 4) {
        imageCalculator("Add create", "es_add_12", "ES_C3");
        rename("es_add_123");
        safeClose("es_add_12");
        imageCalculator("Add create", "es_add_123", "ES_C4");
        rename("empty_space_combined");
        safeClose("es_add_123");
    } else {
        imageCalculator("Add create", "es_add_12", "ES_C3");
        rename("empty_space_combined");
        safeClose("es_add_12");
    }
    logStep("  - Combined " + CFG_CHANNELS_N + " channels using Add");
    
    // Backup for intensity measurement
    selectImage("empty_space_combined");
    run("Duplicate...", "title=es_gray_backup");
    
    selectImage("empty_space_combined");
    setThreshold(0, CFG_EMPTY_SPACE_INTENSITY_MAX);
    setOption("BlackBackground", true);
    run("Convert to Mask");
    
    // Size filter
    run("Analyze Particles...", "size=" + CFG_EMPTY_SPACE_MIN_SIZE + "-Infinity show=Masks");
    if (isOpen("Mask of empty_space_combined")) {
        selectImage("Mask of empty_space_combined");
        rename("empty_space_mask");
    } else {
        newImage("empty_space_mask", "8-bit black", imgWidth, imgHeight, 1);
    }
    safeClose("empty_space_combined");
    
    // Measure QC on combined backup
    selectImage("empty_space_mask");
    var emptyStatArea, emptyStatMean;
    getStatistics(emptyStatArea, emptyStatMean);
    emptySpaceArea = emptyStatArea * emptyStatMean / 255;
    
    if (emptySpaceArea > 0) {
        selectImage("empty_space_mask");
        run("Create Selection");
        if (selectionType() != -1) {
            selectImage("es_gray_backup");
            run("Restore Selection");
            var esArea, esMean, esMin, esMax, esStd;
            getStatistics(esArea, esMean, esMin, esMax, esStd);
            emptySpaceMeanIntensity = esMean;
            run("Select None");
        }
    }
    safeClose("es_gray_backup");
    logStep("  - Empty space (threshold 0-" + CFG_EMPTY_SPACE_INTENSITY_MAX + "): area=" + d2s(emptySpaceArea, 0) + " px, mean intensity: " + d2s(emptySpaceMeanIntensity, 2));
}

// Close ES grayscale images (no longer needed)
safeClose("ES_C1");
safeClose("ES_C2");
safeClose("ES_C3");
safeClose("ES_C4");

// Create tissue mask (inverse of empty space) for mean intensity measurements
selectImage("empty_space_mask");
run("Duplicate...", "title=tissue_mask");
run("Invert");
logStep("  - Created tissue_mask (excludes empty space) for intensity measurements");

//============================================================================
// PREPROCESS AMYLOID CHANNEL (C2)
//============================================================================
// STEPS 4–5 (Amyloid): median filter (r=1) + rolling-ball background
// subtraction (r=25), then MAX projection (Step 3 again, this time on the
// preprocessed amyloid stack).
logStep("Processing Amyloid channel (C2)...");
selectImage("C2-A1");
run("Median...", "radius=" + CFG_MEDIAN_RADIUS_AMYLOID + " stack");
run("Subtract Background...", "rolling=" + CFG_BG_ROLLING_AMYLOID + " sliding stack");
run("Z Project...", "projection=[Max Intensity]");
rename("MAX_C2-A1");
logStep("  - Median filter: r=" + CFG_MEDIAN_RADIUS_AMYLOID);
logStep("  - Background subtraction: rolling=" + CFG_BG_ROLLING_AMYLOID + ", sliding");
logStep("  - Created MAX projection");

// Measure mean C2 intensity excluding empty space
selectImage("tissue_mask");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("MAX_C2-A1");
    run("Restore Selection");
    getStatistics(_, meanIntensity_C2);
    run("Select None");
} else {
    selectImage("MAX_C2-A1");
    getStatistics(_, meanIntensity_C2);
}
logStep("  - Mean intensity of MAX_C2-A1 (excl. empty): " + d2s(meanIntensity_C2, 2));
safeClose("C2-A1");

//============================================================================
// CREATE CLEAN MAXIMUM COLLAGEN-4 INTENSITY IMAGE
//============================================================================
// STEPS 4–5 (ColIV): median filter (r=2) + rolling-ball background
// subtraction (r=5), then MAX projection (Step 3) on the preprocessed
// ColIV stack.
logStep("Processing Collagen-4 channel (C3)...");
selectImage("C3-A1");
run("Median...", "radius=" + CFG_MEDIAN_RADIUS_COLLAGEN + " stack");
run("Subtract Background...", "rolling=" + CFG_BG_ROLLING_COLLAGEN + " stack");
run("Z Project...", "projection=[Max Intensity]");
var maxProjectionID = getImageID();
logStep("  - Background subtraction: rolling=" + CFG_BG_ROLLING_COLLAGEN);
logStep("  - Created MAX projection");

// Measure mean C3 intensity excluding empty space
selectImage("tissue_mask");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("MAX_C3-A1");
    run("Restore Selection");
    getStatistics(_, meanIntensity_C3);
    run("Select None");
} else {
    selectImage("MAX_C3-A1");
    getStatistics(_, meanIntensity_C3);
}
logStep("  - Mean intensity of MAX_C3-A1 (excl. empty): " + d2s(meanIntensity_C3, 2));

//============================================================================
// CREATE COLLAGEN4 MASK USING AUTO LOCAL THRESHOLD
//============================================================================
// STEP 6: ColIV — Phansalkar adaptive local thresholding (r=50)
logStep("Creating initial collagen4 mask...");
selectImage("MAX_C3-A1");
run("Duplicate...", "title=auto_local_collagen4");
run("Auto Local Threshold", "method=Phansalkar radius=" + CFG_PHANSALKAR_RADIUS_COLLAGEN + " parameter_1=" + CFG_PHANSALKAR_K1_COLLAGEN + " parameter_2=" + CFG_PHANSALKAR_K2_COLLAGEN + " white");
run("Analyze Particles...", "size=" + CFG_COLLAGEN_MIN_SIZE + "-Infinity circularity=0.00-1.00 show=Nothing add");
var initialCollagenROIs = roiManager("count");
logStep("  - Phansalkar: r=" + CFG_PHANSALKAR_RADIUS_COLLAGEN + " k1=" + CFG_PHANSALKAR_K1_COLLAGEN + " k2=" + CFG_PHANSALKAR_K2_COLLAGEN);
logStep("  - Min size filter: " + CFG_COLLAGEN_MIN_SIZE + " px2");
logStep("  - Initial ROIs detected: " + initialCollagenROIs);

//============================================================================
// PROCESS NUCLEI CHANNEL AND CREATE NUCLEI MASK
//============================================================================
// STEPS 4–5 (DAPI): rolling-ball (r=15) + median filter (r=1), then MAX
//                   projection (Step 3) on the preprocessed DAPI stack.
// STEP 7 (DAPI Phansalkar thresholding, r=5): produces the nuclei mask.
//
// v9.0: Wrapped in CFG_BLEEDTHROUGH_CORRECTION toggle. When disabled, the
// DAPI segmentation is skipped entirely and the XS round collagen validation
// falls back to Yen intensity thresholding (same path used when fewer than
// CFG_COLLAGEN_MIN_BLED_NUCLEI bled nuclei are detected).
if (CFG_BLEEDTHROUGH_CORRECTION) {
    logStep("Processing Nuclei channel (C4)...");
    selectImage("C4-A1");
    run("Subtract Background...", "rolling=" + CFG_BG_ROLLING_NUCLEI + " stack");
    run("Median...", "radius=" + CFG_MEDIAN_RADIUS_NUCLEI + " stack");
    run("Z Project...", "projection=[Max Intensity]");
    run("Auto Local Threshold", "method=Phansalkar radius=" + CFG_PHANSALKAR_RADIUS_NUCLEI + " parameter_1=" + CFG_PHANSALKAR_K1_NUCLEI + " parameter_2=" + CFG_PHANSALKAR_K2_NUCLEI + " white");
    run("Analyze Particles...", "size=" + CFG_NUCLEI_MIN_SIZE + "-Infinity pixel show=Masks");
    rename("mask_of_nuclei");
    var dapiMaxID = getImageID();
    logStep("  - Background subtraction: rolling=" + CFG_BG_ROLLING_NUCLEI);
    logStep("  - Phansalkar: r=" + CFG_PHANSALKAR_RADIUS_NUCLEI + " k1=" + CFG_PHANSALKAR_K1_NUCLEI + " k2=" + CFG_PHANSALKAR_K2_NUCLEI);
    logStep("  - Min size filter: " + CFG_NUCLEI_MIN_SIZE + " px2");
} else {
    if (CFG_CHANNELS_N == 3) {
        logStep("DAPI processing skipped (3-channel mode, no DAPI channel available).");
    } else {
        logStep("DAPI bleed-through correction DISABLED.");
        logStep("  - Skipping DAPI channel processing and nuclei segmentation.");
    }
    logStep("  - XS round collagen validation will use Yen intensity fallback.");
}

//============================================================================
// CLEANUP INTERMEDIATE IMAGES
//============================================================================
safeClose("MAX_C4-A1");
safeClose("C4-A1");
safeClose("C3-A1");

//============================================================================
// PHASE 3: ARTIFACT REMOVAL AND VALIDATION (Steps 10–14)
//============================================================================
// Steps 10-13 below process the ColIV ROI set (identification of bleed-
// through, classification, intensity validation, and morphological closing).
// Step 14 (SMA 5-level adaptive intensity filtering) is part of this phase
// in the figure but is implemented later in the code, inside the
// "CREATE SMA MASK" section, because it requires both the validated ColIV
// mask and the initial SMA mask as inputs.

//============================================================================
// FIND BLED-THROUGH NUCLEI
//============================================================================
// STEP 10: Identify DAPI bleed-through into ColIV (nuclear overlap filtering)
// STEPS 11, 12, 13 are also implemented within this section — see submarkers
// at the points where each begins.
//
// v9.0: Entire nuclei-overlap-based filtering is wrapped in
// CFG_BLEEDTHROUGH_CORRECTION toggle. When disabled, skip straight to the
// morphological closing fallback (same path used when no bled nuclei found).
logStep("Identifying bled-through nuclei...");
var totalROIs = roiManager("count");

// Pre-declare these so they're available to downstream XS round logic
// regardless of which branch executes below.
var nArtifacts = 0;
var thresholdValue = 0;

if (totalROIs == 0) {
    logStep("WARNING: No blob-like structures found in auto_local_collagen4. Creating empty collagen4 mask.");
    newImage("collagen4_mask", "8-bit black", getWidth(), getHeight(), 1);
} else {
  if (CFG_BLEEDTHROUGH_CORRECTION) {
    logStep("  - Total ROIs before nuclei filtering: " + totalROIs);
    
    run("Set Measurements...", "mean redirect=None decimal=3");
    selectImage("mask_of_nuclei");
    roiManager("Deselect");
    run("Clear Results");
    roiManager("Multi-Measure measure_all");
    
    var toDelete = newArray(totalROIs);
    var deleteCount = 0;
    var meanIntensity = 0;
    
    for (i = 0; i < totalROIs; i++) {
        meanIntensity = getResult("Mean", i);
        var overlapPercent = (meanIntensity / 255) * 100;
        if (overlapPercent < CFG_NUCLEI_OVERLAP_THRESH) {
            toDelete[deleteCount] = i;
            deleteCount++;
        }
    }
    run("Clear Results");
    
    if (deleteCount > 0) {
        toDelete = Array.trim(toDelete, deleteCount);
        roiManager("select", toDelete);
        roiManager("delete");
    }
    
    var keptROIs = roiManager("count");
    logStep("  - Overlap threshold: " + CFG_NUCLEI_OVERLAP_THRESH + "%");
    logStep("  - ROIs with low nuclei overlap (removed from bled-through pool): " + deleteCount);
    logStep("  - ROIs with high nuclei overlap (classified as bled-through): " + keptROIs);
    
    w = getWidth();
    h = getHeight();
    newImage("PRE_BLED_NUCLEI", "8-bit black", w, h, 1);
    setForegroundColor(255, 255, 255);
    if (roiManager("count") > 0) {
        roiManager("fill");
    }
    run("Select None");
    
    selectImage("PRE_BLED_NUCLEI");
    run("Analyze Particles...", "size=" + CFG_BLED_NUCLEI_MIN_SIZE + "-" + CFG_BLED_NUCLEI_MAX_SIZE + " show=Masks");
    rename("BLED_NUCLEI");
    
    logStep("  - Bled nuclei size filter: " + CFG_BLED_NUCLEI_MIN_SIZE + "-" + CFG_BLED_NUCLEI_MAX_SIZE + " px2");
    
    safeClose("PRE_BLED_NUCLEI");
    roiManager("reset");
    run("Clear Results");
    safeClose("mask_of_nuclei");
    
    //============================================================================
    // CALCULATE 90TH PERCENTILE BRIGHTNESS IN BLED-THROUGH NUCLEI
    //============================================================================
    logStep("Calculating noise threshold from bled-through nuclei...");
    selectImage("BLED_NUCLEI");
    run("Analyze Particles...", "size=0-Infinity pixel show=Nothing add");
    
    nArtifacts = roiManager("count");
    
    if (nArtifacts > 0) {
        logStep("  - Detected " + nArtifacts + " bleed-through nuclei regions");
        
        selectImage("MAX_C3-A1");
        roiManager("Deselect");
        roiManager("Combine");
        
        getHistogram(values, counts, 256);
        
        var totalPixels = 0;
        for (i = 0; i < 256; i++) {
            totalPixels += counts[i];
        }
        
        logStep("  - Total pixels in bled nuclei: " + totalPixels);
        
        var targetPixelCount = floor(totalPixels * CFG_NOISE_PERCENTILE);
        var cumulativeCount = 0;
        
        for (i = 0; i < 256; i++) {
            cumulativeCount += counts[i];
            if (cumulativeCount >= targetPixelCount) {
                thresholdValue = i;
                break;
            }
        }
        
        logStep("  - " + (CFG_NOISE_PERCENTILE * 100) + "th percentile intensity: " + thresholdValue);
        
        run("Select None");
        run("Duplicate...", "title=pre_clean_collagen4_mask");
        var statsArea, statsMean, statsMin, statsMax;
        getStatistics(statsArea, statsMean, statsMin, statsMax);
        logStep("  - Image intensity range: " + statsMin + " - " + statsMax);
        
        setThreshold(thresholdValue, 255);
        setOption("BlackBackground", true);
        run("Convert to Mask");
        run("Morphological Filters", "operation=Closing element=Diamond radius=" + CFG_CLOSING_RADIUS_SMALL);
        run("Analyze Particles...", "size=" + CFG_COLLAGEN_MIN_SIZE + "-Infinity show=Masks");
        rename("clean_collagen4_mask");
        
        logStep("  - Applied threshold: " + thresholdValue + "-255");
        logStep("  - Morphological closing radius: " + CFG_CLOSING_RADIUS_SMALL);
        
    } else {
        logStep("  - No bleed-through nuclei found. Skipping intensity filtering.");
        thresholdValue = 0;
        selectImage("auto_local_collagen4");
        run("Morphological Filters", "operation=Closing element=Diamond radius=" + CFG_CLOSING_RADIUS_LARGE);
        rename("clean_collagen4_mask");
    }
  } else {
    // CFG_BLEEDTHROUGH_CORRECTION is false: use fallback morphological closing
    // of auto_local_collagen4 to produce clean_collagen4_mask. Leave nArtifacts
    // and thresholdValue at 0 so the XS round logic uses Yen intensity fallback.
    logStep("  - Total ROIs before filtering: " + totalROIs);
    logStep("  - Skipping nuclei-overlap filtering (correction disabled).");
    selectImage("auto_local_collagen4");
    run("Morphological Filters", "operation=Closing element=Diamond radius=" + CFG_CLOSING_RADIUS_LARGE);
    rename("clean_collagen4_mask");
    logStep("  - clean_collagen4_mask built via morphological closing only (radius=" + CFG_CLOSING_RADIUS_LARGE + ")");
  }
    
    //============================================================================
    // QUALITY CONTROL: Filter ROIs by overlap with clean mask
    //============================================================================
    // STEP 11: Size/shape classification of ColIV ROIs (5 sizes × 4 shapes)
    // STEP 12: Yen intensity validation per category
    // Both happen within the loop below: each ROI is classified into one of
    // 17 size×shape categories (Step 11), then validated either by spatial
    // overlap with the cleaned mask (XS round only) or by Yen intensity
    // threshold with the per-category factor (Step 12).
    logStep("Validating collagen ROIs against clean mask...");
    
    roiManager("reset");
    selectImage("auto_local_collagen4");
    run("Analyze Particles...", "size=" + CFG_COLLAGEN_MIN_SIZE + "-Infinity circularity=0.00-1.00 show=Nothing add composite");
    
    var totalROIs = roiManager("count");
    if (totalROIs == 0) {
        logStep("WARNING: No ROIs found in auto_local_collagen4. Creating empty mask.");
        newImage("collagen4_mask", "8-bit black", getWidth(), getHeight(), 1);
    } else {
        logStep("  - Total ROIs before filtering: " + totalROIs);
        
        // Build category lookup arrays (17 categories: 4 sizes × 4 shapes + XL)
        // Index: sizeCode*4 + shapeCode for sizes 0-3; index 16 for XL
        // sizeCode: 0=XS, 1=S, 2=M, 3=L, 4=XL
        // shapeCode: 0=round, 1=ambiguous, 2=elongated, 3=linear
        var catNames = newArray("XS_round","XS_ambig","XS_elong","XS_linear","S_round","S_ambig","S_elong","S_linear","M_round","M_ambig","M_elong","M_linear","L_round","L_ambig","L_elong","L_linear","XL");
        var yenFactors = newArray(17);
        yenFactors[0]=CFG_YEN_XS_ROUND; yenFactors[1]=CFG_YEN_XS_AMBIG; yenFactors[2]=CFG_YEN_XS_ELONG; yenFactors[3]=CFG_YEN_XS_LINEAR;
        yenFactors[4]=CFG_YEN_S_ROUND;  yenFactors[5]=CFG_YEN_S_AMBIG;  yenFactors[6]=CFG_YEN_S_ELONG;  yenFactors[7]=CFG_YEN_S_LINEAR;
        yenFactors[8]=CFG_YEN_M_ROUND;  yenFactors[9]=CFG_YEN_M_AMBIG;  yenFactors[10]=CFG_YEN_M_ELONG; yenFactors[11]=CFG_YEN_M_LINEAR;
        yenFactors[12]=CFG_YEN_L_ROUND; yenFactors[13]=CFG_YEN_L_AMBIG; yenFactors[14]=CFG_YEN_L_ELONG; yenFactors[15]=CFG_YEN_L_LINEAR;
        yenFactors[16]=CFG_YEN_XL;
        
        var sizeBounds = newArray(CFG_COLLAGEN_XS_MAX, CFG_COLLAGEN_S_MAX, CFG_COLLAGEN_M_MAX, CFG_COLLAGEN_L_MAX);
        var circLinearBounds = newArray(CFG_CIRC_XS_LINEAR, CFG_CIRC_S_LINEAR, CFG_CIRC_M_LINEAR, CFG_CIRC_L_LINEAR);
        var circElongBounds = newArray(CFG_CIRC_XS_ELONG, CFG_CIRC_S_ELONG, CFG_CIRC_M_ELONG, CFG_CIRC_L_ELONG);
        var circRoundBounds = newArray(CFG_CIRC_XS_ROUND, CFG_CIRC_S_ROUND, CFG_CIRC_M_ROUND, CFG_CIRC_L_ROUND);
        
        var catPass = newArray(17);
        var catFail = newArray(17);
        for (var ci2 = 0; ci2 < 17; ci2++) { catPass[ci2] = 0; catFail[ci2] = 0; }
        
        // Measure ROI areas and circularities
        run("Set Measurements...", "area shape redirect=None decimal=3");
        run("Clear Results");
        selectImage("auto_local_collagen4");
        for (i = 0; i < totalROIs; i++) {
            roiManager("select", i);
            roiManager("measure");
        }
        var roiAreas = newArray(maxOf(totalROIs, 1));
        var roiCircs = newArray(maxOf(totalROIs, 1));
        for (i = 0; i < totalROIs; i++) {
            roiAreas[i] = getResult("Area", i);
            roiCircs[i] = getResult("Circ.", i);
        }
        run("Clear Results");
        run("Select None");
        
        // Classify ROIs into post-filter (XS+S) and pre-filter (M+L+XL) groups
        var postIndices = newArray(maxOf(totalROIs, 1));
        var postCount = 0;
        var preIndices = newArray(maxOf(totalROIs, 1));
        var preCount = 0;
        for (i = 0; i < totalROIs; i++) {
            if (roiAreas[i] <= CFG_COLLAGEN_S_MAX) {
                postIndices[postCount] = i;
                postCount++;
            } else {
                preIndices[preCount] = i;
                preCount++;
            }
        }
        if (postCount > 0) postIndices = Array.trim(postIndices, postCount);
        if (preCount > 0) preIndices = Array.trim(preIndices, preCount);
        
        // Count per size for log
        var nXS = 0; var nS = 0; var nM = 0; var nL = 0; var nXL = 0;
        for (i = 0; i < totalROIs; i++) {
            if (roiAreas[i] <= CFG_COLLAGEN_XS_MAX) nXS++;
            else if (roiAreas[i] <= CFG_COLLAGEN_S_MAX) nS++;
            else if (roiAreas[i] <= CFG_COLLAGEN_M_MAX) nM++;
            else if (roiAreas[i] <= CFG_COLLAGEN_L_MAX) nL++;
            else nXL++;
        }
        logStep("  - XS(<=" + CFG_COLLAGEN_XS_MAX + "): " + nXS + " | S(<=" + CFG_COLLAGEN_S_MAX + "): " + nS + " | M(<=" + CFG_COLLAGEN_M_MAX + "): " + nM + " | L(<=" + CFG_COLLAGEN_L_MAX + "): " + nL + " | XL(>" + CFG_COLLAGEN_L_MAX + "): " + nXL);
        logStep("  - Post-filter group (XS+S): " + postCount + " | Pre-filter group (M+L+XL): " + preCount);
        
        //=== COMPUTE YEN MEAN INTENSITY ===
        logStep("Computing Yen mean intensity on C3...");
        roiManager("reset");
        
        selectImage("MAX_C3-A1");
        run("Duplicate...", "title=Yen_Collagen4");
        setAutoThreshold("Yen dark");
        setOption("BlackBackground", true);
        run("Convert to Mask");
        
        selectImage("Yen_Collagen4");
        run("Analyze Particles...", "size=0-Infinity show=Nothing add");
        var nYenC4ROIs = roiManager("count");
        var mean_intensity_yen_c4 = 0;
        
        if (nYenC4ROIs > 0) {
            selectImage("MAX_C3-A1");
            roiManager("Deselect");
            roiManager("Combine");
            var yenC4Area, yenC4Min, yenC4Max;
            getStatistics(yenC4Area, mean_intensity_yen_c4, yenC4Min, yenC4Max);
            run("Select None");
            logStep("  - Yen mask: " + nYenC4ROIs + " regions, mean intensity=" + d2s(mean_intensity_yen_c4, 2));
        } else {
            logStep("  - WARNING: No Yen ROIs on C3. Keeping all objects.");
        }
        safeClose("Yen_Collagen4");
        
        // Re-extract collagen ROIs
        roiManager("reset");
        selectImage("auto_local_collagen4");
        run("Analyze Particles...", "size=" + CFG_COLLAGEN_MIN_SIZE + "-Infinity circularity=0.00-1.00 show=Nothing add composite");
        
        //================================================================
        //=== POST-FILTER GROUP (XS + S): filter then close ===
        //================================================================
        logStep("Filtering post-filter group (XS+S)...");
        
        // Determine if XS round uses overlap or Yen backup
        // v9.0: useOverlap requires both that the correction is enabled AND that
        // enough bled nuclei were detected to calibrate the overlap threshold.
        // When correction is disabled, nArtifacts=0 so the overlap method
        // wouldn't be reachable anyway — the explicit check on
        // CFG_BLEEDTHROUGH_CORRECTION in the log messages is just for clarity.
        var useOverlap = (nArtifacts >= CFG_COLLAGEN_MIN_BLED_NUCLEI);
        if (useOverlap) {
            logStep("  - XS round: overlap method (bled nuclei=" + nArtifacts + " >= " + CFG_COLLAGEN_MIN_BLED_NUCLEI + ")");
            // Multi-Measure on clean_collagen4_mask for overlap values
            run("Set Measurements...", "mean redirect=None decimal=3");
            selectImage("clean_collagen4_mask");
            roiManager("Deselect");
            run("Clear Results");
            roiManager("Multi-Measure measure_all");
        } else if (CFG_BLEEDTHROUGH_CORRECTION) {
            logStep("  - XS round: Yen BACKUP (bled nuclei=" + nArtifacts + " < " + CFG_COLLAGEN_MIN_BLED_NUCLEI + ")");
        } else {
            logStep("  - XS round: Yen BACKUP (bleed-through correction disabled)");
        }
        
        var postPassIndices = newArray(maxOf(postCount, 1));
        var postPassCount = 0;
        
        for (var pi = 0; pi < postCount; pi++) {
            var idx = postIndices[pi];
            var area_val = roiAreas[idx];
            var circ_val = roiCircs[idx];
            
            // Determine size code (0=XS, 1=S)
            var sizeCode = 1;
            if (area_val <= CFG_COLLAGEN_XS_MAX) { sizeCode = 0; }
            
            // Determine shape code (0=round, 1=ambig, 2=elong, 3=linear)
            var shapeCode = 3;
            if (circ_val >= circRoundBounds[sizeCode]) { shapeCode = 0; }
            else if (circ_val >= circElongBounds[sizeCode]) { shapeCode = 1; }
            else if (circ_val >= circLinearBounds[sizeCode]) { shapeCode = 2; }
            
            var catIdx = sizeCode * 4 + shapeCode;
            var passes = 0;
            
            // XS round with enough bled nuclei → overlap method
            if (sizeCode == 0 && shapeCode == 0 && useOverlap) {
                meanIntensity = getResult("Mean", idx);
                var overlapPct = (meanIntensity / 255) * 100;
                if (overlapPct >= CFG_COLLAGEN_OVERLAP_THRESH) { passes = 1; }
            } else {
                // Yen intensity filter
                if (mean_intensity_yen_c4 > 0) {
                    var yenThresh = mean_intensity_yen_c4 * yenFactors[catIdx];
                    selectImage("MAX_C3-A1");
                    roiManager("select", idx);
                    var pfArea, pfMean, pfMin, pfMax, pfStd;
                    getStatistics(pfArea, pfMean, pfMin, pfMax, pfStd);
                    if (pfMax >= yenThresh) { passes = 1; }
                } else {
                    passes = 1;  // No Yen threshold — keep
                }
            }
            
            if (passes) {
                postPassIndices[postPassCount] = idx;
                postPassCount++;
                catPass[catIdx] = catPass[catIdx] + 1;
            } else {
                catFail[catIdx] = catFail[catIdx] + 1;
            }
        }
        if (useOverlap) { run("Clear Results"); }
        run("Select None");
        
        // Log post-filter results
        for (var ci2 = 0; ci2 < 8; ci2++) {
            if (catPass[ci2] + catFail[ci2] > 0) {
                var methodStr = "Yen*" + d2s(yenFactors[ci2], 2);
                if (ci2 == 0 && useOverlap) { methodStr = "overlap>=" + CFG_COLLAGEN_OVERLAP_THRESH + "%"; }
                logStep("  - " + catNames[ci2] + ": pass=" + catPass[ci2] + " fail=" + catFail[ci2] + " (" + methodStr + ")");
            }
        }
        
        //=== Create post-filter mask (XS+S survivors) → AND → close ===
        // STEP 13: Morphological closing to reconnect fragments
        // Applied here to XS+S after per-bin Yen filtering (post-filter group),
        // and below at line ~1980 to M+L+XL before per-bin filtering
        // (pre-filter group). See Section 6 (f) of the CFG block for the
        // rationale behind the post-vs-pre ordering.
        w = getWidth();
        h = getHeight();
        newImage("postfilter_collagen4_raw", "8-bit black", w, h, 1);
        if (postPassCount > 0) {
            postPassIndices = Array.trim(postPassIndices, postPassCount);
            if (postPassCount > 1) {
                roiManager("select", postPassIndices);
                roiManager("Combine");
            } else {
                roiManager("select", postPassIndices[0]);
            }
            setForegroundColor(255, 255, 255);
            run("Fill", "slice");
            run("Select None");
            imageCalculator("AND", "postfilter_collagen4_raw", "auto_local_collagen4");
        }
        selectImage("postfilter_collagen4_raw");
        run("Morphological Filters", "operation=Closing element=Diamond radius=" + CFG_COLLAGEN_CLOSING_POST);
        rename("postfilter_collagen4_mask");
        safeClose("postfilter_collagen4_raw");
        logStep("  - Post-filter mask: " + postPassCount + " survivors, closed(r=" + CFG_COLLAGEN_CLOSING_POST + ")");
        
        //================================================================
        //=== CREATE PRE-FILTER RAW MASK (M+L+XL) while original ROIs still in manager ===
        //================================================================
        logStep("Creating pre-filter raw mask (M+L+XL)...");
        newImage("prefilter_collagen4_raw", "8-bit black", w, h, 1);
        if (preCount > 0) {
            if (preCount > 1) {
                roiManager("select", preIndices);
                roiManager("Combine");
            } else {
                roiManager("select", preIndices[0]);
            }
            setForegroundColor(255, 255, 255);
            run("Fill", "slice");
            run("Select None");
            imageCalculator("AND", "prefilter_collagen4_raw", "auto_local_collagen4");
        }
        
        //================================================================
        //=== PRE-FILTER GROUP (M+L+XL): close → re-extract → filter ===
        //================================================================
        logStep("Processing pre-filter group (M+L+XL)...");
        selectImage("prefilter_collagen4_raw");
        run("Morphological Filters", "operation=Closing element=Diamond radius=" + CFG_COLLAGEN_CLOSING_PRE);
        rename("prefilter_collagen4_closed");
        safeClose("prefilter_collagen4_raw");
        logStep("  - Pre-filter mask closed (radius=" + CFG_COLLAGEN_CLOSING_PRE + ")");
        
        // Re-extract ROIs from closed mask
        roiManager("reset");
        selectImage("prefilter_collagen4_closed");
        run("Analyze Particles...", "size=0-Infinity show=Nothing add composite");
        var nPreClosedROIs = roiManager("count");
        logStep("  - Pre-filter closed ROIs: " + nPreClosedROIs);
        
        // Measure area + circularity of re-extracted ROIs
        var preAreas = newArray(maxOf(nPreClosedROIs, 1));
        var preCircs = newArray(maxOf(nPreClosedROIs, 1));
        if (nPreClosedROIs > 0) {
            run("Set Measurements...", "area shape redirect=None decimal=3");
            run("Clear Results");
            selectImage("prefilter_collagen4_closed");
            for (var mi = 0; mi < nPreClosedROIs; mi++) {
                roiManager("select", mi);
                roiManager("measure");
            }
            for (var mi = 0; mi < nPreClosedROIs; mi++) {
                preAreas[mi] = getResult("Area", mi);
                preCircs[mi] = getResult("Circ.", mi);
            }
            run("Clear Results");
        }
        
        // Filter each re-extracted ROI by Yen with category-appropriate factor
        var preToDelete = newArray(maxOf(nPreClosedROIs, 1));
        var preDeleteCount = 0;
        
        if (nPreClosedROIs > 0 && mean_intensity_yen_c4 > 0) {
            for (var mi = 0; mi < nPreClosedROIs; mi++) {
                // Determine size code (2=M, 3=L, 4=XL)
                var preSizeCode = 4;
                if (preAreas[mi] <= CFG_COLLAGEN_M_MAX) { preSizeCode = 2; }
                else if (preAreas[mi] <= CFG_COLLAGEN_L_MAX) { preSizeCode = 3; }
                
                // Determine category index
                var preCatIdx = 16;  // XL default
                if (preSizeCode < 4) {
                    var preShapeCode = 3;  // linear default
                    if (preCircs[mi] >= circRoundBounds[preSizeCode]) { preShapeCode = 0; }
                    else if (preCircs[mi] >= circElongBounds[preSizeCode]) { preShapeCode = 1; }
                    else if (preCircs[mi] >= circLinearBounds[preSizeCode]) { preShapeCode = 2; }
                    preCatIdx = preSizeCode * 4 + preShapeCode;
                }
                
                var preYenThresh = mean_intensity_yen_c4 * yenFactors[preCatIdx];
                
                selectImage("MAX_C3-A1");
                roiManager("select", mi);
                var prArea, prMean, prMin, prMax, prStd;
                getStatistics(prArea, prMean, prMin, prMax, prStd);
                
                if (prMax < preYenThresh) {
                    preToDelete[preDeleteCount] = mi;
                    preDeleteCount++;
                    catFail[preCatIdx] = catFail[preCatIdx] + 1;
                } else {
                    catPass[preCatIdx] = catPass[preCatIdx] + 1;
                }
            }
            run("Select None");
            
            if (preDeleteCount > 0) {
                preToDelete = Array.trim(preToDelete, preDeleteCount);
                roiManager("select", preToDelete);
                roiManager("delete");
            }
        }
        
        var prePassCount = roiManager("count");
        
        // Log pre-filter results
        for (var ci2 = 8; ci2 < 17; ci2++) {
            if (catPass[ci2] + catFail[ci2] > 0) {
                logStep("  - " + catNames[ci2] + ": pass=" + catPass[ci2] + " fail=" + catFail[ci2] + " (Yen*" + d2s(yenFactors[ci2], 2) + ")");
            }
        }
        
        // Create filtered pre-filter mask
        newImage("prefilter_collagen4_mask", "8-bit black", w, h, 1);
        if (prePassCount > 0) {
            roiManager("Deselect");
            roiManager("Combine");
            setForegroundColor(255, 255, 255);
            run("Fill", "slice");
            run("Select None");
            imageCalculator("AND", "prefilter_collagen4_mask", "prefilter_collagen4_closed");
        }
        safeClose("prefilter_collagen4_closed");
        
        //================================================================
        //=== COMBINE post-filter + pre-filter → collagen4_mask ===
        //================================================================
        imageCalculator("OR create", "postfilter_collagen4_mask", "prefilter_collagen4_mask");
        rename("collagen4_mask");
        safeClose("postfilter_collagen4_mask");
        safeClose("prefilter_collagen4_mask");
        
        // Log final summary
        var finalColArea, finalColMean;
        selectImage("collagen4_mask");
        getStatistics(finalColArea, finalColMean);
        var finalColPixels = finalColArea * finalColMean / 255;
        var totalPass = postPassCount + prePassCount;
        logStep("  - Final collagen4_mask: post(" + postPassCount + ") + pre(" + prePassCount + ") = " + totalPass + ", area=" + d2s(finalColPixels, 0) + " px");
    }
}

//============================================================================
// CLEANUP
//============================================================================
run("Clear Results");
roiManager("reset");
safeClose("BLED_NUCLEI");
safeClose("pre_collagen4_mask-Closing");
safeClose("pre_collagen4_mask");
safeClose("pre_clean_collagen4_mask");
safeClose("clean_collagen4_mask");
safeClose("pre_clean_collagen4_mask-Closing");
safeClose("auto_local_collagen4-Closing");
safeClose("mask_of_nuclei");
safeClose("postfilter_collagen4_raw");
safeClose("postfilter_collagen4_raw-Closing");
safeClose("postfilter_collagen4_mask");
safeClose("prefilter_collagen4_raw");
safeClose("prefilter_collagen4_raw-Closing");
safeClose("prefilter_collagen4_closed");
safeClose("prefilter_collagen4_mask");
safeClose("Yen_Collagen4");
roiManager("reset");

//============================================================================
// CREATE SMA MASK
//============================================================================
// STEPS 4–5 (SMA): rolling-ball (r=25) + median filter (r=2), then MAX
//                  projection (Step 3) on the preprocessed SMA stack.
// STEP 8: SMA — Phansalkar thresholding (r=15) + diamond closing (r=2),
//         producing the initial SMA mask. (Phansalkar call follows below.)
// STEP 14: SMA 5-level adaptive intensity filtering (ColIV overlap-dependent
//          thresholds). Implemented inside this section because it requires
//          both the initial SMA mask (just produced) and the validated ColIV
//          mask (from Phase 3) as inputs. The 5-level classification logic
//          and per-level threshold application appear further down within
//          this section.
logStep("Processing SMA channel (C1)...");
selectImage("C1-A1");
run("Subtract Background...", "rolling=" + CFG_BG_ROLLING_SMA + " stack");
run("Median...", "radius=" + CFG_MEDIAN_RADIUS_SMA + " stack");
run("Z Project...", "projection=[Max Intensity]");
var maxProjectionSMA = getImageID();
logStep("  - Background subtraction: rolling=" + CFG_BG_ROLLING_SMA);
logStep("  - Created MAX projection");

// Measure mean C1 intensity excluding empty space (for QC and adaptive threshold)
selectImage("tissue_mask");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("MAX_C1-A1");
    run("Restore Selection");
    getStatistics(_, meanIntensity_C1);
    run("Select None");
} else {
    selectImage("MAX_C1-A1");
    getStatistics(_, meanIntensity_C1);
}
logStep("  - Mean intensity of MAX_C1-A1 (excl. empty): " + d2s(meanIntensity_C1, 2));
safeClose("tissue_mask");

// Compute per-level adaptive SMA threshold factors
if (CFG_ADAPTIVE_SMA_ENABLED) {
    sma_factor_L5 = computeAdaptiveFactor(meanIntensity_C1, CFG_ADAPTIVE_L5_INTENSITY_LOW, CFG_ADAPTIVE_L5_INTENSITY_HIGH, CFG_ADAPTIVE_L5_FACTOR_MIN, CFG_ADAPTIVE_L5_FACTOR_MAX);
    sma_factor_L4 = computeAdaptiveFactor(meanIntensity_C1, CFG_ADAPTIVE_L4_INTENSITY_LOW, CFG_ADAPTIVE_L4_INTENSITY_HIGH, CFG_ADAPTIVE_L4_FACTOR_MIN, CFG_ADAPTIVE_L4_FACTOR_MAX);
    sma_factor_L3 = computeAdaptiveFactor(meanIntensity_C1, CFG_ADAPTIVE_L3_INTENSITY_LOW, CFG_ADAPTIVE_L3_INTENSITY_HIGH, CFG_ADAPTIVE_L3_FACTOR_MIN, CFG_ADAPTIVE_L3_FACTOR_MAX);
    sma_factor_L2 = computeAdaptiveFactor(meanIntensity_C1, CFG_ADAPTIVE_L2_INTENSITY_LOW, CFG_ADAPTIVE_L2_INTENSITY_HIGH, CFG_ADAPTIVE_L2_FACTOR_MIN, CFG_ADAPTIVE_L2_FACTOR_MAX);
    sma_factor_L1 = computeAdaptiveFactor(meanIntensity_C1, CFG_ADAPTIVE_L1_INTENSITY_LOW, CFG_ADAPTIVE_L1_INTENSITY_HIGH, CFG_ADAPTIVE_L1_FACTOR_MIN, CFG_ADAPTIVE_L1_FACTOR_MAX);
    logStep("  - ADAPTIVE SMA factors (meanC1=" + d2s(meanIntensity_C1, 2) + "):");
    logStep("    L5=" + d2s(sma_factor_L5, 4) + " L4=" + d2s(sma_factor_L4, 4) + " L3=" + d2s(sma_factor_L3, 4) + " L2=" + d2s(sma_factor_L2, 4) + " L1=" + d2s(sma_factor_L1, 4));
} else {
    sma_factor_L5 = CFG_SMA_THRESHOLD_FACTOR;
    sma_factor_L4 = CFG_SMA_THRESHOLD_FACTOR;
    sma_factor_L3 = CFG_SMA_THRESHOLD_FACTOR;
    sma_factor_L2 = CFG_SMA_THRESHOLD_FACTOR;
    sma_factor_L1 = CFG_SMA_THRESHOLD_FACTOR;
    logStep("  - Using fixed SMA factor for all levels: " + d2s(CFG_SMA_THRESHOLD_FACTOR, 4));
}

logStep("Creating Yen threshold for SMA intensity calibration...");
selectImage("MAX_C1-A1");
run("Duplicate...", "title=Yen_SMA");
setAutoThreshold("Yen dark");
setOption("BlackBackground", true);
run("Convert to Mask");

roiManager("reset");
selectImage("Yen_SMA");
run("Analyze Particles...", "size=0-Infinity show=Nothing add");

var nYenROIs = roiManager("count");
var mean_intensity_yen = 0;

if (nYenROIs > 0) {
    logStep("  - Yen mask detected " + nYenROIs + " bright SMA regions");
    
    selectImage("MAX_C1-A1");
    roiManager("Deselect");
    roiManager("Combine");
    
    var yenArea, yenMin, yenMax;
    getStatistics(yenArea, mean_intensity_yen, yenMin, yenMax);
    
    // Compute per-level intensity thresholds
    sma_thresh_L5 = mean_intensity_yen * sma_factor_L5;
    sma_thresh_L4 = mean_intensity_yen * sma_factor_L4;
    sma_thresh_L3 = mean_intensity_yen * sma_factor_L3;
    sma_thresh_L2 = mean_intensity_yen * sma_factor_L2;
    sma_thresh_L1 = mean_intensity_yen * sma_factor_L1;
    
    logStep("  - Yen Mean: " + d2s(mean_intensity_yen, 2));
    logStep("  - Intensity thresholds: L5=" + d2s(sma_thresh_L5, 2) + " L4=" + d2s(sma_thresh_L4, 2) + " L3=" + d2s(sma_thresh_L3, 2) + " L2=" + d2s(sma_thresh_L2, 2) + " L1=" + d2s(sma_thresh_L1, 2));
    
    run("Select None");
} else {
    logStep("WARNING: No Yen ROIs detected. Using default threshold of 0 for all levels.");
}

logStep("Creating Auto Local Threshold mask for SMA...");
roiManager("reset");
selectImage("MAX_C1-A1");
run("Duplicate...", "title=auto_local_SMA");
run("Auto Local Threshold", "method=Phansalkar radius=" + CFG_PHANSALKAR_RADIUS_SMA + " parameter_1=" + CFG_PHANSALKAR_K1_SMA + " parameter_2=" + CFG_PHANSALKAR_K2_SMA + " white");

// Morphological closing bridges fragmented vessel wall segments before ROI extraction.
selectImage("auto_local_SMA");
run("Duplicate...", "title=SMA_for_closing");
run("Morphological Filters", "operation=Closing element=Diamond radius=" + CFG_SMA_CLOSING_RADIUS);
rename("auto_local_SMA_closed");
safeClose("SMA_for_closing");
logStep("  - SMA closing: Diamond r=" + CFG_SMA_CLOSING_RADIUS);

selectImage("auto_local_SMA_closed");
run("Analyze Particles...", "size=" + CFG_SMA_MIN_SIZE + "-Infinity show=Nothing add composite");

var nAutoROIs = roiManager("count");
logStep("  - Phansalkar: r=" + CFG_PHANSALKAR_RADIUS_SMA + " k1=" + CFG_PHANSALKAR_K1_SMA + " k2=" + CFG_PHANSALKAR_K2_SMA);
logStep("  - Min size filter: " + CFG_SMA_MIN_SIZE + " px2");
logStep("  - Auto Local SMA ROIs before level-based filtering: " + nAutoROIs);

// --- 5-LEVEL COLLAGEN OVERLAP + INTENSITY FILTER ---
// Each ROI is classified by collagen overlap percentage (L5→L1, highest first),
// then filtered by the level-specific intensity threshold applied to max pixel.
var sma_total_removed = 0;

if (nAutoROIs > 0 && isOpen("auto_local_collagen4")) {
    logStep("Classifying SMA ROIs by collagen-4 overlap level and filtering by intensity...");
    logStep("  - Overlap levels: L5>=" + CFG_SMA_OL_LEVEL5 + "% L4>=" + CFG_SMA_OL_LEVEL4 + "% L3>=" + CFG_SMA_OL_LEVEL3 + "% L2>=" + CFG_SMA_OL_LEVEL2 + "% L1>=" + CFG_SMA_OL_LEVEL1 + "%");
    
    var toDelete = newArray(nAutoROIs);
    var deleteCount = 0;
    var olArea, olMean, olMin, olMax, olStd;
    var roiArea, roiMean, roiMin, maxPixel, roiStd;
    
    // Track per-level stats
    var classified_L1 = 0; var classified_L2 = 0; var classified_L3 = 0; var classified_L4 = 0; var classified_L5 = 0;
    
    for (var sma_i = 0; sma_i < nAutoROIs; sma_i++) {
        // Measure collagen overlap
        selectImage("auto_local_collagen4");
        roiManager("select", sma_i);
        getStatistics(olArea, olMean, olMin, olMax, olStd);
        var overlapPct = (olMean / 255) * 100;
        
        // Measure max pixel intensity on SMA channel
        selectImage("MAX_C1-A1");
        roiManager("select", sma_i);
        getStatistics(roiArea, roiMean, roiMin, maxPixel, roiStd);
        
        // Classify into level (highest first) and apply level threshold
        var levelThresh = 0;
        var levelName = "";
        if (overlapPct >= CFG_SMA_OL_LEVEL5) {
            levelThresh = sma_thresh_L5; levelName = "L5"; classified_L5++;
        } else if (overlapPct >= CFG_SMA_OL_LEVEL4) {
            levelThresh = sma_thresh_L4; levelName = "L4"; classified_L4++;
        } else if (overlapPct >= CFG_SMA_OL_LEVEL3) {
            levelThresh = sma_thresh_L3; levelName = "L3"; classified_L3++;
        } else if (overlapPct >= CFG_SMA_OL_LEVEL2) {
            levelThresh = sma_thresh_L2; levelName = "L2"; classified_L2++;
        } else {
            levelThresh = sma_thresh_L1; levelName = "L1"; classified_L1++;
        }
        
        // Apply intensity filter for this level
        if (levelThresh > 0 && maxPixel < levelThresh) {
            toDelete[deleteCount] = sma_i;
            deleteCount++;
            // Track removal per level
            if (overlapPct >= CFG_SMA_OL_LEVEL5) { sma_removed_L5++; }
            else if (overlapPct >= CFG_SMA_OL_LEVEL4) { sma_removed_L4++; }
            else if (overlapPct >= CFG_SMA_OL_LEVEL3) { sma_removed_L3++; }
            else if (overlapPct >= CFG_SMA_OL_LEVEL2) { sma_removed_L2++; }
            else { sma_removed_L1++; }
        }
    }
    
    run("Select None");
    
    // Delete failed ROIs
    if (deleteCount > 0) {
        toDelete = Array.trim(toDelete, deleteCount);
        roiManager("select", toDelete);
        roiManager("delete");
    }
    
    sma_total_removed = deleteCount;
    
    // Compute kept counts
    sma_kept_L5 = classified_L5 - sma_removed_L5;
    sma_kept_L4 = classified_L4 - sma_removed_L4;
    sma_kept_L3 = classified_L3 - sma_removed_L3;
    sma_kept_L2 = classified_L2 - sma_removed_L2;
    sma_kept_L1 = classified_L1 - sma_removed_L1;
    
    logStep("  - Classification: L5=" + classified_L5 + " L4=" + classified_L4 + " L3=" + classified_L3 + " L2=" + classified_L2 + " L1=" + classified_L1);
    logStep("  - Removed by intensity: L5=" + sma_removed_L5 + " L4=" + sma_removed_L4 + " L3=" + sma_removed_L3 + " L2=" + sma_removed_L2 + " L1=" + sma_removed_L1 + " (total=" + sma_total_removed + ")");
    logStep("  - Kept: L5=" + sma_kept_L5 + " L4=" + sma_kept_L4 + " L3=" + sma_kept_L3 + " L2=" + sma_kept_L2 + " L1=" + sma_kept_L1 + " (total=" + roiManager("count") + ")");
    
} else if (!isOpen("auto_local_collagen4")) {
    // No collagen mask available — use L1 threshold for all ROIs
    logStep("WARNING: auto_local_collagen4 not available. Applying L1 threshold to all ROIs.");
    if (nAutoROIs > 0 && sma_thresh_L1 > 0) {
        var toDelete = newArray(nAutoROIs);
        var deleteCount = 0;
        var roiArea, roiMean, roiMin, maxPixel, roiStd;
        for (var sma_i = 0; sma_i < nAutoROIs; sma_i++) {
            selectImage("MAX_C1-A1");
            roiManager("select", sma_i);
            getStatistics(roiArea, roiMean, roiMin, maxPixel, roiStd);
            if (maxPixel < sma_thresh_L1) {
                toDelete[deleteCount] = sma_i;
                deleteCount++;
            }
        }
        run("Select None");
        if (deleteCount > 0) {
            toDelete = Array.trim(toDelete, deleteCount);
            roiManager("select", toDelete);
            roiManager("delete");
        }
        sma_total_removed = deleteCount;
        sma_kept_L1 = roiManager("count");
        logStep("  - Applied L1 threshold: removed=" + deleteCount + " kept=" + sma_kept_L1);
    }
}

// Create SMA mask from surviving ROIs
sma_roi_count_for_qc = roiManager("count");

if (sma_roi_count_for_qc > 0) {
    sma_areas_count = sma_roi_count_for_qc;
    sma_areas_for_qc = newArray(sma_areas_count);
    run("Set Measurements...", "area redirect=None decimal=3");
    var tmpArea;
    for (var ai = 0; ai < sma_areas_count; ai++) {
        roiManager("select", ai);
        getStatistics(tmpArea);
        sma_areas_for_qc[ai] = tmpArea;
    }
    run("Select None");
    
    w = getWidth();
    h = getHeight();
    newImage("SMA_selection_mask", "8-bit black", w, h, 1);
    roiManager("Deselect");
    roiManager("Combine");
    setForegroundColor(255, 255, 255);
    run("Fill", "slice");
    run("Select None");
    imageCalculator("AND create", "auto_local_SMA", "SMA_selection_mask");
    rename("SMA_mask");
    safeClose("SMA_selection_mask");
    logStep("  - Created SMA_mask with " + sma_roi_count_for_qc + " validated ROIs (AND with threshold to preserve holes)");
} else {
    newImage("SMA_mask", "8-bit black", getWidth(), getHeight(), 1);
    logStep("WARNING: No SMA ROIs survived filtering. Created empty SMA_mask.");
    sma_roi_count_for_qc = 0;
    sma_areas_count = 0;
}

//============================================================================
// CLEANUP
//============================================================================
run("Clear Results");
roiManager("reset");
safeClose("Yen_SMA");
safeClose("C1-A1");

//============================================================================
// v8.8: FINAL SIZE FILTER ON SMA_MASK
//============================================================================
// The AND step with the unclosed Phansalkar mask can fragment ROIs that were
// bridged by morphological closing. This filter removes any resulting fragments
// below CFG_SMA_MIN_SIZE, ensuring the config value refers to the final mask.
if (isOpen("SMA_mask")) {
    selectImage("SMA_mask");
    roiManager("reset");
    run("Analyze Particles...", "size=" + CFG_SMA_MIN_SIZE + "-Infinity show=Nothing add");
    var nSmaFinal = roiManager("count");
    if (nSmaFinal > 0) {
        newImage("SMA_mask_sizefiltered", "8-bit black", getWidth(), getHeight(), 1);
        selectImage("SMA_mask_sizefiltered");
        if (nSmaFinal > 1) { roiManager("Deselect"); roiManager("Combine"); }
        else { roiManager("select", 0); }
        setForegroundColor(255, 255, 255);
        run("Fill", "slice");
        run("Select None");
        // AND with original to preserve internal structure (holes, irregular edges)
        imageCalculator("AND", "SMA_mask_sizefiltered", "SMA_mask");
        safeClose("SMA_mask");
        selectImage("SMA_mask_sizefiltered");
        rename("SMA_mask");
    } else {
        safeClose("SMA_mask");
        newImage("SMA_mask", "8-bit black", getWidth(), getHeight(), 1);
    }
    roiManager("reset");
    logStep("  - SMA final size filter (>=" + CFG_SMA_MIN_SIZE + " px2): " + nSmaFinal + " ROIs retained");
}

//============================================================================
// PHASE 4: COMPARTMENT MASK GENERATION (Steps 15–18)
//============================================================================
// Steps 15-18 produce the four mutually exclusive vascular and tissue
// compartments: SMA-negative basement membrane vessels (capillaries), SMA-
// positive arterioles (with vessel lumens subtracted), large vessel walls
// (ColIV rim around SMA), and non-vascular tissue. The lumen-subtraction
// of arterioles (Step 16) is implemented in the SMA-hole-detection section
// below, while the capillary partition (Step 15) and the large-vessel
// derivation (Step 17) are interleaved further down in a single section
// because both are computed from the same ColIV ROI pool.

//============================================================================
// GRAYSCALE-BASED HOLE DETECTION FOR SMA MASK
//============================================================================
// STEP 16: SMA-positive arterioles = validated SMA mask with lumen subtraction
// Detected lumen holes are subtracted from the SMA mask so that downstream
// area- and intensity-based amyloid load measurements (IntDen / Area) do
// not divide by an area that includes empty lumen pixels.
if (isOpen("SMA_mask") && isOpen("MAX_C1-A1")) {
    logStep("Detecting vessel lumen holes based on grayscale intensity...");
    
    selectImage("SMA_mask");
    run("Duplicate...", "title=SMA_mask_filled");
    run("Fill Holes");
    run("Options...", "iterations=2 count=1 black do=Dilate");
    run("Options...", "iterations=2 count=1 black do=Erode");
    run("Fill Holes");
    
    selectImage("MAX_C1-A1");
    run("Duplicate...", "title=dark_regions");
    
    var holeArea, holeMean, holeMin, holeMax, holeStd;
    getStatistics(holeArea, holeMean, holeMin, holeMax, holeStd);
    var darkThreshold = holeMean * CFG_HOLE_DARK_THRESHOLD;
    
    setThreshold(0, darkThreshold);
    run("Convert to Mask");
    
    imageCalculator("AND create", "dark_regions", "SMA_mask_filled");
    rename("holes_inside_vessels");
    
    // v8.7 FIX: Use ROI-based mask creation instead of show=Masks
    // show=Masks can produce inverted LUT images that corrupt downstream Analyze Particles
    selectImage("holes_inside_vessels");
    roiManager("reset");
    run("Analyze Particles...", "size=" + CFG_MIN_HOLE_SIZE + "-Infinity show=Nothing add");
    var nHoleROIs = roiManager("count");
    newImage("final_holes", "8-bit black", getWidth(), getHeight(), 1);
    if (nHoleROIs > 0) {
        selectImage("final_holes");
        if (nHoleROIs > 1) { roiManager("Deselect"); roiManager("Combine"); }
        else { roiManager("select", 0); }
        setForegroundColor(255, 255, 255);
        run("Fill", "slice");
        run("Select None");
    }
    roiManager("reset");
    
    imageCalculator("Subtract", "SMA_mask", "final_holes");
    
    selectImage("final_holes");
    var finalHoleArea, finalHoleMean;
    getStatistics(finalHoleArea, finalHoleMean);
    var holesArea = finalHoleArea * finalHoleMean / 255;
    if (holesArea > 0) {
        logStep("  - Detected and removed vessel lumen holes (total area: " + d2s(holesArea, 0) + " px2)");
    } else {
        logStep("  - No vessel lumen holes detected");
    }
    
    safeClose("SMA_mask_filled");
    safeClose("dark_regions");
    safeClose("holes_inside_vessels");
    safeClose("final_holes");
}

//============================================================================
// CREATE CAPILLARIES MASK BY EXCLUDING SMA-POSITIVE VESSELS
//============================================================================
// STEP 15: SMA-negative BM vessels = ColIV minus SMA-overlapping ROIs
// STEP 17: Large vessel walls (ColIV ∩ SMA, minus SMA itself) — computed in
//          the same section because both partitions are derived from the
//          same ColIV ROI pool. See submarker further down at the
//          imageCalculator("Subtract", "pre_large_vessels_mask", "SMA_mask")
//          line where Step 17's defining operation begins.
logStep("Creating capillaries mask (excluding SMA-positive vessels)...");

roiManager("reset");
selectImage("collagen4_mask");
run("Analyze Particles...", "size=" + CFG_COLLAGEN_MIN_SIZE + "-Infinity circularity=0.00-1.00 show=Nothing add composite");

var totalROIs = roiManager("count");
if (totalROIs == 0) {
    logStep("WARNING: No collagen4 ROIs found. Creating empty capillaries_mask and large_vessels_mask.");
    newImage("capillaries_mask", "8-bit black", getWidth(), getHeight(), 1);
    newImage("large_vessels_mask", "8-bit black", getWidth(), getHeight(), 1);
} else {
    logStep("  - Total collagen4 ROIs: " + totalROIs);
    
    run("Set Measurements...", "mean redirect=None decimal=3");
    selectImage("SMA_mask");
    roiManager("Deselect");
    run("Clear Results");
    roiManager("Multi-Measure measure_all");
    
    var toDelete = newArray(totalROIs);
    var deleteCount = 0;
    
    for (i = 0; i < totalROIs; i++) {
        meanIntensity = getResult("Mean", i);
        var overlapPercent = (meanIntensity / 255) * 100;
        if (overlapPercent > CFG_SMA_EXCLUSION_THRESH) {
            toDelete[deleteCount] = i;
            deleteCount++;
        }
    }
    
    run("Clear Results");
    
    // --- CREATE LARGE VESSELS MASK from SMA-overlapping collagen4 ROIs ---
    // Large Vessels = collagen4 objects that overlap SMA, with SMA portion subtracted
    if (deleteCount > 0) {
        var smaOverlapIndices = Array.trim(toDelete, deleteCount);
        logStep("  - Creating Large Vessels mask from " + deleteCount + " SMA-overlapping collagen4 ROIs...");
        
        w = getWidth();
        h = getHeight();
        newImage("LV_selection_mask", "8-bit black", w, h, 1);
        if (deleteCount > 1) {
            roiManager("select", smaOverlapIndices);
            roiManager("Combine");
        } else {
            roiManager("select", smaOverlapIndices[0]);
        }
        setForegroundColor(255, 255, 255);
        run("Fill", "slice");
        run("Select None");
        
        // AND with auto_local_collagen4 to preserve holes (same as capillaries approach)
        imageCalculator("AND create", "auto_local_collagen4", "LV_selection_mask");
        rename("pre_large_vessels_mask");
        safeClose("LV_selection_mask");
        
        // STEP 17: Large vessel walls = ColIV rim around arterioles
        //          (ColIV ∩ SMA, minus SMA itself)
        // Subtract SMA_mask to get only the collagen4 rim
        imageCalculator("Subtract", "pre_large_vessels_mask", "SMA_mask");
        rename("pre_large_vessels_sizefilter");
        safeClose("pre_large_vessels_mask");
        
        // Size filter: remove small fragments
        // v8.7 FIX: After show=Masks, normalize the binary mask to guarantee
        // correct pixel values (255=object, 0=background) regardless of LUT state.
        // show=Masks can produce images with inverted LUT that corrupt downstream
        // Analyze Particles calls during amyloid quantification.
        selectImage("pre_large_vessels_sizefilter");
        run("Analyze Particles...", "size=" + CFG_LV_MIN_SIZE + "-Infinity show=Masks");
        if (isOpen("Mask of pre_large_vessels_sizefilter")) {
            selectImage("Mask of pre_large_vessels_sizefilter");
            // Normalize: ensure white objects on black background
            setOption("BlackBackground", true);
            run("Make Binary");
            rename("large_vessels_mask");
        } else {
            newImage("large_vessels_mask", "8-bit black", getWidth(), getHeight(), 1);
        }
        safeClose("pre_large_vessels_sizefilter");
        
        // Check if anything remains after subtraction and size filter
        selectImage("large_vessels_mask");
        var lvStatArea, lvStatMean;
        getStatistics(lvStatArea, lvStatMean);
        var lvTotalArea = lvStatArea * lvStatMean / 255;
        logStep("  - Large Vessels mask area (collagen4 rim - SMA, min " + CFG_LV_MIN_SIZE + " px2): " + d2s(lvTotalArea, 0) + " px");
        
        // Now delete the SMA-overlapping ROIs from ROI Manager for capillaries
        if (deleteCount > 1) {
            roiManager("select", smaOverlapIndices);
        } else {
            roiManager("select", smaOverlapIndices[0]);
        }
        roiManager("delete");
    } else {
        newImage("large_vessels_mask", "8-bit black", getWidth(), getHeight(), 1);
        logStep("  - No SMA-overlapping collagen4 ROIs found. Created empty large_vessels_mask.");
    }
    
    var capillaryCount = roiManager("count");
    logStep("  - SMA exclusion threshold: >" + CFG_SMA_EXCLUSION_THRESH + "%");
    logStep("  - SMA-positive vessels removed: " + deleteCount);
    logStep("  - Capillary ROIs remaining: " + capillaryCount);
    
    if (capillaryCount > 0) {
        w = getWidth();
        h = getHeight();
        newImage("Cap_selection_mask", "8-bit black", w, h, 1);
        roiManager("Deselect");
        roiManager("Combine");
        setForegroundColor(255, 255, 255);
        run("Fill", "slice");
        run("Select None");
        imageCalculator("AND create", "auto_local_collagen4", "Cap_selection_mask");
        rename("capillaries_mask");
        safeClose("Cap_selection_mask");
        logStep("  - Created capillaries_mask with " + capillaryCount + " capillary ROIs (AND with threshold to preserve holes)");
    } else {
        newImage("capillaries_mask", "8-bit black", getWidth(), getHeight(), 1);
        logStep("WARNING: No capillaries passed SMA exclusion. Created empty capillaries_mask.");
    }
}

//============================================================================
// v8.8: FINAL SIZE FILTER ON CAPILLARIES_MASK
//============================================================================
// Same issue as SMA: the AND step with the unclosed Phansalkar mask can fragment
// ROIs below CFG_COLLAGEN_MIN_SIZE. This ensures the config value is enforced.
if (isOpen("capillaries_mask")) {
    selectImage("capillaries_mask");
    roiManager("reset");
    run("Analyze Particles...", "size=" + CFG_COLLAGEN_MIN_SIZE + "-Infinity show=Nothing add");
    var nCapFinal = roiManager("count");
    if (nCapFinal > 0) {
        newImage("cap_mask_sizefiltered", "8-bit black", getWidth(), getHeight(), 1);
        selectImage("cap_mask_sizefiltered");
        if (nCapFinal > 1) { roiManager("Deselect"); roiManager("Combine"); }
        else { roiManager("select", 0); }
        setForegroundColor(255, 255, 255);
        run("Fill", "slice");
        run("Select None");
        // AND with original to preserve internal structure (holes, irregular edges)
        imageCalculator("AND", "cap_mask_sizefiltered", "capillaries_mask");
        safeClose("capillaries_mask");
        selectImage("cap_mask_sizefiltered");
        rename("capillaries_mask");
    } else {
        safeClose("capillaries_mask");
        newImage("capillaries_mask", "8-bit black", getWidth(), getHeight(), 1);
    }
    roiManager("reset");
    logStep("  - Capillaries final size filter (>=" + CFG_COLLAGEN_MIN_SIZE + " px2): " + nCapFinal + " ROIs retained");
}

//============================================================================
// GRAYSCALE-BASED HOLE DETECTION FOR CAPILLARIES MASK
//============================================================================
if (isOpen("capillaries_mask") && isOpen("MAX_C3-A1")) {
    logStep("Detecting capillary lumen holes based on grayscale intensity...");
    
    selectImage("capillaries_mask");
    run("Duplicate...", "title=cap_mask_filled");
    run("Fill Holes");
    run("Options...", "iterations=2 count=1 black do=Dilate");
    run("Options...", "iterations=2 count=1 black do=Erode");
    run("Fill Holes");
    
    selectImage("MAX_C3-A1");
    run("Duplicate...", "title=cap_dark_regions");
    
    var capHoleArea, capHoleMean, capHoleMin, capHoleMax, capHoleStd;
    getStatistics(capHoleArea, capHoleMean, capHoleMin, capHoleMax, capHoleStd);
    var darkThreshold = capHoleMean * CFG_HOLE_DARK_THRESHOLD;
    
    setThreshold(0, darkThreshold);
    run("Convert to Mask");
    
    imageCalculator("AND create", "cap_dark_regions", "cap_mask_filled");
    rename("cap_holes_inside_vessels");
    
    // v8.7 FIX: Use ROI-based mask creation instead of show=Masks
    selectImage("cap_holes_inside_vessels");
    roiManager("reset");
    run("Analyze Particles...", "size=" + CFG_MIN_HOLE_SIZE + "-Infinity show=Nothing add");
    var nCapHoleROIs = roiManager("count");
    newImage("cap_final_holes", "8-bit black", getWidth(), getHeight(), 1);
    if (nCapHoleROIs > 0) {
        selectImage("cap_final_holes");
        if (nCapHoleROIs > 1) { roiManager("Deselect"); roiManager("Combine"); }
        else { roiManager("select", 0); }
        setForegroundColor(255, 255, 255);
        run("Fill", "slice");
        run("Select None");
    }
    roiManager("reset");
    
    imageCalculator("Subtract", "capillaries_mask", "cap_final_holes");
    
    selectImage("cap_final_holes");
    var capFinalHoleArea, capFinalHoleMean;
    getStatistics(capFinalHoleArea, capFinalHoleMean);
    var capHolesArea = capFinalHoleArea * capFinalHoleMean / 255;
    if (capHolesArea > 0) {
        logStep("  - Detected and removed capillary lumen holes (total area: " + d2s(capHolesArea, 0) + " px2)");
    } else {
        logStep("  - No capillary lumen holes detected");
    }
    
    safeClose("cap_mask_filled");
    safeClose("cap_dark_regions");
    safeClose("cap_holes_inside_vessels");
    safeClose("cap_final_holes");
}

//============================================================================
// CLEANUP
//============================================================================
run("Clear Results");
roiManager("reset");
// collagen4_mask kept open for non-vascular tissue subtraction
safeClose("auto_local_SMA");
safeClose("auto_local_SMA_closed");
safeClose("auto_local_collagen4");

//============================================================================
// PHASE 5: AMYLOID QUANTIFICATION (Steps 19–20)
//============================================================================
// Phase 5 produces the binary amyloid mask (Step 19) and overlays it on
// each of the four compartment masks to compute per-compartment amyloid
// load (Step 20). Step 20 occurs four times in the code, once per
// compartment, because each overlay-and-measure pair generates its own
// log window output.

//============================================================================
// CREATE AMYLOID MASK (GLOBAL THRESHOLD)
//============================================================================
// STEP 19: Fixed global threshold on preprocessed C2 → binary amyloid mask
logStep("Creating amyloid mask (threshold >= " + CFG_AMYLOID_THRESHOLD + ")...");
selectImage("MAX_C2-A1");
run("Duplicate...", "title=amyloid_mask");
setThreshold(CFG_AMYLOID_THRESHOLD, 255);
setOption("BlackBackground", true);
run("Convert to Mask");

// Count amyloid objects for QC
roiManager("reset");
run("Set Measurements...", "area mean redirect=None decimal=3");
selectImage("amyloid_mask");
run("Analyze Particles...", "size=0-Infinity show=Nothing display");

var amyloid_count = nResults;
var amyloid_mean_size = 0;
var amyloid_max_size = 0;
var amyloid_stddev_size = 0;
var amyloid_mean_intensity = 0;
var amyloid_max_intensity = 0;

if (amyloid_count > 0) {
    run("Summarize");
    var totalRows = nResults;
    if (totalRows >= 4) {
        amyloid_mean_size = getResult("Area", totalRows - 4);
        amyloid_stddev_size = getResult("Area", totalRows - 3);
        amyloid_max_size = getResult("Area", totalRows - 1);
    } else if (totalRows == 1) {
        amyloid_mean_size = getResult("Area", 0);
        amyloid_stddev_size = 0;
        amyloid_max_size = getResult("Area", 0);
    }
    run("Clear Results");

    selectImage("amyloid_mask");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("MAX_C2-A1");
        run("Restore Selection");
        var amyQcArea, amyQcMin, amyQcMax, amyQcStd;
        getStatistics(amyQcArea, amyloid_mean_intensity, amyQcMin, amyloid_max_intensity, amyQcStd);
        run("Select None");
    }
} else {
    run("Clear Results");
}

roiManager("reset");

print("--- AMYLOID ---");
print("  Object Count: " + amyloid_count);
print("  Mean Size: " + d2s(amyloid_mean_size, 2) + " px2");
print("  Max Size: " + d2s(amyloid_max_size, 2) + " px2");
print("  StdDev Size: " + d2s(amyloid_stddev_size, 2) + " px2");
print("  Mean Pixel Intensity (C2): " + d2s(amyloid_mean_intensity, 2));
print("  Max Pixel Intensity (C2): " + d2s(amyloid_max_intensity, 2));

//============================================================================
// COMPREHENSIVE QC LOGGING
//============================================================================
logStep("Generating QC statistics...");

if (!isOpen("QC_Log")) {
    run("New... ", "name=QC_Log type=Table");
    print("[QC_Log]", "\\Headings:Image_Name\t" +
          "Cap_Count\tCap_MeanArea\tCap_StdDevArea\tCap_MaxArea\tCap_MeanIntensity\tCap_StdDevIntensity\t" +
          "SMA_Count\tSMA_MeanArea\tSMA_StdDevArea\tSMA_MaxArea\tSMA_MeanIntensity\tSMA_StdDevIntensity\t" +
          "LV_Count\tLV_MeanArea\tLV_StdDevArea\tLV_MaxArea\tLV_MeanIntensity\tLV_StdDevIntensity\t" +
          "Amy_Count\tAmy_MeanSize\tAmy_MaxSize\tAmy_StdDevSize\tAmy_MeanIntensity\tAmy_MaxIntensity\t" +
          "EmptySpace_Area\tEmptySpace_MeanIntensity\t" +
          "MeanInt_C1\tMeanInt_C2\tMeanInt_C3\t" +
          "SMA_Factor_L5\tSMA_Factor_L4\tSMA_Factor_L3\tSMA_Factor_L2\tSMA_Factor_L1\t" +
          "SMA_Kept_L5\tSMA_Kept_L4\tSMA_Kept_L3\tSMA_Kept_L2\tSMA_Kept_L1\tSMA_TotalRemoved");
}

print("================================================================================");
print("QC SUMMARY: " + originalName);
print("================================================================================");

//==================================================================================
// QC FOR CAPILLARIES
//==================================================================================
var cap_ObjectCount = 0;
var cap_MeanArea = 0;
var cap_StdDevArea = 0;
var cap_MaxArea = 0;
var cap_MeanIntensity = 0;
var cap_StdDevIntensity = 0;

roiManager("reset");
selectImage("capillaries_mask");
run("Select None");  // v9.0: clear any stale selection inherited from prior image
run("Set Measurements...", "area redirect=None decimal=3");
run("Analyze Particles...", "size=0-Infinity show=Nothing display");

cap_ObjectCount = nResults;

if (cap_ObjectCount > 0) {
    run("Summarize");
    var totalRows = nResults;
    if (totalRows >= 4) {
        var meanRow = totalRows - 4;
        var sdRow = totalRows - 3;
        var maxRow = totalRows - 1;
        cap_MeanArea = getResult("Area", meanRow);
        cap_StdDevArea = getResult("Area", sdRow);
        cap_MaxArea = getResult("Area", maxRow);
    } else if (totalRows == 1) {
        cap_MeanArea = getResult("Area", 0);
        cap_StdDevArea = 0;
        cap_MaxArea = getResult("Area", 0);
    }
}
run("Clear Results");

if (cap_ObjectCount > 0) {
    roiManager("reset");
    selectImage("capillaries_mask");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("MAX_C3-A1");
        run("Restore Selection");
        var capQcArea, capQcMin, capQcMax;
        getStatistics(capQcArea, cap_MeanIntensity, capQcMin, capQcMax, cap_StdDevIntensity);
        run("Select None");
    }
}

print("--- CAPILLARIES ---");
print("  Object Count: " + cap_ObjectCount);
print("  Mean Area: " + d2s(cap_MeanArea, 2) + " px2");
print("  StdDev Area: " + d2s(cap_StdDevArea, 2) + " px2");
print("  Max Area: " + d2s(cap_MaxArea, 2) + " px2");
print("  Mean Intensity (C3): " + d2s(cap_MeanIntensity, 2));
print("  StdDev Intensity (C3): " + d2s(cap_StdDevIntensity, 2));

//==================================================================================
// QC FOR SMA
//==================================================================================
// v9.0: Count and measure from the post-v8.8-size-filter SMA_mask, the same
// mask the amyloid analysis operates on. Previously this section read from
// sma_areas_for_qc[] which captured the count BEFORE the v8.8 size filter,
// causing the QC log to disagree with the SMA_Amyloid window. The upstream
// sma_areas_for_qc[] / sma_roi_count_for_qc variables are still populated
// during level-classification but are no longer used here. This matches how
// the capillary QC counter works.
var sma_ObjectCount = 0;
var sma_MeanArea = 0;
var sma_StdDevArea = 0;
var sma_MaxArea = 0;
var sma_MeanIntensity = 0;
var sma_StdDevIntensity = 0;

roiManager("reset");
selectImage("SMA_mask");
run("Select None");  // v9.0: clear any stale selection inherited from prior image
run("Set Measurements...", "area redirect=None decimal=3");
run("Analyze Particles...", "size=0-Infinity show=Nothing display");

sma_ObjectCount = nResults;

if (sma_ObjectCount > 0) {
    run("Summarize");
    var smaTotalRows = nResults;
    if (smaTotalRows >= 4) {
        var smaMeanRow = smaTotalRows - 4;
        var smaSdRow = smaTotalRows - 3;
        var smaMaxRow = smaTotalRows - 1;
        sma_MeanArea = getResult("Area", smaMeanRow);
        sma_StdDevArea = getResult("Area", smaSdRow);
        sma_MaxArea = getResult("Area", smaMaxRow);
    } else if (smaTotalRows == 1) {
        sma_MeanArea = getResult("Area", 0);
        sma_StdDevArea = 0;
        sma_MaxArea = getResult("Area", 0);
    }
}
run("Clear Results");

if (sma_ObjectCount > 0 && isOpen("SMA_mask")) {
    roiManager("reset");
    selectImage("SMA_mask");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("MAX_C1-A1");
        run("Restore Selection");
        var smaQcArea, smaQcMin, smaQcMax;
        getStatistics(smaQcArea, sma_MeanIntensity, smaQcMin, smaQcMax, sma_StdDevIntensity);
        run("Select None");
    }
}

roiManager("reset");

print("--- SMA POSITIVE ---");
print("  Object Count: " + sma_ObjectCount);
print("  Mean Area: " + d2s(sma_MeanArea, 2) + " px2");
print("  StdDev Area: " + d2s(sma_StdDevArea, 2) + " px2");
print("  Max Area: " + d2s(sma_MaxArea, 2) + " px2");
print("  Mean Intensity (C1): " + d2s(sma_MeanIntensity, 2));
print("  StdDev Intensity (C1): " + d2s(sma_StdDevIntensity, 2));
print("================================================================================");

//==================================================================================
// QC FOR LARGE VESSELS
//==================================================================================
var lv_ObjectCount = 0;
var lv_MeanArea = 0;
var lv_StdDevArea = 0;
var lv_MaxArea = 0;
var lv_MeanIntensity = 0;
var lv_StdDevIntensity = 0;

roiManager("reset");
selectImage("large_vessels_mask");
run("Select None");  // v9.0 BUGFIX: clear stale selection inherited from prior SMA QC code that operated on MAX_C1-A1. Without this, AP only counts within the inherited selection bounds, severely undercounting LVs.
run("Set Measurements...", "area redirect=None decimal=3");
run("Analyze Particles...", "size=0-Infinity show=Nothing display");

lv_ObjectCount = nResults;

if (lv_ObjectCount > 0) {
    run("Summarize");
    var totalRows = nResults;
    if (totalRows >= 4) {
        var meanRow = totalRows - 4;
        var sdRow = totalRows - 3;
        var maxRow = totalRows - 1;
        lv_MeanArea = getResult("Area", meanRow);
        lv_StdDevArea = getResult("Area", sdRow);
        lv_MaxArea = getResult("Area", maxRow);
    } else if (totalRows == 1) {
        lv_MeanArea = getResult("Area", 0);
        lv_StdDevArea = 0;
        lv_MaxArea = getResult("Area", 0);
    }
}
run("Clear Results");

if (lv_ObjectCount > 0 && isOpen("large_vessels_mask")) {
    selectImage("large_vessels_mask");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("MAX_C3-A1");
        run("Restore Selection");
        var lvQcArea, lvQcMin, lvQcMax;
        getStatistics(lvQcArea, lv_MeanIntensity, lvQcMin, lvQcMax, lv_StdDevIntensity);
        run("Select None");
    }
}

roiManager("reset");

print("--- LARGE VESSELS ---");
print("  Object Count: " + lv_ObjectCount);
print("  Mean Area: " + d2s(lv_MeanArea, 2) + " px2");
print("  StdDev Area: " + d2s(lv_StdDevArea, 2) + " px2");
print("  Max Area: " + d2s(lv_MaxArea, 2) + " px2");
print("  Mean Intensity (C3): " + d2s(lv_MeanIntensity, 2));
print("  StdDev Intensity (C3): " + d2s(lv_StdDevIntensity, 2));
print("================================================================================");

// NOTE: QC_Log data row is printed after non-vascular tissue analysis (needs emptySpace values)

//============================================================================
// AMYLOID OVERLAP ANALYSIS - CAPILLARIES (Log Window 1)
//============================================================================
// STEP 20 (Capillaries): Overlay amyloid mask on the capillary compartment
// → per-ROI measurements (size, total amyloid area, mean area, count,
// density, integrated density). Step 20 recurs in three more sections
// below for SMA/arterioles, large vessels, and non-vascular tissue.
logStep("Amyloid overlap analysis: Capillaries...");

if (!isOpen("Capillaries_Amyloid")) {
    run("New... ", "name=Capillaries_Amyloid type=Table");
    print("[Capillaries_Amyloid]", "\\Headings:Image_ROI\tCapillary_size\tTotal_amyloid_area\tAmyloid_mean_area\tAmyloid_count\tAmyloid_density\tAmyloid_IntDen");
}

roiManager("reset");
selectImage("capillaries_mask");
run("Select None");  // v9.0: clear any stale selection inherited from prior image
run("Analyze Particles...", "size=0-Infinity show=Nothing add composite");
var nCapROIs = roiManager("count");
logStep("  - Capillary ROIs for amyloid analysis: " + nCapROIs);

for (var ci = 0; ci < nCapROIs; ci++) {
    // 1. Capillary wall area (composite ROI excludes lumen holes)
    selectImage("capillaries_mask");
    roiManager("select", ci);
    var capStatArea, capStatMean, capStatMin, capStatMax;
    getStatistics(capStatArea, capStatMean, capStatMin, capStatMax);
    var capSize = capStatArea * capStatMean / 255;

    // 2. Total amyloid area within capillary wall
    selectImage("amyloid_mask");
    roiManager("select", ci);
    var amyStatArea, amyStatMean, amyStatMin, amyStatMax;
    getStatistics(amyStatArea, amyStatMean, amyStatMin, amyStatMax);
    var totalAmyloidArea = amyStatArea * amyStatMean / 255;

    // 3. Particle count + intensity in overlap region
    var amyloidCount = 0;
    var amyloidMeanArea = 0;
    var amyloidDensity = 0;
    var amyloidIntDen = 0;

    if (totalAmyloidArea > 0) {
        // FIX: Clear selection BEFORE duplicate to get full-size copy
        selectImage("amyloid_mask");
        run("Select None");
        run("Duplicate...", "title=temp_amy_cap");
        roiManager("select", ci);
        run("Clear Outside");
        run("Select None");

        run("Clear Results");
        run("Set Measurements...", "area redirect=None decimal=3");
        run("Analyze Particles...", "size=0-Infinity show=Nothing display");
        amyloidCount = nResults;
        if (amyloidCount > 0) {
            amyloidMeanArea = totalAmyloidArea / amyloidCount;
        }
        run("Clear Results");

        // Intensity on MAX_C2-A1 — overlap pixels only
        selectImage("temp_amy_cap");
        run("Create Selection");
        if (selectionType() != -1) {
            selectImage("MAX_C2-A1");
            run("Restore Selection");
            var olArea, olMean, olMin, olMax, olStd;
            getStatistics(olArea, olMean, olMin, olMax, olStd);
            amyloidDensity = olMean;
            amyloidIntDen = olMean * olArea;
            run("Select None");
        }
        safeClose("temp_amy_cap");
    }

    var idx = ci + 1;
    var padded = "" + idx;
    if (idx < 10) padded = "00" + idx;
    else if (idx < 100) padded = "0" + idx;

    print("[Capillaries_Amyloid]", originalName + ", Capillary_" + padded + "\t" +
          d2s(capSize, 2) + "\t" + d2s(totalAmyloidArea, 2) + "\t" +
          d2s(amyloidMeanArea, 2) + "\t" + amyloidCount + "\t" +
          d2s(amyloidDensity, 2) + "\t" + d2s(amyloidIntDen, 2));
}

roiManager("reset");
logStep("  - Capillary amyloid analysis complete: " + nCapROIs + " ROIs processed");

//============================================================================
// AMYLOID OVERLAP ANALYSIS - SMA / ARTERIOLES (Log Window 2)
//============================================================================
// STEP 20 (SMA / arterioles): per-ROI overlay and measurements.
logStep("Amyloid overlap analysis: SMA objects...");

if (!isOpen("SMA_Amyloid")) {
    run("New... ", "name=SMA_Amyloid type=Table");
    print("[SMA_Amyloid]", "\\Headings:Image_ROI\tSMA_size\tTotal_amyloid_area\tAmyloid_mean_area\tAmyloid_count\tAmyloid_density\tAmyloid_IntDen");
}

roiManager("reset");
selectImage("SMA_mask");
run("Select None");  // v9.0: clear any stale selection inherited from prior image
run("Analyze Particles...", "size=0-Infinity show=Nothing add composite");
var nSmaROIs = roiManager("count");
logStep("  - SMA ROIs for amyloid analysis: " + nSmaROIs);

for (var si = 0; si < nSmaROIs; si++) {
    // 1. SMA object area (composite ROI excludes lumen holes)
    selectImage("SMA_mask");
    roiManager("select", si);
    var smaStatArea, smaStatMean, smaStatMin, smaStatMax;
    getStatistics(smaStatArea, smaStatMean, smaStatMin, smaStatMax);
    var smaSize = smaStatArea * smaStatMean / 255;

    // 2. Total amyloid area within SMA object
    selectImage("amyloid_mask");
    roiManager("select", si);
    var amySmaArea, amySmaMean, amySmaMin, amySmaMax;
    getStatistics(amySmaArea, amySmaMean, amySmaMin, amySmaMax);
    var totalAmyloidArea = amySmaArea * amySmaMean / 255;

    // 3. Particle count + intensity in overlap region
    var amyloidCount = 0;
    var amyloidMeanArea = 0;
    var amyloidDensity = 0;
    var amyloidIntDen = 0;

    if (totalAmyloidArea > 0) {
        // FIX: Clear selection BEFORE duplicate to get full-size copy
        selectImage("amyloid_mask");
        run("Select None");
        run("Duplicate...", "title=temp_amy_sma");
        roiManager("select", si);
        run("Clear Outside");
        run("Select None");

        run("Clear Results");
        run("Set Measurements...", "area redirect=None decimal=3");
        run("Analyze Particles...", "size=0-Infinity show=Nothing display");
        amyloidCount = nResults;
        if (amyloidCount > 0) {
            amyloidMeanArea = totalAmyloidArea / amyloidCount;
        }
        run("Clear Results");

        selectImage("temp_amy_sma");
        run("Create Selection");
        if (selectionType() != -1) {
            selectImage("MAX_C2-A1");
            run("Restore Selection");
            var olArea, olMean, olMin, olMax, olStd;
            getStatistics(olArea, olMean, olMin, olMax, olStd);
            amyloidDensity = olMean;
            amyloidIntDen = olMean * olArea;
            run("Select None");
        }
        safeClose("temp_amy_sma");
    }

    var idx = si + 1;
    var padded = "" + idx;
    if (idx < 10) padded = "00" + idx;
    else if (idx < 100) padded = "0" + idx;

    print("[SMA_Amyloid]", originalName + ", SMA_" + padded + "\t" +
          d2s(smaSize, 2) + "\t" + d2s(totalAmyloidArea, 2) + "\t" +
          d2s(amyloidMeanArea, 2) + "\t" + amyloidCount + "\t" +
          d2s(amyloidDensity, 2) + "\t" + d2s(amyloidIntDen, 2));
}

roiManager("reset");
logStep("  - SMA amyloid analysis complete: " + nSmaROIs + " ROIs processed");

//============================================================================
// AMYLOID OVERLAP ANALYSIS - LARGE VESSELS (Log Window 3)
//============================================================================
// STEP 20 (Large vessels): per-ROI overlay and measurements.
logStep("Amyloid overlap analysis: Large Vessels...");

if (!isOpen("Large_Vessels_Amyloid")) {
    run("New... ", "name=Large_Vessels_Amyloid type=Table");
    print("[Large_Vessels_Amyloid]", "\\Headings:Image_ROI\tLV_size\tTotal_amyloid_area\tAmyloid_mean_area\tAmyloid_count\tAmyloid_density\tAmyloid_IntDen");
}

roiManager("reset");
selectImage("large_vessels_mask");
run("Select None");  // v9.0: clear any stale selection inherited from prior image
run("Analyze Particles...", "size=0-Infinity show=Nothing add composite");
var nLvROIs = roiManager("count");
logStep("  - Large Vessel ROIs for amyloid analysis: " + nLvROIs);

for (var li = 0; li < nLvROIs; li++) {
    // 1. Large vessel area (composite ROI)
    selectImage("large_vessels_mask");
    roiManager("select", li);
    var lvStatArea, lvStatMean, lvStatMin, lvStatMax;
    getStatistics(lvStatArea, lvStatMean, lvStatMin, lvStatMax);
    var lvSize = lvStatArea * lvStatMean / 255;

    // 2. Total amyloid area within large vessel
    selectImage("amyloid_mask");
    roiManager("select", li);
    var amyLvArea, amyLvMean, amyLvMin, amyLvMax;
    getStatistics(amyLvArea, amyLvMean, amyLvMin, amyLvMax);
    var totalAmyloidArea = amyLvArea * amyLvMean / 255;

    // 3. Particle count + intensity in overlap region
    var amyloidCount = 0;
    var amyloidMeanArea = 0;
    var amyloidDensity = 0;
    var amyloidIntDen = 0;

    if (totalAmyloidArea > 0) {
        selectImage("amyloid_mask");
        run("Select None");
        run("Duplicate...", "title=temp_amy_lv");
        roiManager("select", li);
        run("Clear Outside");
        run("Select None");

        run("Clear Results");
        run("Set Measurements...", "area redirect=None decimal=3");
        run("Analyze Particles...", "size=0-Infinity show=Nothing display");
        amyloidCount = nResults;
        if (amyloidCount > 0) {
            amyloidMeanArea = totalAmyloidArea / amyloidCount;
        }
        run("Clear Results");

        selectImage("temp_amy_lv");
        run("Create Selection");
        if (selectionType() != -1) {
            selectImage("MAX_C2-A1");
            run("Restore Selection");
            var olArea, olMean, olMin, olMax, olStd;
            getStatistics(olArea, olMean, olMin, olMax, olStd);
            amyloidDensity = olMean;
            amyloidIntDen = olMean * olArea;
            run("Select None");
        }
        safeClose("temp_amy_lv");
    }

    var idx = li + 1;
    var padded = "" + idx;
    if (idx < 10) padded = "00" + idx;
    else if (idx < 100) padded = "0" + idx;

    print("[Large_Vessels_Amyloid]", originalName + ", LV_" + padded + "\t" +
          d2s(lvSize, 2) + "\t" + d2s(totalAmyloidArea, 2) + "\t" +
          d2s(amyloidMeanArea, 2) + "\t" + amyloidCount + "\t" +
          d2s(amyloidDensity, 2) + "\t" + d2s(amyloidIntDen, 2));
}

roiManager("reset");
logStep("  - Large Vessels amyloid analysis complete: " + nLvROIs + " ROIs processed");

//============================================================================
// NON-VASCULAR TISSUE ANALYSIS (Log Window 4)
//============================================================================
// STEP 18: Non-vascular tissue = image minus all vascular masks minus
//          empty space. Computed by subtracting the four other
//          compartment masks from a full-image white mask.
// STEP 20 (Non-vascular tissue): the amyloid overlay and measurements
//          for this compartment occur further down in this section,
//          after the non-vascular mask itself is built.
logStep("Creating non-vascular tissue mask...");
// empty_space_mask was created earlier (adaptive per-channel or static combined strategy)

// Create non-vascular tissue mask = full image - collagen4 - SMA - large vessels - empty space
logStep("Creating non-vascular tissue mask (full image - collagen4 - SMA - large vessels - empty space)...");
newImage("nonvascular_mask", "8-bit white", imgWidth, imgHeight, 1);
imageCalculator("Subtract", "nonvascular_mask", "collagen4_mask");
imageCalculator("Subtract", "nonvascular_mask", "SMA_mask");
imageCalculator("Subtract", "nonvascular_mask", "large_vessels_mask");
imageCalculator("Subtract", "nonvascular_mask", "empty_space_mask");
safeClose("collagen4_mask");

// Log non-vascular tissue stats
selectImage("nonvascular_mask");
var parStatArea, parStatMean;
getStatistics(parStatArea, parStatMean);
var nonvascularSize = parStatArea * parStatMean / 255;
logStep("  - Non-vascular tissue area (image - collagen4 - SMA - large vessels - empty space): " + d2s(nonvascularSize, 0) + " px");

// Step 5: Find amyloid within non-vascular tissue
logStep("Amyloid overlap analysis: Non-vascular tissue...");

if (!isOpen("NonVascular_Amyloid")) {
    run("New... ", "name=NonVascular_Amyloid type=Table");
    print("[NonVascular_Amyloid]", "\\Headings:Image_name\tNonVascular_size\tTotal_Amyloid_Area\tAmyloid_mean_area\tAmyloid_count\tAmyloid_density\tAmyloid_IntDen");
}

var nv_amyloidCount = 0;
var nv_amyloidMeanArea = 0;
var nv_amyloidDensity = 0;
var nv_amyloidIntDen = 0;
var nv_totalAmyloidArea = 0;

// AND amyloid with non-vascular tissue
selectImage("amyloid_mask");
run("Select None");
selectImage("nonvascular_mask");
run("Select None");
imageCalculator("AND create", "amyloid_mask", "nonvascular_mask");
rename("amy_in_nonvascular");

var apStatArea, apStatMean;
getStatistics(apStatArea, apStatMean);
nv_totalAmyloidArea = apStatArea * apStatMean / 255;
logStep("  - Amyloid pixels in non-vascular tissue: " + d2s(nv_totalAmyloidArea, 0));

if (nv_totalAmyloidArea > 0) {
    // Count particles
    selectImage("amy_in_nonvascular");
    run("Clear Results");
    run("Set Measurements...", "area redirect=None decimal=3");
    run("Analyze Particles...", "size=0-Infinity show=Nothing display");
    nv_amyloidCount = nResults;
    if (nv_amyloidCount > 0) {
        nv_amyloidMeanArea = nv_totalAmyloidArea / nv_amyloidCount;
    }
    run("Clear Results");

    // Intensity measurement on overlap pixels only
    selectImage("amy_in_nonvascular");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("MAX_C2-A1");
        run("Restore Selection");
        var olArea, olMean, olMin, olMax, olStd;
        getStatistics(olArea, olMean, olMin, olMax, olStd);
        nv_amyloidDensity = olMean;
        nv_amyloidIntDen = olMean * olArea;
        run("Select None");
    }
}

safeClose("amy_in_nonvascular");
safeClose("nonvascular_mask");

logStep("  - Non-vascular amyloid: count=" + nv_amyloidCount + " totalArea=" + d2s(nv_totalAmyloidArea, 2) + " density=" + d2s(nv_amyloidDensity, 2));

print("[NonVascular_Amyloid]", originalName + "\t" +
      d2s(nonvascularSize, 2) + "\t" + d2s(nv_totalAmyloidArea, 2) + "\t" +
      d2s(nv_amyloidMeanArea, 2) + "\t" + nv_amyloidCount + "\t" +
      d2s(nv_amyloidDensity, 2) + "\t" + d2s(nv_amyloidIntDen, 2));

// QC_Log data row (after all values including emptySpace are computed)
print("[QC_Log]", originalName + "\t" + 
      cap_ObjectCount + "\t" + d2s(cap_MeanArea, 2) + "\t" + d2s(cap_StdDevArea, 2) + "\t" + 
      d2s(cap_MaxArea, 2) + "\t" + d2s(cap_MeanIntensity, 2) + "\t" + d2s(cap_StdDevIntensity, 2) + "\t" +
      sma_ObjectCount + "\t" + d2s(sma_MeanArea, 2) + "\t" + d2s(sma_StdDevArea, 2) + "\t" + 
      d2s(sma_MaxArea, 2) + "\t" + d2s(sma_MeanIntensity, 2) + "\t" + d2s(sma_StdDevIntensity, 2) + "\t" +
      lv_ObjectCount + "\t" + d2s(lv_MeanArea, 2) + "\t" + d2s(lv_StdDevArea, 2) + "\t" + 
      d2s(lv_MaxArea, 2) + "\t" + d2s(lv_MeanIntensity, 2) + "\t" + d2s(lv_StdDevIntensity, 2) + "\t" +
      amyloid_count + "\t" + d2s(amyloid_mean_size, 2) + "\t" + d2s(amyloid_max_size, 2) + "\t" +
      d2s(amyloid_stddev_size, 2) + "\t" + d2s(amyloid_mean_intensity, 2) + "\t" + d2s(amyloid_max_intensity, 2) + "\t" +
      d2s(emptySpaceArea, 2) + "\t" + d2s(emptySpaceMeanIntensity, 2) + "\t" +
      d2s(meanIntensity_C1, 2) + "\t" + d2s(meanIntensity_C2, 2) + "\t" + d2s(meanIntensity_C3, 2) + "\t" +
      d2s(sma_factor_L5, 4) + "\t" + d2s(sma_factor_L4, 4) + "\t" + d2s(sma_factor_L3, 4) + "\t" + d2s(sma_factor_L2, 4) + "\t" + d2s(sma_factor_L1, 4) + "\t" +
      sma_kept_L5 + "\t" + sma_kept_L4 + "\t" + sma_kept_L3 + "\t" + sma_kept_L2 + "\t" + sma_kept_L1 + "\t" + sma_total_removed);

//==================================================================================
// v8.9: VALIDATION EXPORT — binary masks for spatial matching
//==================================================================================
// Saves the four compartment masks and the amyloid mask as binary TIFFs.
// All spatial analysis (connected component labeling, ROI-level amyloid load,
// matching to manual annotations) is performed downstream in R.
// Masks are saved while still in batch mode — no batch toggle needed.
if (CFG_VALIDATION_MODE) {
    logStep("Exporting validation masks...");
    
    File.makeDirectory(CFG_VALIDATION_EXPORT_DIR);
    
    // Strip .tif extension for base filename
    var valBaseName = replace(originalName, ".tif", "");
    
    // --- Binary capillaries mask ---
    selectImage("capillaries_mask");
    run("Select None");
    run("Duplicate...", "title=val_cap_export");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_capillaries_mask.tif");
    close();
    logStep("  - Capillaries binary mask exported");
    
    // --- Binary SMA mask ---
    selectImage("SMA_mask");
    run("Select None");
    run("Duplicate...", "title=val_sma_export");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_sma_mask.tif");
    close();
    logStep("  - SMA binary mask exported");
    
    // --- Binary large vessels mask ---
    selectImage("large_vessels_mask");
    run("Select None");
    run("Duplicate...", "title=val_lv_export");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_largevessels_mask.tif");
    close();
    logStep("  - Large vessels binary mask exported");
    
    // --- Binary amyloid mask ---
    selectImage("amyloid_mask");
    run("Select None");
    run("Duplicate...", "title=val_amy_export");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_amyloid_mask.tif");
    close();
    logStep("  - Amyloid binary mask exported");
    
    // --- Raw MAX projections (no median filter, no background subtraction) ---
    // These are what the observer saw. Used for intensity analysis in R.
    selectImage("RAW_MAX_C1");
    run("Select None");
    run("Duplicate...", "title=val_raw_c1");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_RAW_MAX_C1_SMA.tif");
    close();
    logStep("  - RAW MAX C1 (SMA) exported");
    
    selectImage("RAW_MAX_C2");
    run("Select None");
    run("Duplicate...", "title=val_raw_c2");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_RAW_MAX_C2_Amyloid.tif");
    close();
    logStep("  - RAW MAX C2 (Amyloid) exported");
    
    selectImage("RAW_MAX_C3");
    run("Select None");
    run("Duplicate...", "title=val_raw_c3");
    run("Remove Overlay");
    saveAs("Tiff", CFG_VALIDATION_EXPORT_DIR + valBaseName + "_RAW_MAX_C3_ColIV.tif");
    close();
    logStep("  - RAW MAX C3 (ColIV) exported");
    
    logStep("Validation export complete: " + CFG_VALIDATION_EXPORT_DIR);
}

//==================================================================================
// v8.9: DIAGNOSTIC IMAGE — flat RGB with burned-in mask outlines + calibration
//==================================================================================
// Saves a single flat RGB TIFF: SMA=green, Amyloid=red, ColIV=blue, with
// compartment mask outlines drawn on top. Spatial calibration is set so that
// observer .roi files can be loaded directly on top in Fiji.
// v9.0: wrapped in CFG_VALIDATION_MODE — DIAG is only useful for validation
// workflows (overlaying external .roi annotations), not routine analysis.
if (CFG_VALIDATION_MODE) {
    logStep("Generating diagnostic image...");
    
    // Duplicate MAX projections (originals needed for QC visualization below)
    selectImage("MAX_C1-A1");
    run("Duplicate...", "title=diag_SMA");
    selectImage("MAX_C2-A1");
    run("Duplicate...", "title=diag_Amyloid");
    selectImage("MAX_C3-A1");
    run("Duplicate...", "title=diag_ColIV");
    
    // Merge: c1=red=Amyloid, c2=green=SMA, c3=blue=ColIV
    run("Merge Channels...", "c1=diag_Amyloid c2=diag_SMA c3=diag_ColIV create");
    rename("Diagnostic_Composite");
    
    // Flatten to RGB so .roi overlays are channel-independent
    run("RGB Color");
    rename("DIAG_flat");
    
    // --- Draw mask outlines ---
    
    // Capillaries = blue (0, 100, 255)
    selectImage("capillaries_mask");
    run("Select None");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("DIAG_flat");
        run("Restore Selection");
        setForegroundColor(0, 100, 255);
        run("Draw", "slice");
        run("Select None");
    }
    
    // SMA = green (0, 255, 0)
    selectImage("SMA_mask");
    run("Select None");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("DIAG_flat");
        run("Restore Selection");
        setForegroundColor(0, 255, 0);
        run("Draw", "slice");
        run("Select None");
    }
    
    // Large vessels = cyan (0, 220, 220)
    selectImage("large_vessels_mask");
    run("Select None");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("DIAG_flat");
        run("Restore Selection");
        setForegroundColor(0, 220, 220);
        run("Draw", "slice");
        run("Select None");
    }
    
    // Amyloid = red (255, 50, 50)
    selectImage("amyloid_mask");
    run("Select None");
    run("Create Selection");
    if (selectionType() != -1) {
        selectImage("DIAG_flat");
        run("Restore Selection");
        setForegroundColor(255, 50, 50);
        run("Draw", "slice");
        run("Select None");
    }
    
    // Set spatial calibration so observer .roi coordinates match
    selectImage("DIAG_flat");
    run("Set Scale...", "distance=1 known=" + d2s(0.568, 4) + " unit=um");
    
    // Save
    var diagBaseName = replace(originalName, ".tif", "");
    saveAs("Tiff", CFG_QC_OUTPUT_DIR + "DIAG_" + diagBaseName + ".tif");
    logStep("  - Diagnostic image saved: DIAG_" + diagBaseName + ".tif");
    logStep("  - Background: R=Amyloid, G=SMA, B=ColIV");
    logStep("  - Outlines: Blue=Capillaries, Green=SMA, Cyan=LargeVessels, Red=Amyloid");
    logStep("  - Calibration: 0.568 um/pixel");
    close();
}

//==================================================================================
// QC VISUALIZATION — using Draw approach (avoids ROI Manager color issues)
//==================================================================================
logStep("Generating QC visualization...");

// --- CAPILLARIES PANEL: Blue capillary outlines + Red amyloid outlines ---
selectImage("MAX_C3-A1");
run("Duplicate...", "title=QC_Capillaries_flat");
run("RGB Color");

// Draw blue capillary outlines
selectImage("capillaries_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_Capillaries_flat");
    run("Restore Selection");
    setForegroundColor(0, 100, 255);
    run("Draw", "slice");
    run("Select None");
}

// Draw red amyloid outlines
selectImage("amyloid_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_Capillaries_flat");
    run("Restore Selection");
    setForegroundColor(255, 0, 0);
    run("Draw", "slice");
    run("Select None");
}

// Draw purple empty space outlines
selectImage("empty_space_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_Capillaries_flat");
    run("Restore Selection");
    setForegroundColor(180, 0, 255);
    run("Draw", "slice");
    run("Select None");
}

// Draw teal large vessel outlines
selectImage("large_vessels_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_Capillaries_flat");
    run("Restore Selection");
    setForegroundColor(0, 192, 192);
    run("Draw", "slice");
    run("Select None");
}

// Text labels
selectImage("QC_Capillaries_flat");
setFont("SansSerif", 24, "bold");
setColor(0, 100, 255);
drawString("Capillaries", 10, 30);
setColor(0, 192, 192);
drawString("Large Vessels", 10, 58);
setColor(255, 0, 0);
drawString("Amyloid", 10, 86);
setColor(180, 0, 255);
drawString("Empty space", 10, 114);

// --- SMA PANEL: Green SMA outlines + Red amyloid outlines ---
selectImage("MAX_C1-A1");
run("Duplicate...", "title=QC_SMA_flat");
run("RGB Color");

// Draw green SMA outlines
selectImage("SMA_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_SMA_flat");
    run("Restore Selection");
    setForegroundColor(0, 255, 0);
    run("Draw", "slice");
    run("Select None");
}

// Draw red amyloid outlines
selectImage("amyloid_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_SMA_flat");
    run("Restore Selection");
    setForegroundColor(255, 0, 0);
    run("Draw", "slice");
    run("Select None");
}

// Draw purple empty space outlines
selectImage("empty_space_mask");
run("Select None");
run("Create Selection");
if (selectionType() != -1) {
    selectImage("QC_SMA_flat");
    run("Restore Selection");
    setForegroundColor(180, 0, 255);
    run("Draw", "slice");
    run("Select None");
}

// Text labels
selectImage("QC_SMA_flat");
setFont("SansSerif", 24, "bold");
setColor(0, 255, 0);
drawString("SMA", 10, 30);
setColor(255, 0, 0);
drawString("Amyloid", 10, 58);
setColor(180, 0, 255);
drawString("Empty space", 10, 86);

// --- STITCH ---
run("Combine...", "stack1=QC_Capillaries_flat stack2=QC_SMA_flat");
rename("QC_Composite");

saveAs("Tiff", CFG_QC_OUTPUT_DIR + "QC_" + originalName);
logStep("QC image saved: " + CFG_QC_OUTPUT_DIR + "QC_" + originalName);

close();

//============================================================================
// FINAL CLEANUP
//============================================================================
run("Clear Results");
roiManager("reset");
safeClose("empty_space_mask");

run("Close All");

var totalTime = (getTime() - macroStartTime) / 1000;
print("================================================================================");
print("MACRO COMPLETE: " + originalName);
print("Total processing time: " + d2s(totalTime, 1) + " seconds");
print("================================================================================");
print("");
