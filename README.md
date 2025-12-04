# AssiST — User Guide

This guide explains how to prepare training images, train a convolutional neural network (CNN), optionally explore hyperparameters to choose an optimal configuration, and then infer **MIC** (minimum inhibitory concentration) and **MGC** (maximum growth concentration) from full-plate images.

> Overall workflow  
> **prepare → train → infer (MIC/MGC) → review outputs**
All scripts will be run in MATLAB 

---

## 1) Prerequisites

- **MATLAB** (R2023a or newer recommended).  
- **Toolboxes**
  - Deep Learning Toolbox
  - Image Processing Toolbox
  - Parallel Computing Toolbox *(optional; speeds up exploration)*

### Directory layout used by the project:

```
/prepare/                      # tools to build the training set from plate images
    ├─── full_plate_images/    # user provided: full 96 well plate images for training, in grayscale or RGB (RGB will be converted to grayscale)
    ├─── single_wells/         # auto-created; cropped single wells from full plate images
    ├─── training_sets/        # user provided: cropped single well images in folders corresponding to class 
    ├─── training_setsRotate/  # auto-created; rotated versions of images in training_sets in classes folders 
/train/                        # single-network training
/infer/                        # classify new plates, compute MIC/MGC, export outputs, optional hyperparameter search 
    ├─── full_plate/           # user provided: full 96 well plate images for inferring MIC and MGC, in grayscale or RGB (RGB will be converted to grayscale) 
    ├─── PhenotypeImages/      # auto-created; overlays of per-well phenotypes
    ├─── QAimages/             # auto-created; MIC/MGC overlays on plate
    ├─── MetricImages/         # auto-created; MIC/MGC overlaid on concentration map
``` 
### File index (where things live)
- **Prepare**
  - `prepare/prepareTrainingSet.m` - batch prepare plates
  - `prepare/cropSingleWells.m` - per-plate crop single wells 
  - `prepare/rotateTrainingImages.m` - per single well image rotate and put in class folders 

- **Training**
  - `train/constructSingleNetwork.m` — single-run trainer; saves `singleNetwork.mat`.
  - `train/trainEvaluateNetwork.m` — core model+training routine.
  - `train/loadDataAndPreprocess.m` — splits data & attaches the reader.
  - `train/customImageReader.m` — reads any image → **50×50×1** single in [0,1].
  - `train/explorerCNNparameters.m` — grid search + reporting.
  - `train/extractParameters.m` — maps a linear index to grid params.
  - `train/parsave.m` — saving in `parfor`.
 

- **Inference**
  - `infer/classifyFullPlates.m` — batch classify plates; saves all artifacts.
  - `infer/inferPlatePhenotypes.m` — per-plate infer (crop→classify→MIC/MGC).
  - `infer/inputDialog2.m` — reclassification UI with validation.



---

## 2) Data preparation: Creating training images

Overview: Training uses a standard **folder-per-class** layout. The helper `loadDataAndPreprocess.m` builds `imageDatastore`s with an image reader that converts everything to **grayscale 50×50×1**, single-precision in [0,1].

Directory layout used by this step:

```
/prepare/                    # tools to build the training set from plate images
    ├─── full_plate_images/    # user provided: full 96 well plate images in grayscale or RGB (RGB will be converted to grayscale)
    ├─── single_wells/         # auto-created; cropped single wells from full plate images
    ├─── training_sets/        # user provided: cropped single well images in folders corresponding to class 
        ├─ 1_bubbles/
        ├─ 2_empty/
        ├─ 3_pellet_tiny/
        ├─ 4_cloud_small/
        ├─ 5_pellet_small/
        ├─ 6_cloud/
        └─ 7_pellet/
    ├── training_setsRotate/  # auto-created; rotated versions of images in training_sets in classes folders 
        ├─ 1_bubbles/
        ├─ 2_empty/
        ├─ 3_pellet_tiny/
        ├─ 4_cloud_small/
        ├─ 5_pellet_small/
        ├─ 6_cloud/
        └─ 7_pellet/    
```

### Protocol (you can adapt to your naming / capture setup):

1. All steps will be completed in "prepare" folder 
2. Put raw full-plate photos in a working folder 
3. Run `prepareTrainingSet.m` to find well locations, and crop image to single wells.
4. Create class folders below. Class folder must be named in this format: "#_classname"

```
prepare/training_sets/
├─ 1_bubbles/
├─ 2_empty/
├─ 3_pellet_tiny/
├─ 4_cloud_small/
├─ 5_pellet_small/
├─ 6_cloud/
└─ 7_pellet/
```
5. Move cropped images from "single_wells" folder into the class folders above.   
6. *(Optional)* Augment or balance classes as needed.
7. Run `rotateTrainingImages.m` to rotate the images from "training_sets" to create a larger training set. 
8. Final training set for training neural network is in "training_setsRotate" folder 

---

## 3) Train a single network (recommended path)

Overview: Trains a single network with configuration set by the user 

### Protocol:
1. All steps will be completed in "train" folder  
2. Edit `constructSingleNetwork.m`: At the top of the file, set the hyperparameters you want:

```
% Filters per conv block (enforce n1 <= n2 <= n3)
numFilters1 = 4;
numFilters2 = 4;
numFilters3 = 16;

filterSize  = 5;      % e.g. 3 or 5
learnRate   = 0.01;   % e.g. 1e-2 or 1e-3

trainingLocation = "../prepare/training_setsRotate/";
outputFileName   = "singleNetwork.mat";
```

3. Run `constructSingleNetwork.m`

The script will then:

- Split the dataset into **train/validation** and read using the robust 50×50×1 reader.  
- A small 3‑block CNN is constructed and trained (SGDM optimizer; early validation checks).  
- `singleNetwork.mat` is saved (contains `net` and `info`), along with a confusion chart image for quick sanity‑check of class performance.
- The training code seeds the RNG (`rng(0)`) for reproducibility of splits and initialization.

### Under the hood (architecture summary)

- Input: **50×50×1**
- 3× blocks of **Conv → BatchNorm → ReLU** (with pooling after the first two blocks)
- **Dropout** before the classifier head
- FC(hidden) → ReLU → FC(number of classes) → Softmax → Classification layer
- Advanced knobs (optional name‑value pairs): `MiniBatchSize`, `L2`, `MaxEpochs`, `ValidationPatience`, `CheckpointPath`.

---

## 4) (Optional) Explore many configurations

If you’re unsure about hyperparameters, use `train/explorerCNNparameters.m`. It:

- Defines a small grid (filter counts, filter sizes, learning rates).  
- Trains each configuration (optionally in **parallel**).  
- Collects **training/validation** metrics, computes **external accuracy** on the validation set, and saves each model as `train/Networks/net_XX.mat`.  
- Produces a **results.mat** and **optimalNetwork.mat** with top model, plus a figure summarizing compute vs. accuracy.

Run:

```matlab
cd train
explorerCNNparameters
load(`optimalNetwork.mat`)
% Edit `constructSingleNetwork.m` to have hyperparameters of optimalNetwork that is found in `optimalNetInfo.Config`
constructSingleNetwork 
```

Outputs include `train/Networks/net_XX.mat` for each candidate, a `results.mat` table, and `train/optimalNetwork.mat` containing the best network and training info. You can then use that network for inference instead of `singleNetwork.mat`.

> Notes:
> - Only monotonically non‑decreasing filter counts (n1 ≤ n2 ≤ n3) are trained.  
> - `parsave` is used to safely save inside `parfor`.  
> - External accuracy helper computes mean accuracy over the validation datastore.

---

## 5) Infer MIC and MGC from new plates

### Required in "infer" folder:
1) Full 96 well plate images (`.png`/`.jpg`): grayscale or RGB accepted, pipeline will convert RGB to grayscale 
2) **(Optional)** Plate map Excel in layout:
    - **Drug names** in cells **B3:M10** (8 rows × 12 columns).
    - **Concentrations** in cells **B14:M21** (same shape, numeric).
    - **Drug types** in cells **B25:M32** (same shape). Must be 'bacteriostatic' or 'bactericidal'. 


### Protocol:
1) All steps will be completed in "infer" folder 
2) Use `classifyFullPlates.m`. Configure the **user definitions** at the top:

```matlab
platemapFileName = 'platemap.xlsx';  % optional: valid platemap - Excel map (see layout above); If MIC/MGC inference is unneeded then put "none"
dataFolder       = './DataFolder/';       % folder of full-plate images (.png/.jpg)
netFileName      = '../train/singleNetwork.mat';  % or path to your optimal network
resultsMatName   = "results.mat";    % output filename

% Orientations to get A1 at top-left and H12 at bottom-right
tfFlip   = true;     % horizontal flip first
tfRotate = true;     % then 90° rotate

% Growth classifications for class 
% classes= [1_bubbles, 2_empty, 3_pellet_tiny, 4_cloud_small, 5_pellet_small, 6_cloud, 7_pellet]
numNoGrowth = [1,2,3]; % number of classes corresponding to no growth 
numRestricted = [4,5]; % number of classes corresponding to restricted growth
numFull = [6,7]; % number of classes corresponding to robust growth

% OPTIONAL: Flag wells with dubious classification
% This is the confidence threshold for probability of class inference.
% Anything below # will be flagged.
pClassThreshold = nan; % Put # between 0-1.Put "nan" if you want no flag.
```

3) Run `classifyFullPlates.m`
You will be asked to **click A1 and H12** on the displayed plate image (upper‑left and bottom‑right wells). Top-left well should be A1 and bottom-right should be H12. If orientation of the plate looks wrong then change tfFlip and tfRotate so that the correct orientation is produced.

The script will then:

1. Compute well centers, crop each well, and classify it with your CNN (robust 50×50×1 input).  
2. Render per‑well overlays (phenotype color and class index).  
3. **Optionally flag** any wells that have dubious classification
4. **Optionally reclassify** any wells via a modal dialog (see next section).  
5. **Optional if valid platemap is provided:** Compute **MGC** (maximum growth concentration) and **MIC** (first concentration with no growth) per drug line, handling restricted growth and skip‑well cases.  
5. Export artifacts:  
   - **`PhenotypeImages/`** — per‑well phenotype overlays (after reclassification). 
   - **Optional if valid platemap is provided**  
      - **`QAimages/`** — plate overlays with **MIC (brown)** and **MGC (black)** markings.  
      - **`MetricImages/`** — concentration map with MIC/MGC markings.  
   - **`results.mat`** — contains:
      - `phenotypeImg` with per‑file details.
        **Optional if valid platemap is provided**  - also contains: 
          - `MICtable` and `MGCtable` (rows=drugs, columns=files)
          - raw arrays `drugMICs`, `drugMGCs`


**Class legend (indexed 1… # of classes):**  
1) bubbles, 2) empty, 3) pellet_tiny, 4) cloud_small, 5) pellet_small, 6) cloud, 7) pellet

---

## 6) Manual reclassification dialog

During inference you’ll see a small window to correct any misclassified wells. It shows the **class legend**, and lets you enter:

- **Well IDs** (e.g., `A1 B3 H12`; comma/space/semicolon delimited), and
- **Class IDs** (`1`–`# of classes`). Enter one value to apply to all wells, or a list matching the well count.

Option:
- **Accepts empties** — leaving **both fields empty** and pressing OK skips reclassification.  

---

## 7) Quick start examples

### A. Train a single network and run inference

```matlab
% 0) Prepare training data under prepare/training_setsRotate/<class>/
% 1) Train
cd train
edit constructSingleNetwork.m % set neural network configuration, trainingLocation, outputFilename 
constructSingleNetwork

% 2) Infer MIC/MGC on new plates
cd ../infer
edit classifyFullPlates.m   % set platemapFileName, dataFolder, netFileName as needed
classifyFullPlates
```

### B. Explore hyperparameters before training

```matlab
% 0) Prepare training data under prepare/training_sets/<class>/
% 1) Explore a small grid of CNNs
cd train
explorerCNNparameters

% 2) Use the best saved network for inference
cd ../infer
netFileName = '../train/optimalNetwork.mat';  % or a specific ./Networks/net_XX.mat
classifyFullPlates
```

---

## 8) Troubleshooting & tips

- **Input size/type errors**: ensure readers produce **50×50×1** singles in [0,1]; when you build your own `imageDatastore`, set `ReadFcn=@customImageReader`.
- **Plate orientation wrong**: toggle `tfFlip`/`tfRotate` so that **A1 is top‑left** and **H12 is bottom‑right** before clicking points.
- **Low validation accuracy**: use the *explore* step, rebalance classes, or tweak `learnRate`, `filterSize`, and filter counts.
- **Parallel errors / saving**: exploration uses `parsave` to save `net`/`info` inside `parfor`. If you see permission issues, check the `./Networks` folder.
- **Reproducibility**: training scripts call `rng(0)` so results are consistent across runs (given the same data split).

---

## 9) Glossary

- **Phenotype classes (# of classes):** visual growth patterns per well.  
- **MIC (Minimum Inhibitory Concentration):** lowest concentration where **no growth** is observed.  
    - Determine based on the EUCAST guidelines
      - bacteriocidal drugs:lowest concentration where **any growth** (restricted or robust growth) is observed
      - bacteriostatic drugs: lowest concentration where **no robust growth** is observed
- **MGC (Maximum Growth Concentration):** highest concentration at which growth is **robust**
    - **If highest concentration well with growth is:**  
        - **Classified as Robust Growth** : concentration of highest concentration with robust growth 
        - **Classified as Restricted Growth** : concentration between lowest concentration with restricted growth and highest concentration with robust growth
- **Skip:** will use well with concentration below skipped well to determine MIC and MGC
    -**Best practice** is to redo the broth microdilution for drug with skipped well
---



