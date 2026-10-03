# DG–CA3 connectoid analysis

MATLAB for MaxWell HD-MEA recordings of paired organoids (connectoids) on one array. 
Data was recorded using Maxwell Biosystems MaxOne system
The data from the Maxwell exports four file folders:

1. Activity scan (map all electrodes) eg. /000222
2. Network scan (spontaneous recording of allocated electrodes) eg. /000223
3. Stimulation of organoid 1, record both eg. /000224
4. Stimulation of organoid 2, record both eg. /000225

| launcher | input | output |
|---|---|---|
| `connectoid_launcher.m` | one spontaneous network scan | spatial units, firing rates, putative monosynaptic connections, inter-organoid synchrony |
| `stim_launcher.m` | the two stimulation recordings of one pair | unstimulated organoid responses |

## Requirements

MATLAB R2025a with the Statistics and Machine Learning Toolbox.

**MaxWell MATLAB toolbox 23.2** is required and is not redistributed here. Can be requested through the Maxwell Knowledge Hub. Put it on your MATLAB path

The Maxwell decompression libraries are required to read the raw voltage traces. This repo operates only on the exported spike times from h5 files. 

## Workflow

### Spontaneous recordings

Open `connectoid_launcher.m`, set the path to the Network Scan folder and run the sections in order.

```matlab
config.DataFile = ''; %path to network scan .h5 file or folder
```

Putative units in this analysis are from spatial clustering of electrodes with firing rates above MinElectrodeFiringRateHz and are not spike sorted. Spike sorting could be implemented if desired. 

A directed pair is accepted on a cross-correlogram peak at short positive lag clearing both a minimum count and a multiple of its own off-peak baseline, 0 lag peaks are rejected (window = +1 to +11 ms). 

### Stimulation recordings

Open `stim_launcher.m`, set the two paths to the Stimulation folders and the roi.csv and run the sections in order. The stimulated organoid is detected from the data. 

```matlab
config.StimulateFirst = ''; %path to first stim file or folder
config.StimulateSecond = ''; %path to second stim file or folder

config.OrganoidRoiFile = ''; %path to ROI file if already ran connectoid_launcher
```

Outputs are computed from the post-stimulus blank window (30-100 ms):

- **the window mean**, % change from that organoid's own pre-stimulus
  baseline, with a connected-or-not call from whether its interval clears the
  null's;
- **% of trials responding**, the share of trials whose own window rate exceeds
  the 95th percentile of that organoid's baseline windows. Its no-effect value
  is 5 by construction, which makes it comparable across connectoids with very
  different absolute rates.

The null is every window position on a millisecond grid across the stimulus-free baseline, resampled at the trial count.

## Output

Each run creates a results folder next to the recording — figures, CSV tables, a MAT
file with every setting, and a copy of the launcher that produced it — and
appends one row to `connectoid_spontaneous.csv` or
`connectoid_stimulation.csv` one level above. Set `config.TableFile` to collect
recordings from different folders into one table. Reanalysing a recording
replaces its row.

## Tests

```matlab
addpath(genpath('analysis'));
runtests('analysis/tests');
```

`testStimulationPipeline` runs the whole stimulation chain on a synthetic
connectoid in a few seconds, without a recording.

## Layout

```text
DG_CA3_Analysis/
├── connectoid_launcher.m      one spontaneous recording
├── stim_launcher.m            one connectoid, both stimulation directions
├── analysis/
│   ├── shared/                used by both launchers
│   ├── stim/                  stimulation analysis and figures
│   ├── spontaneous/           units, connectivity, synchrony
│   └── tests/
└── Maxwell/                   not in this repository — see Requirements
```
