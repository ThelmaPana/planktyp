# Planktyp

Typology of Plankton Communities seen by *In Situ* Imaging in the First 500 m of the Global Ocean

Thelma Panaïotis M2 internship
January 2019 - June 2019
Laboratoire d’Océanographie de Villefranche (UMR 7093)

## TODOs
- [ ] check Trichodesmium below 500 m
- [ ] check high watervolumes
 
## Organisation of repository
### Folders
- `lib` contains scripts called when needed
    - `lib/extract_zoo` contains scripts to extract plankton data
- `data` contains data
    - `data/extract_zoo` contains data needed to run plankton data extraction
    - `data/raw` contains other raw data needed to run the analysis
- `plots` contains generated plots, organized by theme

### Scripts
Scripts order is self-explanatory. 
- `00.run_zoo_extraction.R` run extraction of biological data from EcoTaxa (see https://github.com/jiho/UVP5_images_dataset)
- `01.zoo_data_split.R` split biological data into zooplankton, detritus and particles
- `02.ctd_extraction.R` extract CTD data from Ecopart 
- `03.ctd_process.R` process CTD data (clean, despike, smooth, compute clines)
- `04.satellite_process.R` match satellite data with CTD data
- `05.assign_provinces.R` match profiles with partitioning to test (latitude bands, Longhurt provinces, mesopelagic provinces)
- `06.merge_bio_env.R` match biological and environmental data and cut layers
- `07.analysis_epi.R` analyses for epipelagic layer
- `08.analysis_mesosup.R` analyses for mesopelagic superior layer
- `09.analysis_mesoinf.R` analyses for mesopelagic inferior layer
- `10.analysis_bathy.R` analyses for bathypelagic layer
- `11.additional_plots.R` additional and exploratory plots


## Results
See paper: https://docs.google.com/document/d/1rGzqq_zt_J82GPJp4djq53LQi6AjiLVa5TquGno2_h8
