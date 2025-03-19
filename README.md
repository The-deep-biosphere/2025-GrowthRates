# In situ microbial growth rates in deep-sea sediments

## WORK IN PROGRESS!!!!


This repository contains the integrality of the scripts and files needed to reproduce the data analysis presented in the aforementioned article. Note that the authors are in no case Unix/R/Julia professionals, and the code can certainly be written in a more idiomatic way. Do not hesitate to reach out for further help.

## Data
- Raw sequences (fastq files) are deposited on gene bank under the Bioproject [PRJNA784957](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA784957/).
- Other data files at different stages of the processing are present in the Input directories of each processing steps.
- While it is not possible for us to compute the output of the Julia script for all taxonomic groups, we have made it possible to everyone to select a group and process it to investigate its evolution throughout the core. In the [Input directory](./JuliaScript/Input] within the JuliaScript directory, you will find everything needed:
  - The age model (Julia_age_model.csv).
  - The qPCR data (Julia_qPCR.csv).
  - The absolute abundance of each taxa at each taxonomic level (Julia_TAXONOMICLEVEL_mean.csv), along with its uncertainty (Julia_TAXONOMICLEVEL_sd.csv).

## Data processing
All our pipeline is described in the following links. Note that these steps shoudl be run in order, but we also provide the files needed to run each part separately (see Input directories).
- [The processing of the sequences and picking of OTUs](Pipeline%20explanations.md)
- [The decontamination and replicate processing pipelines](./Decontamination_Pooling/DecontaminationPooling.md)
- [The preparation of the data for the Julia script](./PreparationJulia/PreparationJulia.md)
- [The Julia script to process the uncertainties (the main step of this wokflow)](./JuliaScript/JuliaScript.md)
- Plotting the results



Contact: Sven Le Moine Bauer or Steffen Jørgensen
