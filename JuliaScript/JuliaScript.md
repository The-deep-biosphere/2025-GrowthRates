# Julia script for dealing with uncertainties
So now we can run the main part of this whole pipeline: Dealing with the age and abundance uncertainty. This is done using a Julia script. Note that you will need the Project.toml file in the same directory as your script. As well, the following script is tailored to our dataset (file names, column numbers...) and this should be changed to suit other needs.

## Load needed dependencies and data
First we can set up the local environment and load the needed packages. 
```Julia
using Pkg; Pkg.activate("./") # Create a local environment in the current folder.
Pkg.instantiate() # Install packages from the Project.toml file (this may take a while)
# Now we can use the packages (may also take a while the first time the script runs)
#Pkg.update(); # uncomment and specify if any packages have to be updated.
using CSV # read csv files
using DataInterpolations # interpolate to grid
using DataFrames # process data frame
using DelimitedFiles # read delimited files
using GLM # linear models
using Interpolations # interpolate data
using Loess #
using Random # set seed and resample
using StatsBase # basic statistics
using StatsPlots # boxplots
using Tables # process data tables
using UncertainData # handle uncertain data
using XLSX; # read xlsx files
```

Import and process age model data generated in OxCal v.4.4. The ensemble consists of 4224 realisations (columns) based on random walks and is convergent, meaning that adding another realisation will no alter the age distribution for each depth (rows) above a threshold value.

```Julia
# Import age model
age_model = readdlm("./Input/Julia_age_model.csv", Float64);
# Remove 0 and 174 cm horizons as well as row names, which are in the first column.
age_model_c = age_model[2:174, 2:end]; 
# Get median ages for each sediment horizon.
age_c_median = [quantile(age_model_c[i, :], 0.50) for i in 1:size(age_model_c)[1]];
```

Import gene abundance data. `tot16S` and `fungenes` contain the qPCR data for 16S, hzo, and amoA whereas `order` and `otus` contain normalised amplicon reads multiplied by the number of 16S rRNA gene transcripts. 

```Julia
# Load absolute abundance data
# Total 16S abundance
tot16S = readdlm("./Input/Julia_qPCR.csv", ',', Float64, skipstart = 1);
tot16S = tot16S[:, 1:4];
# Order
order = readdlm("./Input/Julia_order_mean.csv", ',', Float64, skipstart = 1);
order_sd = readdlm("./Input/Julia_order_sd.csv", ',', Float64, skipstart = 1);
# Otus
otus = readdlm("./Input/Julia_species_mean.csv", ',', Float64, skipstart = 1);
otus_sd = readdlm("./Input/Julia_species_sd.csv", ',', Float64, skipstart = 1);
# Functional genes
fungenes = readdlm("./Input/Julia_qPCR.csv", ',', Float64, skipstart = 1);
fungenes = fungenes[:, 5:8];
```


video for KDE
https://www.youtube.com/watch?v=t1PEhjyzxLA

