plop


```Julia
# Set up local environment. 
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
# using Plots; pyplot() # use pyplots backend to enable log scales.
using Random # set seed and resample
using StatsBase # basic statistics
using StatsPlots # boxplots
using Tables # process data tables
using UncertainData # handle uncertain data
using XLSX; # read xlsx files
```
video for KDE
https://www.youtube.com/watch?v=t1PEhjyzxLA
