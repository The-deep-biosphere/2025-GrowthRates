# Julia script for dealing with uncertainties
So now we can run the main part of this whole pipeline: Dealing with the age and abundance uncertainty. This is done using a Julia script. Note that you will need the Project.toml and Manifest.toml files in the same directory as your script. As well, the following script is tailored to our dataset (file names, column numbers...) and this should be changed to suit other needs.

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

## Define necessary parameters and functions
A number of functions are set up and parameters decided to streamline down-stream analysis.

First let us set up the parameters. `n_realisations` is the number of realisations comprising the age model ensemble and used to preallocate appropriately sized output vectors. `n_resamples` sets the number of resamples per sediment depth per age model realisation. `bin_width` sets the width of each time grid cell in years. `age_grid` creates a continuous grid from the lowest to the highest recorded ages in the imported age model by the width of `bin_width`. `lq` and `uq` are the lower and upper quantiles that our uncertainty will represent.

```Julia
n_realisations = size(age_model_c)[2];
n_resamples = 50;
bin_width = 50;
age_grid = floor(Int, minimum(age_model_c)):bin_width:ceil(Int, maximum(age_model_c));
lq = 0.25;
uq = 0.75;
```

We will need a local function for distributing data to grid.
```Julia
function to_grid(t, vals, tgrid, extrapolation_bc = NaN)
    Interpolations.LinearInterpolation(t, vals, extrapolation_bc = extrapolation_bc).(tgrid)
    end;
```

And finally the `get_KDEs()` function. This is the workhorse function of the script. It does a certain amount of things:
- First, it takes an age model realisation, abundance data as well as specifications as to the number of resamples to do and how to distribute them to a time grid.
- Second, it will iterate over each realization of the age model (for loop), and do:
  - It checks whether elements in each realisation are strictly increasing and adds increments if not. If this was not the case, there would be issues when interpolating the abundance data onto the time grid later on.
  - It then resamples the microbial abundance at for each of the 173 samples. 
  - Based on the age value (specific to the given realization), the resampled data is mapped by interpolation onto and age grid.
- Third, it removes empty elements.
- Finally, it creates Kernel Density Estimates (KDEs) to from the data.

KDEs is a very powerful method that smooths discrete data points into a continuous probability distribution, estimating the underlying density function. For more theory on the topic, one can watch this [video](https://www.youtube.com/watch?v=t1PEhjyzxLA).

```Julia
function get_KDEs(agemodel, taxon, n_realisations = n_realisations, n_resamples = n_resamples, bin_width = bin_width)
    # Set parameters
    age_grid = floor(Int, minimum(agemodel)):bin_width:ceil(Int, maximum(agemodel))
    M = zeros(length(age_grid), n_realisations, n_resamples);

    # Run for loop
    for i in 1:n_realisations
        # select age model ensemble
        age_model_i = agemodel[:, i]

        # add incremental value to duplicate age values
        for i in 1:size(age_model_i)[1]
        dups = findall(x -> x == age_model_i[i], age_model_i)
            if size(dups)[1] > 1
                age_model_i[dups[2:end]] = age_model_i[dups[2:end]] .+ 0.0001
            end
        end

        # Add incremental value if diff is negative (only occurs around turbidite due to rounding effects)
        for i in 1:size(age_model_i)[1]-1
            di = diff([age_model_i[i], age_model_i[i+1]]) 
            if di[1] <= 0
                age_model_i[i+1] = age_model_i[i+1] - 100 * di[1]
            end
        end

        # perform resampling
        for j in 1:n_resamples
            realisation_i = resample(taxon, TruncateMinimum(0))
            intp_realisation_i = to_grid(age_model_i, realisation_i, age_grid)
            M[:, i, j] = intp_realisation_i
        end
    end

    # Flatten the 3D array M to a 2D array Mvec
    Mvec = zeros(length(age_grid), n_realisations*n_resamples);
    for i in 1:size(Mvec)[1]
        Mvec[i,:] = M[i,:,:]
    end
    M = 0

    # Remove NaNs resulting from empty elements
    Nvec = [Mvec[i, isnan.(Mvec[i, :]) .== 0] for i in 1:size(Mvec)[1]];
    Mvec = 0

    # Remove the first bin if empty
    if size(Nvec[1])[1] == 0
        Nvec = Nvec[2:end]
    end;

    # Convert to and return Kernel Density Estimates
    [UncertainValue(UnivariateKDE, x) for x in Nvec]
end;
```

## Producing the KDEs
Now we are ready to produce the data. As mentionned earlier, any dataset which contains a mean and a standard deviation can be used. Here we will only focus on the OTU_3 shown in Figure 2.

```Julia
# We can set a seed to allow replication.
Random.seed!(123)
# Select the abundance data for OTU_3
OTU_3 = [otus[:, 50] otus_sd[:, 50]]
# Abundance data is transformed to probability distributions, i.e., uncertain values, assuming that the uncertainty is normally distributed.
OTU_3 = [UncertainValue(Normal, OTU_3[row, 1], OTU_3[row, 2], trunc_lower = 0.0001) for row in 1:size(OTU_3, 1)];
# Compute the KDEs
KDEs_OTU_3 = get_KDEs(age_model_c, OTU_3);
```

From the output, we can obtain the median and quantiles.
```Julia
meds_OTU_3 = median.(KDEs_OTU_3)
lqs_OTU_3 = quantile.(KDEs_OTU_3, lq)
uqs_OTU_3 = quantile.(KDEs_OTU_3, uq);

# Remove objects to save memory
OTU_3 = 0; KDEs_OTU_3 = 0;
```
We also want to replace any negative value by 1 due to the log scale on future plots.
```Julia
meds_OTU_3[meds_OTU_3 .< 1] .= 1;
uqs_OTU_3[uqs_OTU_3 .< 1] .= 1;
lqs_OTU_3[lqs_OTU_3 .< 1] .= 1;
```
Good, now we can export the data for plotting in R (last script).
```Julia
df = DataFrame(
    Age = age_grid[1:end-1],
    OTU_3_Median = meds_OTU_3,
    OTU_3_LowerQuantile = lqs_OTU_3,
    OTU_3_UpperQuantile = uqs_OTU_3,
)
CSV.write("./FILENAME.csv", df);
```
## Computing growth/decay rates
Here we are doing this manually, because many parts of the plots are not reliable (remember all these oscillations of values between 0 and ~100 in figure 2) we cannot just ask where is the highest growth/decay rate in the whole plot. So first lets choose the section that we want to investigation. For the groath of OTU_3, this is roughly the section between 35K and 37K years old. We need then to look at which sections of the age grid this corresponds. For getting this info, run the following code:
```Julia
ageID = for (i, val) in enumerate(age_grid)
    println(i, ": ", val)
end
show(ageID)
```
In our case this corresponds to positions 703-743. Let us plot these part:
```Julia
# Define the interval of investigation
sed_age_itvl = 703:743

# Set up the plot
plot_OTU_3 = plot(
    age_grid[sed_age_itvl], meds_OTU_3[sed_age_itvl], 
    ribbon = (meds_OTU_3[sed_age_itvl] .- lqs_OTU_3[sed_age_itvl], 
              uqs_OTU_3[sed_age_itvl] .- meds_OTU_3[sed_age_itvl]), 
    yscale = :log10, ylabel = "Copies g⁻¹ sed.", title = "OTU_3", 
    fillalpha = 0.4, fc = :green, c = :green, legend = false
)
```
![image](https://github.com/user-attachments/assets/0a1125ea-345d-4121-9705-ecc9520d80fb)

Now we can look at the doubling time for the same period:
```Julia
growth_rates = 1 ./ ((log2.(meds_Sediminis[2:end]) .- log2.(meds_Sediminis[1:end-1])) ./ bin_width)

plot(age_grid[sed_age_itvl], growth_rates[sed_age_itvl], 
     ylims = (-1000, 1000),
     ylabel = "years", xlabel = "", title = "Doubling times", 
     legend = false, c = :red)
```
![image](https://github.com/user-attachments/assets/40ce1639-2cf2-419c-8ff7-cfb92a9d1e8d)

This plot is a bit counterintuitive to interpret. We are looking at doubling times, meaning that anything positive means growth, and anything negative means decay. However, the closest to 0 means growth or decay happens the fastest. Therefore, the jumps from extremely positive values to extremely negative values represent a transition from a slow growth to a slow decay. In our case we want to find the lowest positive value around 35500 years old. Let's print the doubling time. By looking again at our age grid we can refine a bit better the interval of investigation to 710-720.
```Julia
println(1 ./ ((log2.(meds_Sediminis[710:720]) .- log2.(meds_Sediminis[709:719])) ./ bin_width))
```
In position 4, we get 62 years doubling time. One could argue that there is a 38 years in position 1, but looking at the plot it looks a bit like a first hic up before the real growth.
Voila. This approach can be applied to each curve, both for growth and decay (then one would look for the smallest negative value).
