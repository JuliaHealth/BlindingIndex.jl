# BlindingIndex.jl

Compute James' and Bang's Blinding Indices.

The code is based on R package [BI](https://github.com/marcschwartz/BI) by Marc Schwartz and Nate Mercaldo.

## Usage

```julia
using BlindingIndex
x = [48 22; 4 30; 330 319]
bi(x)
```

Output:

```
BI James     Estimate     Std. Error   LCL      UCL      2-Sided
Overall      0.896        0.011        0.875    0.918   

BI Bang      Estimate     Std. Error   LCL      UCL      2-Sided
Treatment    0.115        0.018        0.08     0.15    
Placebo      0.022        0.019        -0.016   0.06    
```

## Arguments

- `x::Matrix{Int64}`: 2×2 or 3×2 integer matrix of cross-tabulated counts
- `weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1]`: use default 1996 James weights for 3×2 unless alternative weights are specified; correct guesses are assigned a weight of 0, incorrect guesses are assigned a weight of 0.5, don't know guesses are assigned a weight of 1
- `conf::Float64=0.95`: confidence interval
- `alternative::Symbol=:two`: whether two-sided (`:two`) (the default) or one-sided (`:less` or `:greater`) confidence intervals are returned
- `groups::Vector{String}=["Treatment", "Placebo"]`: treatment group names
- `output::Bool=true`: if true, show the output

## Returns

Named tuple containing:
- `bi_james`: James' index (overall estimated value, standard error, lower confidence interval bound, upper confidence interval bound)
- `bi_bang`: Bang's index (estimated value, standard error, lower confidence interval bound, upper confidence interval bound for groups 1 and 2)

## Notes

The format of the cross-tabulated counts and weights is:

```
             Treatment  Placebo
Treatment    xxx        xxx 
Placebo      xxx        xxx 
Don't Know   xxx        xxx
```

If 2×2 cross-tabulated counts matrix is provided, the last row is filled with zeros. 

## Sources

1. James KE, Bloch DA, Lee KK, Kraemer HC, Fuller RK. An index for assessing blindness in a multi-centre clinical trial: disulfiram for alcohol cessation--a VA cooperative study. Stat Med. 1996 Jul 15;15(13):1421-34. doi: 10.1002/(SICI)1097-0258(19960715)15:13<1421::AID-SIM266>3.0.CO;2-H. PMID: 8841652.
2. Bang H, Ni L, Davis CE. Assessment of blinding in clinical trials. Control Clin Trials. 2004 Apr;25(2):143-56. doi: 10.1016/j.cct.2003.10.016. PMID: 15020033.