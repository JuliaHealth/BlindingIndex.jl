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
BI James     Estimate     Std. Error   LCL      UCL      95% CI, 2-sided
Overall      0.896        0.011        0.875    0.918   

BI Bang      Estimate     Std. Error   LCL      UCL      95% CI, 2-sided
Treatment    0.115        0.018        0.08     0.15    
Placebo      0.022        0.019        -0.016   0.06    
```

## Arguments

- `x::Matrix{Int64}`: 2×2 or 3×2 integer matrix of cross-tabulated counts
- `weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1]`: use default 1996 James weights for 3×2 unless alternative weights are specified; correct guesses are assigned a weight of 0, incorrect guesses are assigned a weight of 0.5, don't know guesses are assigned a weight of 1
- `conf::Float64=0.95`: confidence level for the returned confidence intervals
- `alternative::Symbol=:two`: whether two-sided (`:two`) (the default) or one-sided (`:less` or `:greater`) confidence intervals are returned
- `groups::Vector{String}=["Treatment", "Placebo"]`: treatment group names
- `output::Bool=true`: if true, show the output

## Returns

Named tuple containing:
- `bi_james`: James' index (overall estimated value, standard error, lower confidence interval bound, upper confidence interval bound)
- `bi_bang`: Bang's index (estimated values, standard errors, lower confidence interval bounds, upper confidence interval bounds for group 1 and 2)

## Notes

The format of the cross-tabulated counts and weights is:

```
             Treatment  Placebo
Treatment    xxx        xxx 
Placebo      xxx        xxx 
Don't Know   xxx        xxx
```

If 2×2 cross-tabulated counts matrix is provided, the last row is filled with zeros.

## Interpreting the results

In the cross-tabulated counts, rows are the participants' guesses and columns are the arms they were actually assigned to.

James' index summarises blinding for the whole trial and ranges from 0 to 1:
- 1 means complete blinding (every participant answered "don't know")
- 0.5 means participants guessed no better than chance
- 0 means complete unblinding (every participant guessed correctly)

If the upper confidence limit is below 0.5, the trial can be considered unblinded.

Bang's index is calculated separately for each arm and ranges from -1 to 1:
- 0 means guesses in that arm were consistent with random guessing
- 1 means every participant in that arm guessed their assignment correctly
- -1 means every participant in that arm guessed the opposite assignment

A confidence interval that excludes 0 suggests that participants in that arm were more likely to guess correctly (positive values) or incorrectly (negative values) than chance.

In the example above, James' index of 0.896 indicates the trial was well blinded overall. Bang's index is close to 0 in the placebo arm, while the treatment arm shows a slight tendency towards correct guessing (0.115, 95% CI 0.08 to 0.15).

## Sources

1. James KE, Bloch DA, Lee KK, Kraemer HC, Fuller RK. An index for assessing blindness in a multi-centre clinical trial: disulfiram for alcohol cessation--a VA cooperative study. Stat Med. 1996 Jul 15;15(13):1421-34. doi: 10.1002/(SICI)1097-0258(19960715)15:13<1421::AID-SIM266>3.0.CO;2-H. PMID: 8841652.
2. Bang H, Ni L, Davis CE. Assessment of blinding in clinical trials. Control Clin Trials. 2004 Apr;25(2):143-56. doi: 10.1016/j.cct.2003.10.016. PMID: 15020033.