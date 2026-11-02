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

### Input layout

Rows of `x` are the participants' guesses and columns are the arms they were actually assigned to. In the example above, 48 participants in the treatment arm guessed "treatment", 4 guessed "placebo" and 330 answered "don't know".

### James' index

James' index is a single value for the whole trial, ranging from 0 to 1. With the default weights:
- 0: every participant guessed their arm correctly (complete unblinding)
- 0.5: correct and incorrect guesses are balanced and nobody answered "don't know", which is what random guessing would produce
- values towards 1: more participants answered incorrectly or "don't know". Because the index is corrected for chance, it does not always reach exactly 1 even when nobody guesses correctly

Higher values suggest better blinding. Unblinding may be claimed if the upper limit of the two-sided confidence interval is below 0.5 [3].

### Bang's index

Bang's index is calculated separately for each arm and ranges from -1 to 1. It equals the number of correct guesses minus the number of incorrect guesses, divided by the number of participants in that arm, including those who answered "don't know" [2]:
- 0: correct and incorrect guesses are balanced, as expected under random guessing
- 1: every participant in the arm guessed correctly (complete unblinding)
- -1: every participant in the arm guessed incorrectly (opposite guessing)

Values between -0.2 and 0.2 have been used as an ad hoc threshold for successful blinding [3, 4]. Bang et al. preferred one-sided confidence intervals for this index [2]; use `alternative = :greater` to test for unblinding, where a lower limit above 0 suggests participants guessed their arm more often than chance.

A pattern where the treatment arm is positive and the placebo arm is negative by a similar amount (for example 0.4 and -0.4) means most participants in both arms believed they received the active treatment. Bang describes this as "wishful thinking" and does not consider it a failure of blinding [4].

### Example

For the example above:
- James' index is 0.896 (95% CI 0.875 to 0.918). The upper limit is well above 0.5, so the trial-level index gives no evidence of unblinding.
- Bang's index in the treatment arm is 0.115, with a one-sided lower limit of 0.086 (`alternative = :greater`). This is above 0, so participants receiving treatment guessed their arm slightly more often than chance, although the value is within the -0.2 to 0.2 range.
- Bang's index in the placebo arm is 0.022, with a one-sided lower limit of -0.010. There is no evidence of unblinding in this arm.

Most participants (649 of 753) answered "don't know", and this is the main reason James' index is close to 1.

### Limitations

- "Don't know" answers increase James' index. Random guessing with half of the participants answering "don't know" gives a value of 0.75, not 0.5.
- James' index can hide opposite patterns in the two arms. A trial where Bang's index is 0.6 in the treatment arm and -0.6 in the placebo arm still gives a James' index of 0.5. Report Bang's index for each arm alongside James' index.
- Confidence intervals use a normal approximation and are not limited to the range of the index. With small samples they can extend beyond 0 and 1 for James' index or beyond -1 and 1 for Bang's index.
- If every participant answers "don't know", both indices are undefined and the package currently returns `NaN`.
- `bibang` accepts a `weights` argument for consistency with `bijames`, but Bang's index does not use weights.
- Participants may guess correctly because of a noticeable effect or side effect of the treatment, not because the blinding procedure failed. The indices measure participants' beliefs and cannot tell these situations apart.

## Sources

1. James KE, Bloch DA, Lee KK, Kraemer HC, Fuller RK. An index for assessing blindness in a multi-centre clinical trial: disulfiram for alcohol cessation--a VA cooperative study. Stat Med. 1996 Jul 15;15(13):1421-34. doi: 10.1002/(SICI)1097-0258(19960715)15:13<1421::AID-SIM266>3.0.CO;2-H. PMID: 8841652.
2. Bang H, Ni L, Davis CE. Assessment of blinding in clinical trials. Control Clin Trials. 2004 Apr;25(2):143-56. doi: 10.1016/j.cct.2003.10.016. PMID: 15020033.
3. Schwartz M, Mercaldo N. BI: Blinding Indices for Randomized, Controlled, Clinical Trials. R package version 1.2.0. https://CRAN.R-project.org/package=BI
4. Bang H. Random guess and wishful thinking are the best blinding scenarios. Contemp Clin Trials Commun. 2016 Aug 15;3:117-121. doi: 10.1016/j.conctc.2016.05.003. PMID: 27822568