module BlindingIndex

export bi
export bijames
export bibang
export addmargins
export nrow
export ncol

using Distributions

function addmargins(x::AbstractMatrix)
    n = zeros(Int64, size(x) .+ 1)
    n[1:size(x, 1), 1:size(x, 2)] = x
    n[end, 1:(end - 1)] = sum(x, dims=1)
    n[1:(end - 1), end] = sum(x, dims=2)
    n[end, end] = sum(x)
    return n
end

nrow(x::AbstractMatrix) = size(x, 1)
ncol(x::AbstractMatrix) = size(x, 2)

"""
    bi(x; <keyword arguments>)

Compute James' and Bang's Blinding Indices.

# Arguments

- `x::Matrix{Int64}`: 2×2 or 3×2 integer matrix of cross-tabulated counts
- `weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1]`: use default 1996 James weights (correct guesses are assigned a weight of 0, incorrect guesses are assigned a weight of 0.5, don't know guesses are assigned a weight of 1) unless alternative weights are specified
- `conf::Float64=0.95`: confidence level for the returned confidence intervals
- `alternative::Symbol=:two`: type of returned confidence interval (two-sided, `:two`) or one-sided (`:less` or `:greater`)
- `groups::Vector{String}=["Treatment", "Placebo"]`: treatment group names
- `output::Bool=true`: if true, show the output

# Returns

Named tuple containing:
- `bi_james`: James' index: overall estimated value, standard error, lower confidence interval bound, upper confidence interval bound
- `bi_bang`: Bang's index: estimated values, standard errors, lower confidence interval bounds, upper confidence interval bounds for group 1 and 2

# Notes

The format of the cross-tabulated counts and weights is:

             Treatment  Placebo
Treatment    xxx        xxx 
Placebo      xxx        xxx 
Don't Know   xxx        xxx

If 2×2 cross-tabulated counts matrix is provided (treatment and placebo rows), the last row (don't know) is filled with zeros.

# Sources

1. James KE, Bloch DA, Lee KK, Kraemer HC, Fuller RK. An index for assessing blindness in a multi-centre clinical trial: disulfiram for alcohol cessation--a VA cooperative study. Stat Med. 1996 Jul 15;15(13):1421-34. doi: 10.1002/(SICI)1097-0258(19960715)15:13<1421::AID-SIM266>3.0.CO;2-H. PMID: 8841652.
2. Bang H, Ni L, Davis CE. Assessment of blinding in clinical trials. Control Clin Trials. 2004 Apr;25(2):143-56. doi: 10.1016/j.cct.2003.10.016. PMID: 15020033.
"""
function bi(x::Matrix{Int64}; weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1], conf::Float64=0.95, alternative::Symbol=:two, groups::Vector{String}=["Treatment", "Placebo"], output::Bool=true)

    @assert alternative in [:two, :less, :greater] "alternative must be :two, :less, or :greater"
    @assert length(groups) == 2 "groups must contain 2 names"

    size(x) == (2, 2) && (x = vcat(x, [0 0]))
    @assert size(x) == (3, 2) "x must be 2×2 or 3×2 integer matrix of cross-tabulated counts"

    jbi = bijames(x, weights=weights, conf=conf, alternative=alternative)
    bbi = bibang(x, weights=weights, conf=conf, alternative=alternative)

    if output
        sided = alternative === :two ? "2-sided" : "1-sided"

        println("$(rpad("BI James", 12)) $(rpad("Estimate", 12)) $(rpad("Std. Error", 12)) $(rpad("LCL", 8)) $(rpad("UCL", 8)) $(round(Int64, conf * 100))% CI, $sided")
        println("$(rpad("Overall", 12)) $(rpad(jbi.bi_est, 12)) $(rpad(jbi.bi_se, 12)) $(rpad(jbi.bi_lcl, 8)) $(rpad(jbi.bi_ucl, 8))")
        println()
        println("$(rpad("BI Bang", 12)) $(rpad("Estimate", 12)) $(rpad("Std. Error", 12)) $(rpad("LCL", 8)) $(rpad("UCL", 8)) $(round(Int64, conf * 100))% CI, $sided")
        println("$(rpad(groups[1], 12)) $(rpad(bbi.bi_est[1], 12)) $(rpad(bbi.bi_se[1], 12)) $(rpad(bbi.bi_lcl[1], 8)) $(rpad(bbi.bi_ucl[1], 8))")
        println("$(rpad(groups[2], 12)) $(rpad(bbi.bi_est[2], 12)) $(rpad(bbi.bi_se[2], 12)) $(rpad(bbi.bi_lcl[2], 8)) $(rpad(bbi.bi_ucl[2], 8))")
    end

    return (bi_james=jbi, bi_bang=bbi)

end

"""
    bijames(x; <keyword arguments>)

Compute James' Blinding Indices.

# Arguments

- `x::Matrix{Int64}`: 2×2 or 3×2 integer matrix of cross-tabulated counts
- `weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1]`: use default 1996 James weights (correct guesses are assigned a weight of 0, incorrect guesses are assigned a weight of 0.5, don't know guesses are assigned a weight of 1) unless alternative weights are specified
- `conf::Float64=0.95`: confidence level for the returned confidence intervals
- `alternative::Symbol=:two`: type of returned confidence interval (two-sided, `:two`) or one-sided (`:less` or `:greater`)

# Returns

Named tuple containing:
- `bi_est::Float64`: estimated value
- `bi_se::Float64`: standard error
- `bi_lcl::Float64`: lower confidence interval bound
- `bi_ucl::Float64`: upper confidence interval bound

# Notes

The format of the cross-tabulated counts and weights is:

             Treatment  Placebo
Treatment    xxx        xxx 
Placebo      xxx        xxx 
Don't Know   xxx        xxx

If 2×2 cross-tabulated counts matrix is provided, the last row is filled with zeros.

# Sources

1. James KE, Bloch DA, Lee KK, Kraemer HC, Fuller RK. An index for assessing blindness in a multi-centre clinical trial: disulfiram for alcohol cessation--a VA cooperative study. Stat Med. 1996 Jul 15;15(13):1421-34. doi: 10.1002/(SICI)1097-0258(19960715)15:13<1421::AID-SIM266>3.0.CO;2-H. PMID: 8841652.
"""
function bijames(x::Matrix{Int64}; weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1], conf::Float64=0.95, alternative::Symbol=:two)

    @assert alternative in [:two, :less, :greater] "alternative must be :two, :less, or :greater"

    bi_lcl = 0.0
    bi_ucl = 0.0

    @assert size(weights) == (3, 2) "weights must be a 3 row by 2 column integer matrix specifying alternative James weights for each cell in x"

    @assert conf > 0.0 && conf < 1.0 "conf must be > 0 and < 1"

    if x[1:2, :] == [0, 0] && all(i -> i > 0, x[3, :])
        bi_est = 1
        bi_se = 0
    else
        x1 = addmargins(x)
        p = x1 / maximum(x1)
        pdk = p[nrow(p) - 1, ncol(p)]
        pdo = 0
        pde = 0
        v = 0
        term1_denom = 0
        for i in 1:(nrow(p) - 2)
            for j in 1:(ncol(p) - 1)
                pdo += ((weights[i, j] * p[i, j]) / (1 - pdk))
                pde += ((weights[i, j] * p[i, ncol(p)] * (p[nrow(p), j] - p[nrow(p) - 1, j])) / (1 - pdk) ^ 2)
                term1_denom += weights[i,j] * p[i, ncol(p)] * (p[nrow(p), j] - p[nrow(p) - 1, j])
            end
        end
        kd = (pdo - pde) / pde
        term1_denom = 4 * term1_denom ^ 2
        term1_num = 0
        for i in 1:(nrow(p) - 2)
            for j in 1:(ncol(p) - 1)
                extra = 0
                for r in 1:(ncol(p) - 1)
                    extra += (weights[r, j] * p[r, ncol(p)] + weights[i, r] * (p[nrow(p), r] - p[nrow(p) - 1, r]))
                end
                term1_num += ((p[i, j] * (1 - pdk) ^ 2 * ((1 - pdk) * weights[i, j] - (1 + kd) * extra) ^ 2))
            end
        end
        v1 = term1_num / term1_denom
        v2 = pdk * (1 - pdk) - (1 - pdk) * (1 + kd) * (pdk + ((1 - pdk) * (1 + kd)) / 4)
        v = (v1 + v2) / x1[nrow(x1), ncol(x1)]

        bi_est = (1 + pdk + (1 - pdk) * kd) / 2
        bi_se = sqrt(v)
    end

    alpha = alternative === :two ? (1 - conf) / 2 : 1 - conf
    ci = -quantile(Normal(0.0, 1.0), alpha)

    if alternative === :two
        bi_lcl = bi_est - (ci * bi_se)
        bi_ucl = bi_est + (ci * bi_se)
    elseif alternative === :less
        bi_lcl = 0
        bi_ucl = bi_est + (ci * bi_se)
    elseif alternative === :greater
        bi_lcl = bi_est - (ci * bi_se)
        bi_ucl = 1
    end

    return (bi_est=round(bi_est, digits=3), bi_se=round(bi_se, digits=3), bi_lcl=round(bi_lcl, digits=3), bi_ucl=round(bi_ucl, digits=3))

end

"""
    bibang(x; <keyword arguments>)

Compute Bang's Blinding Indices.

# Arguments

- `x::Matrix{Int64}`: 2×2 or 3×2 integer matrix of cross-tabulated counts
- `weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1]`: use default 1996 James weights (correct guesses are assigned a weight of 0, incorrect guesses are assigned a weight of 0.5, don't know guesses are assigned a weight of 1) unless alternative weights are specified
- `conf::Float64=0.95`: confidence level for the returned confidence intervals
- `alternative::Symbol=:two`: type of returned confidence interval (two-sided, `:two`) or one-sided (`:less` or `:greater`)

# Returns

Named tuple containing:
- `bi_est::Float64`: estimated values for group 1 and 2
- `bi_se::Float64`: standard errors for group 1 and 2
- `bi_lcl::Float64`: lower confidence interval bounds for group 1 and 2
- `bi_ucl::Float64`: upper confidence interval bounds for group 1 and 2

# Notes

The format of the cross-tabulated counts and weights is:

             Treatment  Placebo
Treatment    xxx        xxx 
Placebo      xxx        xxx 
Don't Know   xxx        xxx

If 2×2 cross-tabulated counts matrix is provided, the last row is filled with zeros.

# Sources

1. Bang H, Ni L, Davis CE. Assessment of blinding in clinical trials. Control Clin Trials. 2004 Apr;25(2):143-56. doi: 10.1016/j.cct.2003.10.016. PMID: 15020033.
"""
function bibang(x::Matrix{Int64}; weights::Matrix{Float64}=[0 0.5; 0.5 0; 1 1], conf::Float64=0.95, alternative::Symbol=:two)

    @assert alternative in [:two, :less, :greater] "alternative must be :two, :less, or :greater"

    bi_lcl = zeros(2)
    bi_ucl = zeros(2)

    @assert size(weights) == (3, 2) "weights must be a 3 row by 2 column integer matrix specifying alternative James weights for each cell in x"

    @assert conf > 0.0 && conf < 1.0 "conf must be > 0 and < 1"

    x2 = addmargins(x')

    bi_est = zeros(2)
    bi_se = zeros(2)
    for i in 1:(nrow(x2) - 1)
        bi_est[i] = (2 * (x2[i, i] / (x2[i, 1] + x2[i, 2])) - 1) * 
                    ((x2[i, 1] + x2[i, 2]) / (x2[i, 1] + x2[i, 2] + x2[i, 3]))
        bi_se[i] = sqrt(((x2[i, 1] / x2[i, ncol(x2)]) * (1 - (x2[i, 1] / x2[i, ncol(x2)])) +
                     (x2[i, 2] / x2[i, ncol(x2)]) * (1 - (x2[i, 2] / x2[i, ncol(x2)])) +
                     2 * (x2[i, 1] / x2[i, ncol(x2)]) * (x2[i, 2] / x2[i, ncol(x2)])) / x2[i, ncol(x2)])
    end

    alpha = alternative === :two ? (1 - conf) / 2 : 1 - conf
    ci = -quantile(Normal(0.0, 1.0), alpha)

    if alternative === :two
        bi_lcl = bi_est - (ci * bi_se)
        bi_ucl = bi_est + (ci * bi_se)
    elseif alternative === :less
        bi_lcl = [-1, -1]
        bi_ucl = bi_est + (ci * bi_se)
    elseif alternative === :greater
        bi_lcl = bi_est - (ci * bi_se)
        bi_ucl = [1, 1]
    end

    return (bi_est=round.(bi_est, digits=3), bi_se=round.(bi_se, digits=3), bi_lcl=round.(bi_lcl, digits=3), bi_ucl=round.(bi_ucl, digits=3))

end

end # BlindingIndex