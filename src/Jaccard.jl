using Statistics, LsqFit, Distributions, LinearAlgebra

"""
Return length(intersection(x,y)), length(union(x,y)) for Sets `x` and `y`
without creating (allocating) the union and intersection sets.
"""
function jaccardIdent(x::Set{String}, y::Set{String})::Float64
    _intersect = count(i -> i ∈ y, x)
    _union = length(x) + length(y) - _intersect
    return _intersect / _union
end

"""
Return a generator traversing an upper triangle. Returns (i,j) indices,
intended to be used as:
```
for (i,j) in uppertriangle(mat)
        ...
    end
```
"""
uppertriangle(A) = ((i, j) for j in axes(A, 2) for i in 1:j-1)

# NOTE: pre-allocate a Set() with an estimated state
#s = Set{String}()
# sizehint!(s, 1000)

#=
MAKE SURE TO DECONSTRUCT THIS INTO A NAMEDTUPLE WHEN IT IS USED:
e.g. (; xs, ys, qbottom, qtop, z) = jaccardScores(data)
 =#
"""
Using the Jaccard identity matrix (intersect/union) as input,
calculate the Z-scores and confidence intervals. Returns a `NamedTuple`
(:xs, :ys, :qbottom, :qtop, :zscore), each a `Vector{Int64}`. The output is
in a Tables-compliant format, suitable for downstream plotting.
"""
function jaccardScores(mat::Matrix{Float64}; alpha::Float64=0.95, zscore::Union{Int64, Float64})::NamedTuple
    n = size(mat, 1)
    N = n * (n - 1) ÷ 2

    ys = Vector{Float64}(undef, N)
    xs = Vector{Float64}(undef, N)
    z_scores = Vector{Float64}(undef, N)

    @inbounds for diag in 1:(n-1)
        diag_len = n - diag
        μ = 0.0
        for j in (diag+1):n
            μ += mat[j-diag, j]
        end
        μ /= diag_len

        σ² = 0.0
        for j in (diag+1):n
            σ² += (mat[j-diag, j] - μ)^2
        end
        σ = sqrt(σ² / (diag_len - 1))

        for j in (diag+1):n
            k = (diag - 1) * n - diag * (diag - 1) ÷ 2 + (j - diag)
            val = Float64(mat[j-diag, j])
            ys[k] = val
            xs[k] = Float64(diag)
            z_scores[k] = σ == 0 ? 0.0 : (val - μ) / σ
        end
    end
    return (xs=xs, ys=ys, z_scores=z_scores)
end

"""
Fit prediction bands for plotting
"""
function fit_prediction_bands(xs::Vector{Float64}, ys::Vector{Float64}, alpha::Float64)::NamedTuple
    model(x, p) = exp.(p[1] .+ p[2] .* exp.(-x .* p[3]))

    fit = curve_fit(model, xs, ys, [1.0, 1.0, 1.0])
    fitted = model(xs, fit.param)

    n = length(ys)
    nparams = length(fit.param)
    σ_resid = std(ys .- fitted)
    t_crit = quantile(TDist(n - nparams), 1 - (1 - alpha) / 2)

    J = jacobian(fit)
    leverage = vec(sum((J * inv(J'J)) .* J, dims=2))

    Qbottom = Vector{Float64}(undef, n)
    Qtop = Vector{Float64}(undef, n)
    abs_z = Vector{Float64}(undef, n)  # also got dropped earlier

    @inbounds for k in 1:n
        se = σ_resid * sqrt(1.0 + leverage[k])
        Qbottom[k] = fitted[k] - t_crit * se
        Qtop[k] = fitted[k] + t_crit * se
    end

    return (fitted=fitted, Qbottom=Qbottom, Qtop=Qtop)
end

#= intended usage
(; xs, ys, z_scores) = extract_and_zscore(data)
(; fitted, Qbottom, Qtop) = fit_prediction_bands(xs, ys, alpha)
abs_z = abs.(z_scores)

@vlplot(
    data = (x=xs, y=ys, z=z_scores, az=abs_z, f=fitted, qb=Qbottom, qt=Qtop),
    mark = :point,
    x = {"x:q", title = "Diagonal distance"},
    y = {"y:q", title = "Value"},
    color = {"z:q", scale = {scheme = "redblue"}}
)
=#
