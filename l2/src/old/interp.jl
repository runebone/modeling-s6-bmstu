module Interp

export k1, k2

using Interpolations

include("data.jl")
using .Data: table, ξ, η1, η2

# itp_type = BSpline ∘ Cubic ∘ Line ∘ OnGrid
# itp1 = scale(interpolate(η1, itp_type()), ξ)
# itp2 = scale(interpolate(η2, itp_type()), ξ)

itp1 = LinearInterpolation(ξ, η1)
itp2 = LinearInterpolation(ξ, η2)

k1(T) = itp1(T)
k2(T) = itp2(T)

end
