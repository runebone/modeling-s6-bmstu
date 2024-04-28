module LinearRegression

export m, c

using GLM, DataFrames

include("data.jl")
using .Data: table, ξ, η1, η2

# k(T) = e^(m * log(T) + c)
# log(k(T)) = m * log(T) + c

c = AbstractFloat[0.0, 0.0]
m = AbstractFloat[0.0, 0.0]

df1 = DataFrame(logT = ξ, logkT = η1)
model1 = lm(@formula(logkT ~ logT), df1)
c[1], m[1] = coef(model1)

df2 = DataFrame(logT = ξ, logkT = η2)
model2 = lm(@formula(logkT ~ logT), df2)
c[2], m[2] = coef(model2)

end
