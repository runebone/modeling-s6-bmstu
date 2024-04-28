using GLMakie

include("interp.jl")
using .Interp: k1, k2

x = log.(LinRange(2000, 10000, 100))
y1 = k1.(x)
y2 = k2.(x)

fig = Figure()

ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])

lines!(ax1, x, y1, color = :blue)
lines!(ax1, x, y2, color = :red)

include("linregr.jl")
using .Linregr: m, c

lines!(ax2, x, (x -> m[1] * x + c[1]).(x), color = :blue)
lines!(ax2, x, (x -> m[2] * x + c[2]).(x), color = :red)

function T(z, Tw=2e3, T₀=1e4, p=4)
    @assert 0 ≤ z ≤ 1
    return (Tw - T₀) * z^p + T₀
end

function k(z, v = 1)
    @assert 0 ≤ z ≤ 1
    @assert v == 1 || v == 2
    return exp(m[v] * log(T(z)) + c[v])
end

ax3 = Axis(fig[1:2, 2])

z = LinRange(0, 1, 100)
lines!(ax3, z, (z -> k(z, 1)).(z), color = :blue)
lines!(ax3, z, (z -> k(z, 2)).(z), color = :red)

fig
