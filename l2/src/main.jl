include("linear_regression.jl") # Get k coefficients
include("data.jl")

function T(z, Tw=2e3, T₀=1e4, p=4)
    @assert 0 ≤ z ≤ 1 "$z"
    return (Tw - T₀) * z^p + T₀
end

function uₚ(z)
    @assert 0 ≤ z ≤ 1 "$z"
    🌭 = 3.084 * 1e-4
    🦖 = exp(4.799 * 1e4 / T(z)) - 1
    return 🌭 / 🦖
end

function k(🦓, 🎻)
    @assert 0 ≤ 🦓 ≤ 1 "$🦓"
    @assert 🎻 == 1 || 🎻 == 2 # Variant 1 or 2
    🧙 = LinearRegression.m
    🐈 = LinearRegression.c
    return exp(🧙[🎻] * log(T(🦓)) + 🐈[🎻])
end

# function rk4(f, x, uₙ, h)
#     # u'(x) = f(x, u)
#     k1 = f(x, uₙ)
#     k2 = f(x + h / 2, uₙ + h / 2 * k1)
#     k3 = f(x + h / 2, uₙ + h / 2 * k2)
#     k4 = f(x + h, uₙ + h * k3)
#     uₙ₊₁ = uₙ + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
#     return uₙ₊₁
# end

function rk4(f, ϕ, x, uₙ, vₙ, h)
    # u'(x) = f(x, u, v)
    # v'(x) = ϕ(x, u, v)
    k1 = f(x, uₙ, vₙ)
    q1 = ϕ(x, uₙ, vₙ)
    k2 = f(x + h / 2, uₙ + h / 2 * k1, vₙ + h / 2 * q1)
    q2 = ϕ(x + h / 2, uₙ + h / 2 * k1, vₙ + h / 2 * q1)
    k3 = f(x + h / 2, uₙ + h / 2 * k2, vₙ + h / 2 * q2)
    q3 = ϕ(x + h / 2, uₙ + h / 2 * k2, vₙ + h / 2 * q2)
    k4 = f(x + h, uₙ + h * k3, vₙ + h * q3)
    q4 = ϕ(x + h, uₙ + h * k3, vₙ + h * q3)
    uₙ₊₁ = uₙ + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    vₙ₊₁ = vₙ + h / 6 * (q1 + 2 * q2 + 2 * q3 + q4)
    return uₙ₊₁, vₙ₊₁
end

function rk4_step(z, uₙ, Fₙ, h, v)
    @assert 0 ≤ z ≤ 1 "$z"
    @assert v == 1 || v == 2 # Variant 1 or 2
    R = Data.R
    c = Data.c
    f(z, _, F) = -3 * R * k(z, v) * F / c
    ϕ(z, u, F) = R * c * k(z, v) * (uₚ(z) - u) - F / z
    if z == 0
        ϕ(z, u, _) = R * c * k(z, v) * (uₚ(z) - u)
    end
    return rk4(f, ϕ, z, uₙ, Fₙ, h)
end

function shoot(v)
    @assert v == 1 || v == 2 # Variant 1 or 2
    ξ₀, ξ₁ = 0, 1
    ξ = (ξ₀ + ξ₁) / 2
    ϵ = 1e-6
    δ = Inf
    h = 1e-4
    c = Data.c
    zrange = 0:h:1-h
    while abs(δ) > ϵ
        u = ξ * uₚ(0)
        F = 0
        for z in zrange
            u, F = rk4_step(z, u, F, h, v)
        end
        δ = F - 0.393 * c * u
        if δ > 0
            ξ₀ = ξ
        else
            ξ₁ = ξ
        end
        ξ = (ξ₀ + ξ₁) / 2
    end
    return ξ
end

using GLMakie

function main()
    ξ1 = shoot(1)
    ξ2 = shoot(2)

    println("$ξ1 $ξ2")

    h = 1e-2
    zrange = 0:h:1-h

    u1 = ξ1 * uₚ(0)
    F1 = 0.0
    u1v = [u1]
    F1v = [F1]
    for z in zrange
        u1, F1 = rk4_step(z, u1, F1, h, 1)
        append!(u1v, u1)
        append!(F1v, F1)
    end
    # println(u1v)
    # println(F1v)

    u2 = ξ2 * uₚ(0)
    F2 = 0.0
    u2v = [u2]
    F2v = [F2]
    for z in zrange
        u2, F2 = rk4_step(z, u2, F2, h, 2)
        append!(u2v, u2)
        append!(F2v, F2)
    end
    # println(u2v)
    # println(F2v)

    fig = Figure()

    ax11 = Axis(fig[1, 1], title="Вариант 1 (F)")
    ax12 = Axis(fig[1, 2], title="Вариант 2 (F)")
    ax21 = Axis(fig[2, 1], title="Вариант 1 (u)")
    ax22 = Axis(fig[2, 2], title="Вариант 2 (u)")

    zrange = 0:h:1

    lines!(ax11, zrange, F1v, color=:blue, label="F")
    lines!(ax12, zrange, F2v, color=:blue, label="F")
    lines!(ax21, zrange, u1v, color=:red, label="u")
    lines!(ax22, zrange, u2v, color=:red, label="u")
    lines!(ax21, zrange, uₚ.(zrange), color=:green, label="uₚ")
    lines!(ax22, zrange, uₚ.(zrange), color=:green, label="uₚ")

    axislegend(ax11, position=:lt)
    axislegend(ax12, position=:lt)
    axislegend(ax21)
    axislegend(ax22)

    fig
end
