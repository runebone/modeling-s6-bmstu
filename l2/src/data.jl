module Data

export table, R, c, ξ, η1, η2

const table = Dict(
    :n => 1:9,
    :t => 2e3:1e3:1e4,
    :ktv1 => [
        8.200E-03,
        2.768E-02,
        6.560E-02,
        1.281E-01,
        2.214E-01,
        3.516E-01,
        5.248E-01,
        7.472E-01,
        1.025E+00,
    ],
    :ktv2 => [
        1.600E+00,
        5.400E+00,
        1.280E+01,
        2.500E+01,
        4.320E+01,
        6.860E+01,
        1.024E+02,
        1.458E+02,
        2.000E+02,
    ]
)

# const R = 3.5e-1
const R = 3.5e-3
# const c = 3e10
const c = 299_792_458

const ξ = log.(table[:t])
const η1 = log.(table[:ktv1])
const η2 = log.(table[:ktv2])

end
