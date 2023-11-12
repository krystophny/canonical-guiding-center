using Plots
using NLsolve

include("field_uniform.jl")

const m = 1.0
const qe = 1.0
const c = 1.0
const Jperp = 1.0

# Non-canonical coordinates: z = [x1, x2, x3, vpar]
# Canonical coordinates: q = [x2, x3], p from guiding-center Lagrangian

omc(x) = qe * Bmod(x) / (m * c)
domc(x) = qe * dBmod(x) / (m * c)

h(x) = B(x) / Bmod(x)
dh(x) = dB(x) / Bmod(x) -  dBmod(x)*B(x)' / Bmod(x)^2

q(z) = z[2:3]
dq(z) = [0.0 0.0;
         1.0 0.0;
         0.0 1.0;
         0.0 0.0]

p(z) = m .* z[4] .* h(z[1:3])[2:3] .+ qe / c .* A(z[1:3])[2:3]
function dp(z)
    dp = zeros(4, 2)
    dp[1:3, :] = m*z[4]*dh(z[1:3])[1:3,2:3] + qe/c*dA(z[1:3])[1:3,2:3]
    dp[4, :] = m*h(z[1:3])[2:3]
    return dp
end

dz(z) = inv([dq(z) dp(z)])

H(z) = 0.5 * m * z[4]^2 + Jperp*omc(z[1:3])
dH(z) = [Jperp * domc(z[1:3]); m*z[4]]

# Transformation of canonical particle coordinates to guiding-center coordinates
function Fcan(z, qval, pval)
    Fcan = zeros(4)
    Fcan[1:2] .= qval .- q(z)
    Fcan[3:4] .= pval .- p(z)
    return Fcan
end

function z(q, p, zold)
    sol = nlsolve(x -> Fcan(x, q, p), zold, ftol=1e-13)
    return sol.zero
end

"""
F(z, qold, pold, dt)

Compute the residuals of the midpoint rule for root finding. The output
becomes zero for `z = zmid`, which transforms to the canonical coordinates `qmid` and `pmid` realized at the midpoint `t + dt/2`.

# Arguments
- `z::Vector{Float64}`: Non-canonical coordinate tuple to test.
- `qold::Vector{Float64}`: The position coordinates at the previous time step t.
- `pold::Vector{Float64}`: The momentum coordinates at the previous time step t.
- `dt::Float64`: The time step size.

# Returns
- `Float64`: The residual of the midpoint rule.
"""
function F(z, qold, pold, dt)
    F = zeros(4)
    qmid = q(z)
    pmid = p(z)
    dHval = dH(z)
    dzval = dz(z)
    F[1:2] .= qmid .- (qold .+ 0.5 * dt * (dzval[3:end, :] * dHval))
    F[3:4] .= pmid .- (pold .- 0.5 * dt * (dzval[1:2, :] * dHval))
    return F
end

# Initialize orbit

z0 = [0.2, 0.1, 0.02, 0.2]
qold = q(z0)
pold = p(z0)
dt = 0.1
nt = 100
zs = zeros(4, nt)
zs[:, 1] = z0

for kt = 1:nt-1
    # Find zmid that transforms to midpoint of canonical coordinates
    sol = nlsolve(x -> F(x, qold, pold, dt), zs[:, kt])
    zmid = sol.zero

    # Update canonical coordinates for next iteration
    dzval = dz(zmid)
    dHval = dH(zmid)
    qold += dt * dzval[3:4,:] * dHval
    pold -= dt * dzval[1:2,:] * dHval

    # Store values of non-canonical coordinates for time-step
    zs[:, kt+1] = z(qold, pold, zmid)
end
