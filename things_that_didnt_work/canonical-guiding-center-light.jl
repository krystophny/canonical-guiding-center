using Plots
using NLsolve

# z(1:3) = x
# z(4) = vpar
# z(5) = phi
# z(6) = vperp

B0 = 1.0
m = 1.0
qe = 1.0
c = 1.0

Bmod(x) = B0*(1.0 + x[3]^2)
A(x) = Bmod(x)*[-0.5*x[2], 0.5*x[1], 0]

"""
d/dx_i A_j
"""
function dA(x)
    dA = zeros(3,3)
    dA[1,2] = 0.5*(1.0 + x[3]^2)
    dA[2,1] = -0.5*(1.0 + x[3]^2)
    dA[3,1] = -x[2]*x[3]
    dA[3,2] = x[1]*x[3]
end

B(x) = B0*[-x[1]*x[3], -x[2]*x[3], 1+x[3]^2]
h(x) = B(x)/Bmod(x)
omc(x) = qe*Bmod(x)/(m*c)  # Cyclotron frequency

e1(x) = [1.0, 0.0, 0.0]
e2(x) = [0.0, 1.0, 0.0]

rhovec(z) = sin(z[5])*e1(z[1:3]) - cos(z[5])*e2(z[1:3])
rholen(z) = z[6]/omc(z[1:3])
rho(z) = rholen(z)*rhovec(z)

vperpvec(z) = cos(z[5])*e1(z[1:3]) + sin(z[5])*e2(z[1:3])

q(z) = z[1:3] .+ rho(z)
p(z) = m*z[4]*h(z[1:3]) + m*z[6]*vperpvec(z[1:3]) + qe/c*A(q(z))

H(q,p) = 0.5*((p-qe/c*A(q))'*(p-qe/c*A(q)))/m

dHdq(q,p) = -qe/c*dA(q)*(p-qe/c*A(q))/m
dHdp(q,p) = (p-qe/c*A(q))/m

z0 = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
q0 = q(z0)
p0 = p(z0)

dt = 0.01/omc(z0[1:3])

function F(z, qold, pold)
   F = zeros(6)
   qnew = q(z)
   pnew = p(z)
   qmid = 0.5*(qold + qnew)
   pmid = 0.5*(pold + pnew)
   F[1:3] = qnew - (qold + dt*dHdp(qmid, pmid))
   F[4:6] = pnew - (pold - dt*dHdq(qmid, pmid))
   F
end

nt = 10
zs = zeros(nt,6)
qs = zeros(nt,3)
ps = zeros(nt,3)
zs[1,:] = [0.0, 0.1, 0.1, 0.0, 0.0, 0.1]
qs[1,:] .= q(z0)
ps[1,:] .= p(z0)

zguess = zeros(6)

for kt = 1:nt-1
    zold = zs[kt,:]
    qold = q(zold)
    pold = p(zold)

    zguess .= zold
    zguess[5] += dt*2.0*pi/omc(zold[1:3])

    sol = nlsolve(z -> F(z, qold, pold), zguess, autodiff = :forward)
    zs[kt+1,:] .= sol.zero
    qs[kt+1,:] .= q(zs[kt+1,:])
    ps[kt+1,:] .= p(zs[kt+1,:])
end

## Plot
scatter(zs[:,1], zs[:,2], label="Guiding center")
plot!(qs[:,1], qs[:,2], label="Particle")
xlims!(-1.5,1.5)
ylims!(-1.5,1.5)
