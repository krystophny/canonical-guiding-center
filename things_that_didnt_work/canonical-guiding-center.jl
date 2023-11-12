using Plots

# z(1:3) = x
# z(4) = vpar
# z(5) = phi
# z(6) = Jperp

B0 = 1.0
m = 1.0
qe = 1.0
c = 1.0

A(x) = B0*[-0.5*x[2], 0.5*x[1], 0]
B(x) = B0*[0, 0, 1]
Bmod(x) = sqrt(B(x)'*B(x))
h(x) = B(x)/Bmod(x)
omc(x) = qe*Bmod(x)/(m*c)  # Cyclotron frequency

e1(x) = [1.0, 0.0, 0.0]
e2(x) = [0.0, 1.0, 0.0]

rhovec(z) = sin(z[5])*e1(z[1:3]) - cos(z[5])*e2(z[1:3])
rholen(z) = sqrt(2*z[6]/(m*omc(z[1:3])))
rho(z) = rholen(z)*rhovec(z)

vperpvec(z) = cos(z[5])*e1(z[1:3]) + sin(z[5])*e2(z[1:3])

q(z) = z[1:3] .+ rho(z)
p(z) = m*z[4]*h(z[1:3]) + rholen(z)*omc(z[1:3])*vperpvec(z) + qe/c*A(q(z))

phis = 0:2π/20:2π
zs = zeros(length(phis), 6)
zs[:,5] = phis
zs[:,6] .= 1.0

qvec = reduce(hcat, q.(eachrow(zvec)))

## Plot

plot(qvec[1,:], qvec[2,:], aspect_ratio=:equal, legend=false)

## Derivatives

drho(z) = drholen(z)*rhovec(z) + rholen(z)*drhovec(z)

"""
Vector fields are always indexed as dv_ij = d/dx_i v_j (i=1..3, j=1..3)
Phase space vector fields like dv_ij = d/dz_i v_j (i=1..6, j=1..3)
"""
function de1(x)
    zeros(3,3)
end

function de2(x)
    zeros(3,3)
end

"""
\d\hat{\vrho}(\xset,\phi)&=\left(\sin\phi\frac{\partial\mathbf{e}_{1}(\xset)}{\partial x^{k}}-\cos\phi\frac{\partial\mathbf{e}_{2}(\xset)}{\partial x^{k}}\right)\d x^{k}\\&+(\cos\phi\mathbf{e}_{1}(\xset)+\sin\phi\mathbf{e}_{2}(\xset))\d\phi
"""
function drhovec(z)
    drhovec = zeros(6,3)
    drhovec[1:3, 1:3] = sin(z[5])*de1(z[1:3]) - cos(z[5])*de2(z[1:3])
    rhovec[5, 1:3] = cos(z[5])*e1(z[1:3]) + sin(z[5])*e2(z[1:3])
end

"""
\d\rho_{c}(\xset,J_{\perp})=\d\sqrt{\frac{2J_{\perp}}{m\omega_{c}(\xset)}}
=-\frac{\rho_{c}}{2\omega_{c}(\xset)}\frac{\partial\omega_{c}(\xset)}{\partial x^{k}}\d x^{k}
+ \frac{1}{2J_{\perp}}\d J_{\perp}
"""
function drholen(z)
    drholen = zeros(6)
    drholen[1:3] = -0.5*rholen(z)/omc(z[1:3])*domc(z[1:3])
    drholen[6] = 0.5/z[6]
end

function dq(z)
    dq = zeros(6,3)
    dq[1:3, 1:3] = eye(3)
end

function F(z, zold)
   F = zeros(6)
   qold = q(zold)
   pold = p(zold)
   F[1:3] = z[1:3] + rho(z) - (qdot(0.5*(zold + znew)) + )
   )
   F[4:6] = p(z) - pnew
end
