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

rholen(z) = sqrt(2*z[6]/(m*omc(z[1:3])))
rhovec(z) = sin(z[5])*e1(z[1:3]) + cos(z[5])*e2(z[1:3])
rho(z) = rholen(z)*rhovec(z)

q(z) = z[1:3] .+ rho(z)
p(z) = m*z[4]*h(z[1:3]) + qe/c*A(q(z))

phis = 0:2π/20:2π
zs = zeros(length(phis), 6)
zs[:,5] = phis
zs[:,6] .= 1.0

qvec = reduce(hcat, q.(eachrow(zvec)))

## Plot
using Plots

plot(qvec[1,:], qvec[2,:], aspect_ratio=:equal, legend=false)
