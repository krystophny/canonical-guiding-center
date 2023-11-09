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

e1(x) = [1,0,0]
e2(x) = [0,1,0]

rholen(z) = sqrt(2*z[6]/(m*omc(z[1:3])))
rhovec(z) = sin(z[4])*e1(z[1:3]) + cos(z[4])*e2(z[1:3])
rho(z) = rholen(z).*rhovec(z)

q(z) = z[1:3] + rho(z)
p(z) = m*z[4]*h(z[1:3]) + qe/c*A(q(z))

z = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

println(q(z))
println(p(z))
