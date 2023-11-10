#%%
import numpy as np
from scipy.optimize import root

#from field_uniform import *
from field_bottle import *

def h(x):
    return B(x) / Bmod(x)

def omc(x):
    return qe * Bmod(x) / (m * c)  # Cyclotron frequency

# z = (X1, X2, X3, vpar, phi, vperp)

def rhovec(z):
    return np.sin(z[4]) * e1(z[0:3]) - np.cos(z[4]) * e2(z[0:3])

def rholen(z):
    return z[5] / omc(z[0:3])

def vperpvec(z):
    return np.cos(z[4]) * e1(z[0:3]) + np.sin(z[4]) * e2(z[0:3])

def rho(z):
    return rholen(z) * rhovec(z)

def q(z):
    return z[0:3] + rho(z)

def p(z):
    return m * z[3] * h(z[0:3]) + m * z[5] * vperpvec(z) + qe / c * A(q(z))

def H(q, p):
    return 0.5 * np.dot(p - qe / c * A(q), p - qe / c * A(q)) / m

def dHdq(q, p):
    return -qe / c * np.dot(dA(q), (p - qe / c * A(q))) / m

def dHdp(q, p):
    return (p - qe / c * A(q)) / m

z0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
q0 = q(z0)
p0 = p(z0)

tau = 2*np.pi/omc(z0[0:3])
dt = 1.234567 * tau

# Midpoint rule in particle coordinates including implicit
# guiding-center transformation
def F(z, qold, pold):
    F = np.zeros(6)
    qnew = q(z)
    pnew = p(z)
    qmid = 0.5 * (qold + qnew)
    pmid = 0.5 * (pold + pnew)
    F[0:3] = qnew - (qold + dt * dHdp(qmid, pmid))
    F[3:6] = pnew - (pold - dt * dHdq(qmid, pmid))
    return F

# Midpoint rule in canonical particle coordinates only for comparison
def Fcan(zcan, zold):
    Fcan = np.zeros(6)
    qnew = zcan[0:3]
    pnew = zcan[3:]
    qmid = 0.5 * (zold[0:3] + qnew)
    pmid = 0.5 * (zold[3:]  + pnew)
    Fcan[0:3] = qnew - (qold + dt * dHdp(qmid, pmid))
    Fcan[3:6] = pnew - (pold - dt * dHdq(qmid, pmid))
    return Fcan

nt = 1000
zs = np.zeros((nt, 6))

qs = np.zeros((nt, 3))
ps = np.zeros((nt, 3))
qs_euler = np.zeros((nt, 3))
ps_euler = np.zeros((nt, 3))

zs[0, :] = np.array([0.2, -0.1, 0.8, 0.0, 0.0, 0.2])

qs[0, :] = q(zs[0, :])
ps[0, :] = p(zs[0, :])

qs_euler[0, :] = qs[0, :]
ps_euler[0, :] = ps[0, :]

zs_can = np.zeros((nt, 6))
zs_can[0, 0:3] = qs[0, :]
zs_can[0, 3:]  = ps[0, :]

zguess = np.zeros(6)

for kt in range(1, nt):
    zold = zs[kt - 1, :]
    qold = q(zold)
    pold = p(zold)

    zguess[:] = zold
    zguess[4] += dt * 2.0 * np.pi / omc(zold[0:3])

    sol = root(lambda z: F(z, qold, pold), zguess, tol=1e-13)

    zs[kt, :] = sol.x
    qs[kt, :] = q(zs[kt, :])
    ps[kt, :] = p(zs[kt, :])

    qs_euler[kt, :] = qs[kt - 1, :] + dt * dHdp(qs[kt - 1, :], ps[kt - 1, :])
    ps_euler[kt, :] = ps[kt - 1, :] - dt * dHdq(qs[kt - 1, :], ps[kt - 1, :])

    sol_can = root(lambda z: Fcan(z, zs_can[kt-1,:]), zs_can[kt-1,:], tol=1e-13)
    zs_can[kt,:] = sol_can.x

## Plot
import matplotlib.pyplot as plt

plt.plot(zs[:, 0], zs[:, 1], 'b.', label="Guiding center (SIMPLE)")
plt.plot(qs[:, 0], qs[:, 1], 'g.', label="Particle (SIMPLE)")
plt.plot(qs_euler[:, 0], qs_euler[:, 1], 'r.', label="Particle (Euler)")
plt.plot(zs_can[:, 0], zs_can[:, 1], 'm.', label="Particle (Naive midpoint)")
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.legend()
plt.show()

plt.plot(np.sqrt(zs[:, 0]**2 + zs[:, 1]**2), zs[:,2], 'b.', label="Guiding center (SIMPLE)")
plt.plot(np.sqrt(qs[:, 0]**2 + qs[:, 1]**2), qs[:,2], 'g.', label="Particle (SIMPLE)")
plt.plot(np.sqrt(qs_euler[:, 0]**2 + qs_euler[:, 1]**2), qs_euler[:,2], 'r.', label="Particle (Euler)")
plt.plot(np.sqrt(zs_can[:, 0]**2 + zs_can[:, 1]**2), zs_can[:,2], 'm.', label="Particle (Naive midpoint)")
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.legend()
plt.show()

tnorm = 2*np.pi*np.linspace(0.0,nt*dt,nt)/omc(z0[0:3])
plt.plot(tnorm, zs[:,4],'.-')
plt.xlabel('t/(omc/2*pi)')
plt.show()

plt.plot(tnorm,[H(qs[kt,:], ps[kt,:]) for kt in range(nt)])


# %%
