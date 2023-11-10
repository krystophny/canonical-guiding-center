#%%
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt

#from field_uniform import *
from field_bottle import *

def h(x):
    return B(x) / Bmod(x)

def dh(x):
    return dB(x) / Bmod(x) - B(x) * dBmod(x) / Bmod(x)**2

def omc(x):
    return qe * Bmod(x) / (m * c)  # Cyclotron frequency

def domc(x):
    return qe * dBmod(x) / (m * c)

# z = (X1, X2, X3, vpar, phi, vperp)

def rhovec(z):
    return np.sin(z[4]) * e1(z[0:3]) - np.cos(z[4]) * e2(z[0:3])

def rholen(z):
    return z[5] / omc(z[0:3])

def rho(z):
    return rholen(z) * rhovec(z)

def drhovec(z):
    drhovec = np.zeros((6, 3))
    drhovec[0:3, :] = np.sin(z[4]) * de1(z[0:3]) - np.cos(z[4]) * de2(z[0:3])
    drhovec[4, :] = np.cos(z[4]) * e1(z[0:3]) + np.sin(z[4]) * e2(z[0:3])
    return drhovec

def drholen(z):
    drholen = np.zeros(6)
    drholen[0:3] = - z[5] * omc(z[0:3])**(-2) * domc(z[0:3])
    drholen[5] = 1.0 / omc(z[0:3])
    return drholen

def drho(z):
    return np.outer(drholen(z), rhovec(z)) + rholen(z) * drhovec(z)

def vperpvec(z):
    return np.cos(z[4]) * e1(z[0:3]) + np.sin(z[4]) * e2(z[0:3])

def vperplen(z):
    return z[5]

def vperp(z):
    return vperplen(z)*vperpvec(z)

def dvperpvec(z):
    dvperpvec = np.zeros((6, 3))
    dvperpvec[0:3, :] = np.cos(z[4]) * de1(z[0:3]) + np.sin(z[4]) * de2(z[0:3])
    dvperpvec[4, :] = -np.sin(z[4]) * e1(z[0:3]) + np.cos(z[4]) * e2(z[0:3])

    return dvperpvec

def dvperplen(z):
    dvperplen = np.zeros(6)
    dvperplen[5] = 1.0
    return dvperplen

def dvperp(z):
    return np.outer(dvperplen(z), vperpvec(z)) + vperplen(z) * dvperpvec(z)

def q(z):
    return z[0:3] + rho(z)

def p(z):
    return m * z[3] * h(z[0:3]) + m * vperp(z) + qe / c * A(q(z))

def dq(z):
    dq = np.zeros((6, 3))
    dq[0:3, :] = np.identity(3)
    dq += drho(z)
    return dq

def dp(z):
    dp = np.zeros((6, 3))
    dAval = dA(z[0:3])
    drhoval = drho(z)

    dp[0:3, :] = m * z[3] * dh(z[0:3]) + dAval
    dp[3,:] = m*h(z[0:3])

    dp += m * dvperp(z)

    # Inner derviative
    dp += np.dot(drhoval, dAval)

    return dp

def dzcan(z):
    return np.concatenate((dq(z), dp(z)), axis=1)

def dz(z):
    return np.linalg.inv(dzcan(z))

def H(z):
    return 0.5 * m * (z[3]**2 + z[5]**2)

def dH(z):
    dH = np.zeros(6)
    dH[3] = m*z[3]
    dH[5] = m*z[5]
    return dH

# Transformation of canonical particle coordinates to guiding-center coordinates
def Fcan(z, qval, pval):
    Fcan = np.zeros(6)
    Fcan[0:3] = qval - q(z)
    Fcan[3:6] = pval - p(z)
    return Fcan

def z(q,p,zold):
    sol = root(Fcan, zold, args=(q,p), tol=1e-13)
    return sol.x

# Midpoint rule in particle coordinates including implicit
# guiding-center transformation
def F(z, qold, pold):
    F = np.zeros(6)
    qmid = q(z)
    pmid = p(z)
    dHval = dH(z)
    dzval = dz(z)
    F[0:3] = qmid - (qold + 0.5*dt * np.dot(dzval[3:,:], dHval))
    F[3:6] = pmid - (pold - 0.5*dt * np.dot(dzval[:3,:], dHval))
    return F



nt = 1000
zs = np.zeros((nt, 6))
zs_euler = np.zeros((nt, 6))

z0 = np.array([0.2, -0.1, 0.8, -0.01, 0.0, 0.1])
zs[0, :] = z0
zs_euler[0, :] = z0

tau = 2*np.pi/omc(z0[:3])
dt = 1.234567e-2 * tau

qs = np.zeros((nt, 3))
ps = np.zeros((nt, 3))

qs[0, :] = q(z0)
ps[0, :] = p(z0)

qs_euler = np.zeros((nt, 3))
ps_euler = np.zeros((nt, 3))

qs_euler[0, :] = q(z0)
ps_euler[0, :] = p(z0)

zguess = np.zeros(6)

for kt in range(1, nt):
    zold = zs[kt - 1, :]
    qold = qs[kt - 1, :]
    pold = ps[kt - 1, :]

    # Guess for midpoint
    zguess[:] = zold
    zguess[:3] += 0.5 * dt * zold[3]*h(zold[0:3])
    zguess[4] += 0.5 * dt * omc(zold[0:3])

    sol = root(lambda z: F(z, qold, pold), zguess, tol=1e-13)
    zmid = sol.x

    dzval = dz(zmid)
    dHval = dH(zmid)

    qnew = qold + dt * np.dot(dzval[3:,:], dHval)
    pnew = pold - dt * np.dot(dzval[:3,:], dHval)

    zguess[:] = zmid
    zguess[:3] += 0.5 * dt * zmid[3]*h(zmid[0:3])
    zguess[4] += 0.5 * dt * omc(zmid[0:3])

    znew = z(qnew, pnew, zguess)

    zs[kt, :] = znew
    qs[kt, :] = qnew
    ps[kt, :] = pnew

    dzval_euler = dz(zs_euler[kt - 1, :])
    qs_euler[kt, :] = qs_euler[kt - 1, :] + dt * np.dot(
        dzval_euler[3:,:], dH(zs_euler[kt - 1, :])
    )
    ps_euler[kt, :] = ps_euler[kt - 1, :] - dt * np.dot(
        dzval_euler[:3,:], dH(zs_euler[kt - 1, :])
    )
    zs_euler[kt, :] = z(qs_euler[kt, :], ps_euler[kt, :], zs_euler[kt - 1, :])

## Plot

plt.plot(qs_euler[:, 0], qs_euler[:, 1], 'r,', label="Particle (Euler)")
plt.plot(qs[:, 0], qs[:, 1], 'g,', label="Particle (SIMPLE)")
plt.plot(zs[:, 0], zs[:, 1], 'b,', label="Guiding-Center (SIMPLE)")
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.legend()
plt.show()

plt.plot(np.sqrt(qs_euler[:, 0]**2 + qs_euler[:, 1]**2), qs_euler[:,2], 'r,', label="Particle (Euler)")
plt.plot(np.sqrt(qs[:, 0]**2 + qs[:, 1]**2), qs[:,2], 'g,', label="Particle (SIMPLE)")
plt.plot(np.sqrt(zs[:, 0]**2 + zs[:, 1]**2), qs[:,2], 'b,', label="Guiding-Center (SIMPLE)")
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.legend()
plt.show()

tnorm = 2*np.pi*np.linspace(0.0,nt*dt,nt)/omc(zs[0,0:3])
plt.plot(tnorm, zs[:,4],'.-')
plt.xlabel('t/(omc/2*pi)')
plt.show()

plt.plot(tnorm,[H(zs[kt,:]) for kt in range(nt)])
plt.plot(tnorm,[H(zs_euler[kt,:]) for kt in range(nt)])
plt.ylim(0,0.02)
plt.show()

plt.plot(tnorm,zs[:,5],'.-')

# %%
