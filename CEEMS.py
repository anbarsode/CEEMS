import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Basic data
m1 = np.random.random()
m2 = 1.0
k = 1.0 # F = k*r**p
p = -2.0 # 1 - 4 * np.random.random()

# Initial conditions
'''# Random, adjusted to look pretty on average
r10 = np.random.random((2))-0.5
r20 = -r10 * m1 / m2
v10 = 5.0 * (np.random.random((2))-0.5)
v20 = -v10 * m1 / m2'''

# Circular orbits
r10 = np.array((2.0,0.0))
r20 = -r10 * m1 / m2
v10 = (2/m1*np.linalg.norm(r10-r20)**p)**0.5 * np.array((0.0,1.0))
v20 = -v10 * m1 / m2

# Numerical integration parameters
dt = 0.01
t = 50.0
Nsteps = int(t/dt)

model = 3 # 1 -> Random equivalent equal mass
          # 2 -> Non-gravitational central force, same universe
          # 3 -> Newtonian gravity, same universe
          
figname = './CEEMS_gravity.png'

#####################################################################################

def F(r1, r2, k):
    return k * np.linalg.norm(r1-r2)**(p-1) * (r1-r2)

def transformation(m1, m2, p, k, model):
    M = m1 + m2
    mu = m1 * m2 / M
    Mc = mu**0.6 * M**0.4 # Chirp mass

    if model == 1: A = np.random.random()
    elif model == 2: A = 2 * mu
    elif model == 3: A = 2**0.2 * Mc
    else: print('Error : Invalid model')
    
    gamma = (2.0 * mu / A)**0.5
    beta1 = -gamma/2.0 - (m1/2.0/A - gamma**2.0/4.0)**0.5
    alpha1 = beta1 + gamma
    beta2 = gamma/2.0 - (m2/2.0/A - gamma**2.0/4.0)**0.5
    alpha2 = beta2 - gamma
    J = alpha1*beta2 - alpha2*beta1
    k_new = k / gamma**(p+1)
        
    return A, alpha1, alpha2, beta1, beta2, gamma, J, k_new

def single_step_RK4(r1i, v1i, r2i, v2i, m1, m2, F, k, dt):
    c10 = v1i * dt
    c20 = v2i * dt
    d10 = -F(r1i, r2i, k) / m1 * dt
    d20 = -d10 * m1 / m2
    
    c11 = (v1i + d10 * 0.5) * dt
    c21 = (v2i + d20 * 0.5) * dt
    d11 = -F(r1i + c11 * 0.5, r2i + c21 * 0.5, k) / m1 * dt
    d21 = -d11 * m1 / m2
    
    c12 = (v1i + d11 * 0.5) * dt
    c22 = (v2i + d21 * 0.5) * dt
    d12 = -F(r1i + c12 * 0.5, r2i + c22 * 0.5, k) / m1 * dt
    d22 = -d12 * m1 / m2
    
    c13 = (v1i + d12 * 0.5) * dt
    c23 = (v2i + d22 * 0.5) * dt
    d13 = -F(r1i + c13 * 0.5, r2i + c23 * 0.5, k) / m1 * dt
    d23 = -d13 * m1 / m2
    
    r1f = r1i + (c10 + (c11 + c12) * 2.0 + c13) * (1/6.0)
    r2f = r2i + (c20 + (c21 + c22) * 2.0 + c23) * (1/6.0)
    v1f = v1i + (d10 + (d11 + d12) * 2.0 + d13) * (1/6.0)
    v2f = v2i + (d20 + (d21 + d22) * 2.0 + d23) * (1/6.0)
    return r1f, v1f, r2f, v2f
    
def simulation(r10, v10, r20, v20, m1, m2, F, k, dt, Nsteps):
    r1 = [r10]
    r2 = [r20]
    
    r1t = r10.copy()
    v1t = v10.copy()
    r2t = r20.copy()
    v2t = v20.copy()
    
    for _ in tqdm(range(Nsteps)):
        r1t, v1t, r2t, v2t = single_step_RK4(r1t, v1t, r2t, v2t, m1, m2, F, k, dt)
        r1.append(r1t)
        r2.append(r2t)
        
    return np.array(r1), np.array(r2)

def make_plots(r1, r2, ra, rb, r1ab, r2ab, skip, Nsteps, model, figtitle, figname):
    maxlim = np.max([r1.max(), r2.max()])
    minlim = np.min([r1.min(), r2.min()])
    pad = (maxlim - minlim) * 0.05
    maxlim += pad
    minlim -= pad

    fig,ax = plt.subplots(2,2,figsize=(8,8))
    plt.suptitle(figtitle)

    ax[0][0].scatter(r1[::100,0], r1[::100,1], s = 10/np.arange(Nsteps+1,0,-100)**0.75, label=r'$r_1$', color = 'r')
    ax[0][0].scatter(r2[::100,0], r2[::100,1], s = 10/np.arange(Nsteps+1,0,-100)**0.75, label=r'$r_2$', color = 'b')
    ax[0][0].set_xlabel('x')
    ax[0][0].set_ylabel('y')
    ax[0][0].set_xlim([minlim,maxlim])
    ax[0][0].set_ylim([minlim,maxlim])
    ax[0][0].legend()

    ax[0][1].scatter(r1ab[::100,0], r1ab[::100,1], s = 10/np.arange(Nsteps+1,0,-100)**0.75, label=r'$r_1$ from $r_a,r_b$', color = 'r')
    ax[0][1].scatter(r2ab[::100,0], r2ab[::100,1], s = 10/np.arange(Nsteps+1,0,-100)**0.75, label=r'$r_2$ from $r_a,r_b$', color = 'b')
    ax[0][1].set_xlabel('x')
    ax[0][1].set_ylabel('y')
    ax[0][1].set_xlim([minlim,maxlim])
    ax[0][1].set_ylim([minlim,maxlim])
    ax[0][1].legend()


    maxlim = np.max([ra.max(), rb.max()])
    minlim = np.min([ra.min(), rb.min()])
    pad = (maxlim - minlim) * 0.05
    maxlim += pad
    minlim -= pad

    ax[1][0].scatter(ra[::100,0], ra[::100,1], s = 10/np.arange(Nsteps+1,0,-100)**0.75, label=r'$r_a$', color = 'r')
    ax[1][0].scatter(rb[::100,0], rb[::100,1], s = 10/np.arange(Nsteps+1,0,-100)**0.75, label=r'$r_b$', color = 'b')
    ax[1][0].set_xlabel('x')
    ax[1][0].set_ylabel('y')
    ax[1][0].set_xlim([minlim,maxlim])
    ax[1][0].set_ylim([minlim,maxlim])
    ax[1][0].legend()

    ax[1][1].plot(np.linalg.norm((r1ab - r1), axis=1), label=r'$1$', color = 'r')
    ax[1][1].plot(np.linalg.norm((r2ab - r2), axis=1), label=r'$2$', color = 'b')
    ax[1][1].set_xlabel('timestep')
    ax[1][1].set_ylabel(r'$\left \| \bar{r}-\bar{r}_{eq}\right \|$')
    ax[1][1].set_yscale('log')
    ax[1][1].legend()

    plt.subplots_adjust(left=0.1,right=0.95,top=0.92,bottom=0.1,wspace=0.3,hspace=0.15)
    plt.savefig(figname)
    plt.show()
    
#####################################################################################

r1, r2 = simulation(r10, v10, r20, v20, m1, m2, F, k, dt, Nsteps)

A, alpha1, alpha2, beta1, beta2, gamma, J, k_new = transformation(m1, m2, p, k, model)

ra0 = alpha1 * r10 + alpha2 * r20
rb0 = beta1 * r10 + beta2 * r20
va0 = alpha1 * v10 + alpha2 * v20
vb0 = beta1 * v10 + beta2 * v20

ra, rb = simulation(ra0, va0, rb0, vb0, A, A, F, k_new, dt, Nsteps)

ra12 = alpha1 * r1 + alpha2 * r2
rb12 = beta1 * r1 + beta2 * r2
r1ab = (beta2 * ra - alpha2 * rb) / J
r2ab = -(beta1 * ra - alpha1 * rb) / J

figtitle = r'$Model:%d, m_1=%.2f,m_2=%.2f,p=%.2f,A=%.2f$' % (model,m1,m2,p,A)
make_plots(r1, r2, ra, rb, r1ab, r2ab, 10, Nsteps, model, figtitle, figname)