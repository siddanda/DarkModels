import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sc
import sympy as sy
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import corner
import argparse
from scipy.optimize import fsolve, brentq, least_squares, minimize, root, curve_fit, differential_evolution

import priors

#Units

kpc = 1
Mpc = 10**3 * kpc
pc = 1e-3*kpc
cm = 3.24078e-19*pc
km = 1e5*cm

M_s = 1
kg = 5.027652e-31*M_s
gram = 1e-3*kg

s = 1
yr = 3.154e+7*s
Gyr = 1e9*yr

Ls = 1 #solar luminosity

G=4.3*10**(-3)*pc/M_s#(km/s)**2

#Cosmological Parameters

H0 = 67.8/Mpc #km/s
Omega_L = 0.692
h = H0/(100/Mpc) #/(km/s)

Delta_c = 18*math.pi**2 - 82*Omega_L - 39*Omega_L**2
p_crit = 3*H0**2 / (8*math.pi*G)

#Abundance Matching Parameters

M0 = 12.054
eps0 = -1.441
alp0 = 1.976
bet0 = 0.465
delta0 = 0.436
gamma0 = -1.016

log10_M1Ms = M0
eps = eps0
alp = alp0
bet = bet0
delta = delta0
gamma = 10**gamma0

#Imported Data

df = pd.read_csv('Rotmod.txt', delim_whitespace=True)

dg = pd.read_csv('SPARC_Data.mrt', delim_whitespace=True, index_col=0)
dg = dg.T

Galaxies = []

for i in range(len(dg.columns)):
    GN = (dg.columns)[i]
    Vflat = (dg[GN].values)[14]
    
    if Vflat != 0.0:
        Galaxies.append(GN)

parser = argparse.ArgumentParser()
parser.add_argument("-j", type=int, default=0)
# Extract dictionary of arguments
args = parser.parse_args()
# Extract parsed parameter
j = args.j

Galaxy_Name = Galaxies[j]

df2 = df.loc[df['Galaxy'] == Galaxy_Name]

rad = df2['Rad'].values
v_obs = df2['Vobs'].values
errV = df2['errV'].values

v_disk = df2['Vdisk'].values
v_gas = df2['Vgas'].values
v_bulge = df2['Vbul'].values

SBdisk = df2['SBdisk'].values
SBbul = df2['SBbul'].values

Data = dg[Galaxy_Name].values

if np.all(v_bulge) == 0:
    bulge = False
else:
    bulge = True

Ltot = Data[6]*10**9*Ls
Lbulge = Data[18]*10**9*Ls
MHI = Data[12]*10**9*M_s
Vf = Data[14] #km/s

#-------------------------------------------------------------------------------------------------------#
#-----------------------------------------------Functions-----------------------------------------------#

#----------------------------------------------DC14 Profile----------------------------------------------#

#DC14 Functions

def B(a, b, x):
    return integrate.quad(lambda t: t**(a-1) * (1-t)**(b-1), 0, x)[0]

def getAlpha(X):
    return 2.94 - np.log10((10**(X+2.33))**(-1.08) + (10**(X+2.33))**(2.29))

def getBeta(X):
    return 4.23 + 1.34*X + 0.26*X**2

def getGamma(X):
    return -0.06 + np.log10((10**(X+2.56))**(-0.68) + 10**(X+2.56))

def DC14(r, c_vir, v_vir, upsilon_d, upsilon_b):
    
    r_2 = v_vir/(np.sqrt(Delta_c/2)*c_vir*H0)
    Rvir = c_vir*r_2
    Mvir = v_vir**2 * Rvir/G
    
    M_stellar = upsilon_d*(Ltot - Lbulge) + upsilon_b*Lbulge
        
    X = min(-1.3,np.log10(M_stellar/Mvir))
        
    alpha = getAlpha(X)
    beta = getBeta(X)
    gamma = getGamma(X)
    
    a = (3-gamma)/alpha
    b = (beta-3)/alpha

    R_s = r_2 * ((beta-2)/(2-gamma))**(1/alpha)
    w = 3 - gamma

    f = lambda alpha,beta,gamma,r: sc.hyp2f1(w/alpha, (beta-gamma)/alpha, 1 + w/alpha, - (r/R_s)**alpha)
    
    x = r/r_2
    return v_vir * np.sqrt((c_vir/x)**(1-w) * f(alpha,beta,gamma,r)/f(alpha,beta,gamma,Rvir))

def getModel(c_vir, v_vir, upsilon_d, upsilon_b):
    
        ####################################################
    
    #Rotation Curve for DM Halo

    v_halo = np.zeros(len(rad))

    for i in range(len(rad)):
        v_halo[i] = DC14(rad[i]*kpc, c_vir, v_vir, upsilon_d, upsilon_b)
        
    ####################################################
    
    #Total Rotation Curve (model)

    v_tot = np.sqrt(v_halo**2 + upsilon_d*v_disk**2 + upsilon_b*v_bulge**2 + v_gas**2)
    
    ####################################################

    return v_tot

def getMax(c_vir, v_vir, upsilon_d, upsilon_b):
    
    r_2 = v_vir/(np.sqrt(Delta_c/2)*c_vir*H0)
    Rvir = c_vir*r_2
    Mvir = v_vir**2 * Rvir/G
    
    M_stellar = upsilon_d*(Ltot - Lbulge) + upsilon_b*Lbulge
        
    X = min(-1.3,np.log10(M_stellar/Mvir))
        
    alpha = getAlpha(X)
    beta = getBeta(X)
    gamma = getGamma(X)
    
    a = (3-gamma)/alpha
    b = (beta-3)/alpha
    
    R_s = r_2 * ((beta-2)/(2-gamma))**(1/alpha)

    f = lambda epsilon: B(a, b, epsilon) - alpha * epsilon**a * (1-epsilon)**b
    epsilon_max = brentq(f, 0.01, 1-1e-3)

    x = (epsilon_max/(1-epsilon_max))**(1/alpha)
    Rmax = x*R_s

    Vmax = DC14(Rmax, c_vir, v_vir, upsilon_d, upsilon_b)

    return Vmax, Rmax

#------------------------------------------MCMC Functions------------------------------------------#

# log_likelihood: the likelihood function (in log)
# log_prior: the prior function (in log)
# log_probability: the probability function, returns the probability (in log) of a point in parameter space being the true values

def log_likelihood(theta, x, y, yerr, bulge):
    if bulge == False:
        c_vir, v_vir, upsilon_d = theta
        upsilon_b = 0
    else:
        c_vir, v_vir, upsilon_d, upsilon_b = theta

    model = getModel(c_vir, v_vir, upsilon_d, upsilon_b) 
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_prior(theta):
    
    if bulge == False:
        c_vir, v_vir, upsilon_d = theta
        upsilon_b = 0
    else:
        c_vir, v_vir, upsilon_d, upsilon_b = theta

    r_2 = v_vir/(np.sqrt(Delta_c/2)*c_vir*H0)
    Rvir = c_vir*r_2
    Mvir = v_vir**2 * Rvir/G

    M_stellar = upsilon_d*(Ltot - Lbulge) + upsilon_b*Lbulge
    M_gas = 1.33*MHI

    X = min(-1.3,np.log10(M_stellar/Mvir))

    #-----Vmax/Vf Prior-----#
    try:
        Vmax, Rmax = getMax(c_vir, v_vir, upsilon_d, upsilon_b)
    except:
        return -np.inf

    #-----Concentration-Mass Relation Prior-----#

    logc_NFW = 1.025 - 0.097*np.log10(Mvir*h/(10**12*M_s))
    logc_Prior = logc_NFW + np.log10(1.0 + 0.00001*np.exp(3.4*(X+4.5)))
    prior_c = priors.GaussianPrior(logc_Prior, 0.11)
    log_c = np.log10(c_vir)
    lnp_c = prior_c(log_c)

    xmin, xmax = prior_c.support()

    #-----Mass-to-Light Ratio Priors-----#

    priorYd = priors.GaussianPrior(np.log10(0.5), 0.1)
    lnp_Yd = priorYd(np.log10(upsilon_d))

    priorYb = priors.GaussianPrior(np.log10(0.7), 0.1)
    lnp_Yb = priorYb(np.log10(upsilon_b))

    #-----Abundance Matching Prior-----#

    M_stellar = upsilon_d*(Ltot - Lbulge) + upsilon_b*Lbulge
    x = np.log10(Mvir/M_s) - log10_M1Ms
    log10_MstarM1 = eps - np.log10(10**(-alp*x) + 10**(-bet*x)) + gamma*np.exp(-0.5*(x/delta)**2)

    prior_M = priors.GaussianPrior(log10_MstarM1, 0.3)
    log_M = np.log10(M_stellar/M_s) - log10_M1Ms
    lnp_M = prior_M(log_M)

    #-----Log Prior Calculations-----#
        
    if bulge==False:
        if 0 < c_vir < 1000 and 0 < v_vir < 500 and 1/np.sqrt(2) < Vmax/Vf < np.sqrt(2) and xmin < log_c < xmax:
            return lnp_c + lnp_Yd + lnp_M
        return -np.inf

    if 0 < c_vir < 1000 and 0 < v_vir < 500 and 1/np.sqrt(2) < Vmax/Vf < np.sqrt(2) and xmin < log_c < xmax:
        return lnp_c + lnp_Yd + lnp_Yb + lnp_M
    return -np.inf

def log_probability(theta, x, y, yerr):

    lp = log_prior(theta)
    ll = log_likelihood(theta, x, y, yerr, bulge)
  
    if not (np.isfinite(lp) and np.isfinite(ll)):
        return -np.inf

    return lp + ll

import emcee
import sys
from schwimmbad import MPIPool
with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    
    c_vir = 1
    v_vir = 1
    upsilon_d = 0.5

    if bulge == False:
        upsilon_b = 0
        initial = [c_vir, v_vir, upsilon_d]
    else:
        upsilon_b = 0.7
        initial = [c_vir, v_vir, upsilon_d, upsilon_b]

    lprob = log_probability(initial, rad, v_obs, errV)

    while (not np.isfinite(lprob)) and c_vir < 1000:
        c_vir += 1
        v_vir = 1
        while v_vir < 500 and (not np.isfinite(lprob)):

            v_vir += 1
            if bulge == False:
                initial = [c_vir, v_vir, upsilon_d]
            else:
                initial = [c_vir, v_vir, upsilon_d, upsilon_b]

            lprob = log_probability(initial, rad, v_obs, errV)

    print(Galaxy_Name, c_vir, v_vir, lprob)

    pos = np.array(initial)
    pos = pos + 1e-4 * np.random.randn(40, len(initial))
    nwalkers, ndim = pos.shape

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(rad, v_obs, errV), pool = pool)
    sampler.run_mcmc(pos, 5000, progress=True);

    samples = sampler.get_chain()
    np.save('Outputs/DC14_Outputs/Samples/'+Galaxy_Name+'.npy',samples)

    if bulge == False:
        labels = ["c_vir", "v_vir", "upsilon_d"]
    else:
        labels = ["c_vir", "v_vir", "upsilon_d", "upsilon_b"]
        
    fig_steps, axes = plt.subplots(len(labels), figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("step number")
    fig_steps.savefig('Outputs/DC14_Outputs/Steps/'+Galaxy_Name)
    plt.close(fig_steps)

    tau = sampler.get_autocorr_time(quiet=True)
    flat_samples = sampler.get_chain(flat=True)

    np.save('Outputs/DC14_Outputs/Flat_Samples/'+Galaxy_Name+'.npy', sampler.get_chain(flat=True))

    flat_samples = sampler.get_chain(discard=2500, flat=True)
    fig_corner = corner.corner(
    flat_samples, labels=labels, show_titles=True)

    fig_corner.savefig('Outputs/DC14_Outputs/Corner/'+Galaxy_Name)
    plt.close(fig_corner)

    c_vir = np.percentile(flat_samples[:, 0], 50)
    c2p5 = np.percentile(flat_samples[:, 0], 2.5)
    c97p5 = np.percentile(flat_samples[:, 0], 97.5)

    v_vir = np.percentile(flat_samples[:, 1], 50)
    v2p5 = np.percentile(flat_samples[:, 1], 2.5)
    v97p5 = np.percentile(flat_samples[:, 1], 97.5)

    upsilon_d = np.percentile(flat_samples[:, 2], 50)
    UD2p5 = np.percentile(flat_samples[:, 2], 2.5)
    UD97p5 = np.percentile(flat_samples[:, 2], 97.5)

    if bulge == False:
        upsilon_b = 0
        UB2p5 = 0
        UB97p5 = 0
    else:
        upsilon_b = np.percentile(flat_samples[:, 3], 50)
        UB2p5 = np.percentile(flat_samples[:, 3], 2.5)
        UB97p5 = np.percentile(flat_samples[:, 3], 97.5)

    v_tot = getModel(c_vir, v_vir, upsilon_d, upsilon_b)

    v_d = np.sqrt(upsilon_d)*v_disk
    v_b = np.sqrt(upsilon_b)*v_bulge

    v_halo = np.zeros(len(rad))
    for i in range(len(rad)):
        v_halo[i] = DC14(rad[i]*kpc, c_vir, v_vir, upsilon_d, upsilon_b)

    fig_plot, ax = plt.subplots()
    ax.errorbar(rad, v_obs, yerr=errV, fmt=".k", capsize=0, label = "Total (SPARC Data)")
    ax.plot(rad, v_tot, ":k", label="Total")
    ax.plot(rad, v_halo, "--", label="Halo")
    ax.plot(rad, v_d,"--", label = "Disk")
    ax.plot(rad, v_gas, "--", label = "Gas")
    if bulge == True:
        ax.plot(rad, v_b, "--", label = "Bulge")
    ax.set_xlabel("Radius(kpc)")
    ax.set_ylabel("Circular Velocity (km/s)")
    fig_plot.suptitle('MCMC Best Fit Rotation Curve')
    ax.legend(bbox_to_anchor=(0, 1), loc='center', ncol=1)

    fig_plot.savefig('Outputs/DC14_Outputs/Plots/'+Galaxy_Name)
    plt.close(fig_plot)

    if bulge == False:
        theta_ml = [c_vir, v_vir, upsilon_d]
    else:
        theta_ml = [c_vir, v_vir, upsilon_d, upsilon_b]

    k = ndim #number of parameters
    n = len(rad) #number of data points

    ln_L = log_likelihood(theta_ml, rad, v_obs, errV, bulge)

    BIC = k * np.log(n) - 2*ln_L

    Data = open("Outputs/DC14_Outputs/DC14_Data.txt", "a")
    line1 = str(round(c_vir, 2))
    line2 = str(round(c2p5, 2))
    line3 = str(round(c97p5, 2))

    line4 = str(round(v_vir, 2))
    line5 = str(round(v2p5, 2))
    line6 = str(round(v97p5, 2))

    line7 = str(round(upsilon_d,2))
    line8 = str(round(UD2p5, 2))
    line9 = str(round(UD97p5, 2))

    line10 = str(round(upsilon_b, 2))
    line11 = str(round(UB2p5, 2))
    line12 = str(round(UB97p5, 2))

    line13 = str(round(BIC, 2))

    Data.writelines(['\n', Galaxy_Name, '\t', line1, '\t', line2, '\t', line3, '\t', line4, '\t', line5, '\t', line6, '\t', line7, '\t', line8, '\t', line9,'\t', line10, '\t', line11, '\t', line12, '\t', line13])
    Data.close()