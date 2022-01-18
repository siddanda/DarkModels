import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import brentq, curve_fit
from scipy.integrate import quad
import array
import argparse

import emcee
import corner
import sys
from schwimmbad import MPIPool

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
crossSection_M = 3*cm**2/gram
t_age = 10e9*yr

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
gamma = gamma0

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
#Sections
# 1) Dark Matter Potential
# 2) Isothermal Profile
# 3) NFW Profile
# 4) Stitched Profile
# 5) Stitching Functions
# 6) Hernquist Fit
# 7) Rotation Curve Models
# 8) MCMC Functions
#-------------------------------------------------------------------------------------------------------#

#-----------------------------------------Dark Matter Potential-----------------------------------------#

# Interpolation: Solves the Hernquist-specific SIDM differential equation and interpolates it
# getConstants: Returns the constants relevant to the Hernquist-specific SIDM differential equation
# d_Kap: The differential equation, separated into two first order diff eqs.

def Interpolation(rho_dm0, sigma_0, R0, Mtot):

    r = np.linspace(0.001, 500, 500)
    y = np.zeros(len(r))
    
    for i in range(len(y)):
        y[i] = r[i]/(r[i]+R0)
    
    a1 = getConstants(rho_dm0,sigma_0, R0, Mtot)[2]
    potential = odeint(d_Kap,[-a1,0], y, args = (rho_dm0,sigma_0, R0, Mtot))[..., 1]

    h_r = interp1d(y,potential, kind='cubic', fill_value='extrapolate')
    
    return h_r

def getConstants(rho_dm0, sigma_0, R0, Mtot):
    
    rho_b0 = Mtot/(2*math.pi*R0**3)
    
    pot_b0 = -G*Mtot/R0
    
    a0 = 4*math.pi*G*R0**2*rho_dm0/sigma_0**2
    a1 = -pot_b0/sigma_0**2
    
    return rho_b0, a0, a1

#Solving ODEs

def d_Kap(unknowns,y,rho_dm0,sigma_0, R0, Mtot):
    rho_b0, a0, a1 = getConstants(rho_dm0,sigma_0, R0, Mtot)
    x,z = unknowns
    return np.array((-(2/y)*(a1+x) - a0*np.exp(z)/(1-y)**4, x))

#-----------------------------------------Isothermal Profile-----------------------------------------#

# integrand_iso: integrand (profile * area element) for mass enclosed in a given radius in the SIDM isothermal profile
# rho_iso: the SIDM isothermal profile

def integrand_iso(r, rho_dm0, h_r, R0):
    return (4*math.pi*r**2) * rho_iso(r, rho_dm0, h_r, R0)

def rho_iso(r, rho_dm0, h_r, R0):
    
    y1 = r/(r+R0)
    h_y1 = h_r(y1)
    return rho_dm0*np.exp(h_y1)

#---------------------------------------------NFW Profile--------------------------------------------#

# rho_nfw : NFW Profile
# getRvir: used in root finder to determine the radius at which the AVERAGE density of the CDM is Delta_c*p_crit

def rho_nfw(r, rho_s, R_s):
    return rho_s * (r/R_s)**(-1) * (1+r/R_s)**(-2)

def getRvir(Rvir, rho_s, R_s):
    Mvir = 4*math.pi*rho_s*R_s**3*(np.log(1+Rvir/R_s)-Rvir/(R_s+Rvir))
    return Mvir/((4/3)*math.pi*Rvir**3) - Delta_c*p_crit

#------------------------------------------Stitched Profile------------------------------------------#

# integrand_Stitched: integrand (profile * area element) for mass enclosed in a given radius in the stitched profile
# rho_Stitched: the stitched profile

def integrand_Stitched(r, rho_dm0, h_r, r1, R_s, rho_s, R0):
    return (4*math.pi*r**2) * rho_Stitched(r, rho_dm0, h_r, r1, R_s, rho_s, R0)

def rho_Stitched(r, rho_dm0, h_r, r1, R_s, rho_s, R0):
    
    if r<=r1:
        return rho_iso(r, rho_dm0, h_r, R0)
    return rho_nfw(r, rho_s, R_s)        

#-----------------------------------------Stitching Functions-----------------------------------------#

# getStitch: used in root finder to determine the stitching radius
# get_r1: returns the stitching radius
# equations: used in root finder to determine the NFW parameters at r1
# getNFW: returns NFW parameters that give 1) Continuous enclosed mass and 2) Continuous density profile

def getStitch(r1, rho_dm0, Gamma, h_r, R0):    
    return Gamma*rho_iso(r1, rho_dm0, h_r, R0) - 1

def get_r1(rho_dm0, sigma_0, h_r, R0):
    Gamma = (crossSection_M) * (4/np.sqrt(math.pi)) * (sigma_0*km/s)*t_age
    r1 = brentq(getStitch, 0, 10000, args = (rho_dm0, Gamma, h_r, R0))   
    return r1  #These are in kpc

def equations(x, rho_dm0, h_r, r1, R0):
    
    R_s = r1/x

    M_iso = integrate.quad(integrand_iso, 0, r1, args = (rho_dm0, h_r, R0))[0]
    M_nfw = M_iso #This is the first equation

    rho_s = M_nfw/(4*math.pi*(R_s)**3 * (np.log(np.abs(1+r1/R_s))-r1/(r1+R_s)))
    
    return rho_iso(r1, rho_dm0, h_r, R0) - rho_nfw(r1, rho_s, R_s) #This is the second equation

def getNFW(rho_dm0, h_r, r1, R0):
    x = brentq(equations, 0.001, 1000, args = (rho_dm0, h_r, r1, R0))
    R_s = r1/x
    
    M_iso = integrate.quad(integrand_iso, 0, r1, args = (rho_dm0, h_r, R0))[0]
    M_nfw = M_iso
    rho_s = M_nfw/(4*math.pi*(R_s)**3 * (np.log(1+r1/R_s)-r1/(r1+R_s)))
    
    return R_s, rho_s

#------------------------------------------Hernquist Fit------------------------------------------#
#(Only ambiguity is whether or not to fit Mtot or use SPARC data)

# Mass: mass profile of the Hernquist model
# r0_curveFit: returns the best fit scale radius R0 and total baryon mass Mtot to the SPARC (mass profile) data

def Mass(r, R0, Mtot):
    mass = Mtot/(1+R0/r)**2
    return mass

def r0_curveFit(rad, upsilon_d, upsilon_b):
    xdata = rad
    
    M_data = np.zeros(len(rad))
    for i in range(len(rad)):
        M_data[i] = (upsilon_d*v_disk[i]**2 + upsilon_b*v_bulge[i]**2 + v_gas[i]**2)*rad[i]/G
    
    popt, pcov = curve_fit(Mass, xdata, M_data, bounds = ((0,0),(1000,1e15)))
    
    return popt

#--------------------------------------Rotation Curve Models--------------------------------------#

# getModel: the circular velocity function

def getModel(rho_dm0, sigma_0, upsilon_d, upsilon_b, h_r, r1, R0, R_s, rho_s):
    
    #Rotation Curve for DM Halo
    
    v_halo = np.zeros(len(rad))

    for i in range(len(rad)):

        M = (quad(integrand_Stitched,0,rad[i]*kpc,args=(rho_dm0, h_r, r1, R_s, rho_s, R0)))[0]
        v2 = G*M/rad[i]
        
        v_halo[i] = np.sqrt(v2)
        
    ####################################################
    
    #Total Rotation Curve (model)

    v_tot = np.sqrt(v_halo**2 + upsilon_d*v_disk**2 + upsilon_b*v_bulge**2 + v_gas**2)
    
    ####################################################
    
    return v_tot

#------------------------------------------MCMC Functions------------------------------------------#

# log_likelihood: the likelihood function (in log)

# log_prior_1: the prior function that imposes the flat priors on rho_dm0, sigma_0, upsilon_d, upsilon_b
# log_prior_2: the prior function that imposes the additional priors on Vmax and c_vir

# log_probability: the probability function, returns the probability (in log) of a point in parameter space being the true values

def log_likelihood(theta, h_r, r1, R0, R_s, rho_s, x, y, yerr):
    if bulge == False:
        rho_dm0, sigma_0, upsilon_d = theta
        upsilon_b = 0
    else:
        rho_dm0, sigma_0, upsilon_d, upsilon_b = theta
    
    model = getModel(rho_dm0, sigma_0, upsilon_d, upsilon_b, h_r, r1, R0, R_s, rho_s)
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_prior_1(theta):
    if bulge == False:
        rho_dm0, sigma_0, upsilon_d = theta
        upsilon_b = 0

        Gamma = rho_dm0 * (crossSection_M) * (4/np.sqrt(math.pi)) * (sigma_0*km/s) 
        log_Gamma = np.log10(Gamma * 10*Gyr)
        log_sigma = np.log10(sigma_0)  #/km/s

        priorYd = priors.GaussianPrior(np.log10(0.5), 0.1)
        lnp_Yd = priorYd(np.log10(upsilon_d))

        if np.log10(2) < log_Gamma < 5 and np.log10(2) < log_sigma < np.log10(500):
            return lnp_Yd
        return -np.inf

    else:
        rho_dm0, sigma_0, upsilon_d, upsilon_b = theta

        Gamma = rho_dm0 * (crossSection_M) * (4/np.sqrt(math.pi)) * (sigma_0*km/s) 
        log_Gamma = np.log10(Gamma * 10*Gyr)
        log_sigma = np.log10(sigma_0)  #/km/s

        priorYd = priors.GaussianPrior(np.log10(0.5), 0.1)
        lnp_Yd = priorYd(np.log10(upsilon_d))

        priorYb = priors.GaussianPrior(np.log10(0.7), 0.1)
        lnp_Yb = priorYb(np.log10(upsilon_b))
 
        if np.log10(2) < log_Gamma < 5 and np.log10(2) < log_sigma < np.log10(500):
            return lnp_Yd + lnp_Yb
        return -np.inf

def log_prior_2(theta, Roots):
    if bulge == False:
        rho_dm0, sigma_0, upsilon_d = theta
        upsilon_b = 0
    else:
        rho_dm0, sigma_0, upsilon_d, upsilon_b = theta

    R0, r1, R_s, rho_s, Rvir = Roots

    #-----Vmax/Vf Prior-----#

    Vmax = 1.64*R_s*np.sqrt(G*rho_s)

    #-----Concentration-Mass Relation Prior-----#

    c_vir = Rvir/R_s

    Mvir = 4*math.pi*rho_s*R_s**3*(np.log(1+c_vir) - c_vir/(1+c_vir))
    logc_Prior = 1.025 - 0.097*np.log10(Mvir*h/(10**12*M_s))

    prior_c = priors.GaussianPrior(logc_Prior, 0.11)

    log_c = np.log10(c_vir)
    lnp_c = prior_c(log_c)

    xmin, xmax = prior_c.support()

    #-----Abundance Matching Prior-----#

    M_stellar = upsilon_d*(Ltot - Lbulge) + upsilon_b*Lbulge
    x = np.log10(Mvir) - log10_M1Ms
    log10_MstarM1 = eps - np.log10(10**(-alp*x) + 10**(-bet*x)) + gamma*np.exp(-0.5*(x/delta)**2)

    prior_M = priors.GaussianPrior(log10_MstarM1, 0.2)
    log_M = np.log10(M_stellar/M_s) - log10_M1Ms
    lnp_M = prior_c(log_M)

    #-----Log Prior Calculations-----#

    if 1/np.sqrt(2) < Vmax/Vf < np.sqrt(2) and xmin < log_c < xmax:
        return lnp_c + lnp_M
    return -np.inf

def log_probability(theta, x, y, yerr):

    if bulge == False:
        rho_dm0, sigma_0, upsilon_d = theta
        upsilon_b = 0
    else:
        rho_dm0, sigma_0, upsilon_d, upsilon_b = theta

    lp1 = log_prior_1(theta)
        
    if not np.isfinite(lp1):
        return -np.inf
      
    R0, Mtot = r0_curveFit(rad, upsilon_d, upsilon_b)
    h_r = Interpolation(rho_dm0, sigma_0, R0, Mtot)

    try:
        r1 = get_r1(rho_dm0, sigma_0, h_r, R0)
        R_s, rho_s = getNFW(rho_dm0, h_r, r1, R0)
        Rvir = brentq(getRvir, 0.001, 10000, args = (rho_s, R_s))
    except:# Exception as eProb:
        #print('eProb: ' + str(eProb))
        return -np.inf

    Roots = [R0, r1, R_s, rho_s, Rvir]
    lp2 = log_prior_2(theta, Roots)

    if not np.isfinite(lp2):
        return -np.inf
    
    ll = log_likelihood(theta, h_r, r1, R0, R_s, rho_s, x, y, yerr)

    if not np.isfinite(ll):
        return -np.inf
    
    return ll + lp1 + lp2

#-------------------------------------------------------------------------------------------------------#
#--------------------------------------Parallelisation with MPIPool-------------------------------------#
#-------------------------------------------------------------------------------------------------------#

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

#-------------------------------------------------------------------------------------------------------#
#--------------Initialisation: Starting parameters based on rough preliminary minimisation--------------#
#-------------------------------------------------------------------------------------------------------#

    rho_values = np.array(range(1,150,3))*1e7
    sig_values = np.array(range(2,500,10))
    Yd_values = np.array(range(2,10,2))/10

    if bulge==False:
        upsilon_b = 0
        xx,yy,zz = np.meshgrid(Yd_values, rho_values, sig_values)
        positions = np.vstack([yy.ravel(), zz.ravel(), xx.ravel()]).T
    else:
        upsilon_b = 0.7
        xx,ww,yy,zz = np.meshgrid(Yd_values, upsilon_b, rho_values, sig_values)
        positions = np.vstack([yy.ravel(), zz.ravel(), xx.ravel(), ww.ravel()]).T
 
    i=0
    lprob = log_probability(positions[i], rad, v_obs, errV)

    while not np.isfinite(lprob):
        initial = positions[i]
        lprob = log_probability(initial, rad, v_obs, errV)
        i=i+1

    print(Galaxy_Name, initial, lprob)

#--------------------------------------------------------------------------------------------------------#
#----------------------------MCMC Run: 40 Walkers, 5000 Steps, 2500 Discarded----------------------------#
#Saved:
# - Samples
# - Flat Samples
# - Step Plots
# - Corner Plots
# - Rotation Curve Plots
# - Best Fit Parameters with 1\sigma errors
# - BIC Values
#--------------------------------------------------------------------------------------------------------#

#-----------------------MCMC Run-----------------------#

    pos = np.array(initial)
    pos = pos + 1e-4 * np.random.randn(40, len(initial))
    nwalkers, ndim = pos.shape

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(rad, v_obs, errV), pool = pool)
    sampler.run_mcmc(pos, 5000, progress=True)

#------------------------Samples Step Plots------------------------#

    samples = sampler.get_chain()
    np.save('Outputs/Hrnq_Outputs/Samples/'+Galaxy_Name+'.npy', samples)

    if bulge == False:
        labels = ['rho_dm0', "sigma_0", "upsilon_d"]
    else:
        labels = ['rho_dm0', "sigma_0", "upsilon_d", "upsilon_b"]

    fig_steps, axes = plt.subplots(len(labels), figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("step number")

    fig_steps.savefig('Outputs/Hrnq_Outputs/Steps/'+Galaxy_Name)

    tau = sampler.get_autocorr_time(quiet=True)

#-----------------------Flat Samples and Corner Plots-----------------------#

    flat_samples = sampler.get_chain(flat=True)

    np.save('Outputs/Hrnq_Outputs/Flat_Samples/'+Galaxy_Name+'.npy', sampler.get_chain(flat=True))

    flat_samples = sampler.get_chain(discard=2500, flat=True)
    fig_corner = corner.corner(
    flat_samples, labels=labels, show_titles=True)

    fig_corner.savefig('Outputs/Hrnq_Outputs/Corner/'+Galaxy_Name)
    plt.close(fig_corner)

#------------------------Best Fit Parameters with 1\sigma errors-----------------------#
    
    rho_dm0 = np.percentile(flat_samples[:, 0], 50)
    rho16 = np.percentile(flat_samples[:, 0], 16)
    rho84 = np.percentile(flat_samples[:, 0], 84)

    sigma_0 = np.percentile(flat_samples[:, 1], 50)
    sig16 = np.percentile(flat_samples[:, 1], 16)
    sig84 = np.percentile(flat_samples[:, 1], 84)

    upsilon_d = np.percentile(flat_samples[:, 2], 50)
    Yd16 = np.percentile(flat_samples[:, 2], 16)
    Yd84 = np.percentile(flat_samples[:, 2], 84)

    if bulge == True:
        upsilon_b = np.percentile(flat_samples[:, 3], 50)
        Yb16 = np.percentile(flat_samples[:, 3], 16)
        Yb84 = np.percentile(flat_samples[:, 3], 84)
    else:
        upsilon_b = 0

#------------------------Rotation Curve Plots-----------------------#

    R0, Mtot = r0_curveFit(rad, upsilon_d, upsilon_b)
    h_r = Interpolation(rho_dm0, sigma_0, R0, Mtot)
    r1 = get_r1(rho_dm0, sigma_0, h_r, R0)
    R_s, rho_s = getNFW(rho_dm0, h_r, r1, R0)

    v_tot = getModel(rho_dm0, sigma_0, upsilon_d, upsilon_b, h_r, r1, R0, R_s, rho_s)

    v_halo = np.zeros(len(rad))
    for i in range(len(rad)):
        M = (quad(integrand_Stitched,0,rad[i]*kpc,args=(rho_dm0, h_r, r1, R_s, rho_s, R0)))[0]
        v2 = G*M/rad[i]
        v_halo[i] = np.sqrt(v2)

    v_d = np.sqrt(upsilon_d)*v_disk
    v_b = np.sqrt(upsilon_b)*v_bulge

    fig_plot, ax = plt.subplots()
    ax.errorbar(rad, v_obs, yerr=errV, fmt=".k", capsize=0, label = "Total (SPARC Data)")
    ax.plot(rad, v_tot, ":k", label = "Total")
    ax.plot(rad, v_halo, "--", label = "Halo")
    ax.plot(rad, v_gas, "--", label = "Gas")
    ax.plot(rad, v_d,"--", label = "Disk")
    if bulge == True:
        ax.plot(rad, v_b, "--", label = "Bulge")
    ax.set_xlabel("Radius (kpc)")
    ax.set_ylabel("Circular Velocity (km/s)")
    fig_plot.suptitle(str(Galaxy_Name) + 'MCMC Best Fit Rotation Curve')
    ax.legend(bbox_to_anchor=(0, 1), loc='center', ncol=1)

    fig_plot.savefig('Outputs/Hrnq_Outputs/Plots/'+Galaxy_Name)
    plt.close(fig_plot)

    if bulge == True:
        line10 = str(round(upsilon_b, 2))
        line11 = str(round(Yb16, 2))
        line12 = str(round(Yb84, 2))
        theta_ml = [rho_dm0, sigma_0, upsilon_d, upsilon_b]
    else:
        line10 = str(0.0)
        line11 = str(0.0)
        line12 = str(0.0)
        theta_ml = [rho_dm0, sigma_0, upsilon_d]

#------------------------BIC Values-----------------------#

    k = ndim #number of parameters
    n = len(rad) #number of data points

    ln_L = log_likelihood(theta_ml, h_r, r1, R0, R_s, rho_s, rad, v_obs, errV)

    BIC = k * np.log(n) - 2*ln_L

#------------------------Comprehensive Text File-----------------------#

    Data = open("Outputs/Hrnq_Outputs/Hrnq_Data.txt", "a")
    line1 = str(round(rho_dm0/1e7, 2))
    line2 = str(round(rho16/1e7, 2))
    line3 = str(round(rho84/1e7, 2))
    line4 = str(round(sigma_0, 2))
    line5 = str(round(sig16, 2))
    line6 = str(round(sig84, 2))
    line7 = str(round(upsilon_d,2))
    line8 = str(round(Yd16,2))
    line9 = str(round(Yd84,2))
    line10 = line10
    line11 = line11
    line12 = line12
    line13 = str(round(BIC, 2))
    if r1 > rad[-1]:
        Flag = 'Flag'
    else:
        Flag = ''

    Data.writelines(['\n', Galaxy_Name, '\t', line1, '\t', line2, '\t', line3, '\t', line4, '\t', line5, '\t', line6, '\t', line7, '\t', line8, '\t', line9, '\t', line10, '\t', line11, '\t', line12, '\t', line13, '\t', Flag])
    Data.close()

#-------------------------------------------------------------------------------------------------------#