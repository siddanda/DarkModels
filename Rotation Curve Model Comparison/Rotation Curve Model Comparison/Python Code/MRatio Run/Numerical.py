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

Galaxy_Name = 'ESO079-G014'#Galaxies[j]

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

r = np.round(np.linspace(0, 500, 100001),3)

rad_mid = np.round((rad[0:-1]+rad[1:])/2,3)
r_dMdr = []
for i in range(len(r)):
    if rad_mid[0]<=r[i]<=rad_mid[-1]:
        r_dMdr = np.append(r_dMdr,r[i])

r_after = np.linspace(rad[-1], 500, 100)
r_DP = np.concatenate(([rad[0]],r_dMdr,r_after))

#-------------------------------------------------------------------------------------------------------#
#-----------------------------------------------Functions-----------------------------------------------#
#Sections
# 1) Dark Matter Potential
# 2) Isothermal Profile
# 3) NFW Profile
# 4) Stitched Profile
# 5) Stitching Functions
# 6) Rotation Curve Models
# 8) MCMC Functions
#-------------------------------------------------------------------------------------------------------#

#-----------------------------------------Dark Matter Potential-----------------------------------------#

# Interpolation: Solves the Hernquist-specific SIDM differential equation and interpolates it
# getConstants: Returns the constants relevant to the Hernquist-specific SIDM differential equation
# d_Kap: The differential equation, separated into two first order diff eqs.

def Interpolation(rho_dm0, sigma_0, Density_Function):
    potential = odeint(d,[0,0], r[1:], args = (Density_Function, rho_dm0, sigma_0))[..., 1]   
    psi_tot = interp1d(r[1:], potential, kind='quadratic', fill_value='extrapolate')
    return psi_tot

#Solving ODEs

def d(z, r, Density_Function, rho_dm0, sigma_0):
    x,y = z
    return np.array((-(2/r)*x + 4*math.pi*G*(Density_Function(r) + rho_dm0*np.exp(-y/sigma_0**2)), x))

#-----------------------------------------Isothermal Profile-----------------------------------------#

# integrand_iso: integrand (profile * area element) for mass enclosed in a given radius in the SIDM isothermal profile
# rho_iso: the SIDM isothermal profile

#Isothermal Profile Functions

def integrand_iso(r, rho_dm0, sigma_0, psi_tot):
    return (4*math.pi*r**2) * rho_iso(r, rho_dm0, sigma_0, psi_tot)

def rho_iso(r, rho_dm0, sigma_0, psi_tot):
    return rho_dm0*np.exp(-psi_tot(r)/sigma_0**2)

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

def integrand_Stitched(r, rho_dm0, sigma_0, psi_tot, r1, R_s, rho_s):
    return (4*math.pi*r**2) * rho_Stitched(r, rho_dm0, sigma_0, psi_tot, r1, R_s, rho_s)

def rho_Stitched(r, rho_dm0, sigma_0, psi_tot, r1, R_s, rho_s):
    if r<=r1:
        return rho_iso(r, rho_dm0, sigma_0, psi_tot)
    return rho_nfw(r, rho_s, R_s)

#-----------------------------------------Stitching Functions-----------------------------------------#

# getStitch: used in root finder to determine the stitching radius
# get_r1: returns the stitching radius
# equations: used in root finder to determine the NFW parameters at r1
# getNFW: returns NFW parameters that give 1) Continuous enclosed mass and 2) Continuous density profile

def getStitch(r1, rho_dm0, sigma_0, Gamma, psi_tot):    
    return Gamma*rho_iso(r1, rho_dm0, sigma_0, psi_tot) - 1

def get_r1(rho_dm0, psi_tot, sigma_0):
    Gamma = (crossSection_M) * (4/np.sqrt(math.pi)) * (sigma_0*km/s)*t_age
    r1 = brentq(getStitch, r[0], r[-1], args = (rho_dm0, sigma_0, Gamma, psi_tot))
    return r1  #These are in kpc

def equations(x, rho_dm0, sigma_0, psi_tot, r1):
    
    R_s = r1/x

    M_iso = integrate.quad(integrand_iso, 0, r1, args = (rho_dm0, sigma_0, psi_tot))[0]
    M_nfw = M_iso

    rho_s = M_nfw/(4*math.pi*(R_s)**3 * (np.log(np.abs(1+r1/R_s))-r1/(r1+R_s)))
    
    return rho_iso(r1, rho_dm0, sigma_0, psi_tot) - rho_nfw(r1, rho_s, R_s)

def getNFW(rho_dm0, sigma_0, psi_tot, r1):
    x = brentq(equations, 0.001, 1000, args = (rho_dm0, sigma_0, psi_tot, r1))
    R_s = r1/x
    
    M_iso = integrate.quad(integrand_iso, 0, r1, args = (rho_dm0, sigma_0, psi_tot))[0]
    M_nfw = M_iso

    rho_s = M_nfw/(4*math.pi*(R_s)**3 * (np.log(1+r1/R_s)-r1/(r1+R_s)))
    
    return R_s, rho_s

#--------------------------------------Rotation Curve Models--------------------------------------#

# getModel: the circular velocity function for galaxies

def getModel(rho_dm0, sigma_0, upsilon_d, upsilon_b, psi_tot, r1, R_s, rho_s):

    ####################################################
    
    #Rotation Curve for DM Halo
    
    v_halo = np.zeros(len(rad))

    for i in range(len(rad)):

        M = (quad(integrand_Stitched,0,rad[i]*kpc,args=(rho_dm0, sigma_0, psi_tot, r1, R_s, rho_s)))[0]
        v2 = G*M/rad[i]

        v_halo[i] = np.sqrt(v2)

    ####################################################
    
    #Total Rotation Curve (model)

    v_tot = np.sqrt(v_halo**2 + upsilon_d*v_disk**2 + upsilon_b*v_bulge**2 + v_gas**2)
    
    ####################################################
    
    return v_tot

#--------------------------------------Numerical Functions--------------------------------------#

def getDF(rad, upsilon_d, upsilon_b):
    
    Mreal = (upsilon_d*v_disk**2 + upsilon_b*v_bulge**2 + v_gas**2)*rad/G

    # r < rad[0]:

    p0 = Mreal[0]/(4*math.pi/3 * rad[0]**3)

    # rad[0] < r < rad[-1]:

    dM_dr = np.diff(Mreal)/np.diff(rad)
    dMdrInt = interp1d(rad_mid, dM_dr, kind='linear')
    
    density = 1/(4*math.pi*r_dMdr**2) * dMdrInt(r_dMdr)

    # r > rad[-1]:

    p1 = density[-1]*((rad[-2]+rad[-1])/(2*r_after))**4
    
    # Complete Density Profile

    density = np.concatenate(([p0], density, p1))
    
    Density_Function = interp1d(r_DP, density, kind='linear', fill_value = (p0,density[-1]), bounds_error=False)
    
    return Density_Function
    
#------------------------------------------MCMC Functions------------------------------------------#

# log_likelihood_nb: the likelihood function (in log) for galaxies with no central bulge
# log_likelihood_b: the likelihood function (in log) for galaxies with central bulge

# log_prior_1: the prior function that imposes the flat priors on rho_dm0, sigma_0, upsilon_d, upsilon_b
# log_prior_2: the prior function that imposes the additional priors on c_vir and mass ratios

# log_probability: the probability function, returns the probability (in log) of a point in parameter space being the true values

def log_likelihood(theta, psi_tot, r1, R_s, rho_s, x, y, yerr):
    if bulge == False:
        rho_dm0, sigma_0, upsilon_d = theta
        upsilon_b = 0
    else:
        rho_dm0, sigma_0, upsilon_d, upsilon_b = theta

    model = getModel(rho_dm0, sigma_0, upsilon_d, upsilon_b, psi_tot, r1, R_s, rho_s)
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

    r1, R_s, rho_s, Rvir = Roots

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
    M_gas = 1.33*MHI
    Mratio = (M_stellar + M_gas)/Mvir

    #-----Log Prior Calculations-----#

    if 1/np.sqrt(2) < Vmax/Vf < np.sqrt(2) and xmin < log_c < xmax and Mratio < 0.2:
        return lnp_c
    return -np.inf

def log_probability(theta, rad, v_obs, errV):

    if bulge == False:
        rho_dm0, sigma_0, upsilon_d = theta
        upsilon_b = 0
    else:
        rho_dm0, sigma_0, upsilon_d, upsilon_b = theta

    lp1 = log_prior_1(theta)
        
    if not np.isfinite(lp1):
        return -np.inf

    Density_Function = getDF(rad, upsilon_d, upsilon_b)
    psi_tot = Interpolation(rho_dm0, sigma_0, Density_Function)

    try:
        r1 = get_r1(rho_dm0, psi_tot, sigma_0)
        R_s, rho_s = getNFW(rho_dm0, sigma_0, psi_tot, r1)
        Rvir = brentq(getRvir, 0.001, 10000, args = (rho_s, R_s))
    except:# Exception as eProb:
        #print('eProb: ' + str(eProb))
        return -np.inf

    Roots = [r1, R_s, rho_s, Rvir]
    lp2 = log_prior_2(theta, Roots)

    if not np.isfinite(lp2):
        return -np.inf

    ll = log_likelihood(theta, psi_tot, r1, R_s, rho_s, rad, v_obs, errV)

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
#----------------------------------Initialisation: Starting parameters----------------------------------#
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
    np.save('Outputs/Num_Outputs/Samples/'+Galaxy_Name+'.npy', samples)

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

    fig_steps.savefig('Outputs/Num_Outputs/Steps/'+Galaxy_Name)

    tau = sampler.get_autocorr_time(quiet=True)

#-----------------------Flat Samples and Corner Plots-----------------------#

    np.save('Outputs/Num_Outputs/Flat_Samples/'+Galaxy_Name+'.npy', sampler.get_chain(flat=True))

    flat_samples = sampler.get_chain(discard=2500, flat=True)

    fig_corner = corner.corner(
    flat_samples, labels=labels, show_titles=True)

    fig_corner.savefig('Outputs/Num_Outputs/Corner/'+Galaxy_Name)
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

    Density_Function = getDF(rad, upsilon_d, upsilon_b)
    psi_tot = Interpolation(rho_dm0, sigma_0, Density_Function)
    r1 = get_r1(rho_dm0, psi_tot, sigma_0)
    R_s, rho_s = getNFW(rho_dm0, sigma_0, psi_tot, r1)

    v_tot = getModel(rho_dm0, sigma_0, upsilon_d, upsilon_b, psi_tot, r1, R_s, rho_s)

    v_halo = np.zeros(len(rad))
    for i in range(len(rad)):
        M = (quad(integrand_Stitched,0,rad[i]*kpc,args=(rho_dm0, sigma_0, psi_tot, r1, R_s, rho_s)))[0]
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

    fig_plot.savefig('Outputs/Num_Outputs/Plots/'+Galaxy_Name)
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

    ln_L = log_likelihood(theta_ml, psi_tot, r1, R_s, rho_s, rad, v_obs, errV)

    BIC = k * np.log(n) - 2*ln_L

#------------------------Comprehensive Text File-----------------------#

    Data = open("Outputs/Num_Outputs/Num_Data.txt", "a")
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
    if r1>rad[-1]:
        Flag = 'Flag'
    else:
        Flag = ''

    Data.writelines(['\n', Galaxy_Name, '\t', line1, '\t', line2, '\t', line3, '\t', line4, '\t', line5, '\t', line6, '\t', line7, '\t', line8, '\t', line9, '\t', line10, '\t', line11, '\t', line12, '\t', line13, '\t', Flag])
    Data.close()

#-------------------------------------------------------------------------------------------------------#