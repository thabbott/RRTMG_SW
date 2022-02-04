import rrtmg_sw_wrapper as sw
import numpy as np
import matplotlib.pyplot as plt

# Constants
cp = 1e3
g = 9.80
Lv = 2.5e6
Rv = 460e0
Rd = 287e0
e0 = 610e0
T0 = 273e0
Mair = 29e0
Mh2o = 18e0
Mco2 = 44e0
h = 6.626e-34
kB = 1.381e-23
c = 3.00e8
NUM_BANDS = 14
NUM_AER = 6


# Atmospheric profiles
def calc_qsat(T, p):
    esat = e0*np.exp(-Lv/Rv*(1/T - 1/T0))
    return esat*Rd/(p*Rv)

def calc_profiles(Ts, Gamma, zt, RH, 
        n = 1000, ps = 101325e0, ztop = 25e3):
    
    # Set height grid
    zf = np.linspace(0, ztop, n+1)
    zc = 0.5*(zf[1:] + zf[:-1])
    
    # Calculate temperatures
    Tt = Ts - Gamma*zt
    Tf = Ts - Gamma*zf
    Tf[zf >= zt] = Tt
    Tc = Ts - Gamma*zc
    Tc[zc >= zt] = Tt

    # Calculate pressures
    pt = ps*((Ts - Gamma*zt)/Ts)**(g/(Rd*Gamma))
    pf = ps*(Tf/Ts)**(g/(Rd*Gamma))
    pf[zf >= zt] = pt*np.exp(-g*(zf[zf >= zt] - zt)/(Rd*Tt))
    pc = ps*(Tc/Ts)**(g/(Rd*Gamma))
    pc[zc >= zt] = pt*np.exp(-g*(zc[zc >= zt] - zt)/(Rd*Tt))

    # Calculate humidity
    qc = RH*calc_qsat(Tc, pc)

    # Return profiles
    return (zc, zf, pc, pf, Tc, Tf, qc)

# Planck function
def calc_B(nu, T):
    return 2*h*c**2*nu**3/(np.exp(h*c*nu/(kB*T)) - 1)

# RRTMG interface
def rrtmg_zeros(shape):
    return np.zeros(shape, dtype = float, order = 'F')

def run_rrtmg_sw(Ts, Gamma, zt, pco2, pch4, RH, 
        n = 1000, ps = 101325e0, ztop = 25e3):

    # Generate profiles
    zc, zf, pc, pf, Tc, Tf, q = calc_profiles(
            Ts, Gamma, zt, RH, n = n, ps = ps, ztop = ztop)

    # Set inputs
    ncol = 1
    nlay = n
    icld = 0
    iaer = 0
    play = rrtmg_zeros((ncol,nlay))
    play[0,:] = pc/1e2
    plev = rrtmg_zeros((ncol,nlay+1))
    plev[0,:] = pf/1e2
    tlay = rrtmg_zeros((ncol,nlay))
    tlay[0,:] = Tc
    tlev = rrtmg_zeros((ncol,nlay+1))
    tlev[0,:] = Tf
    tsfc = rrtmg_zeros((ncol,))
    tsfc[0] = Ts
    h2ovmr = rrtmg_zeros((ncol,nlay))
    h2ovmr[0,:] = q*Rv/Rd
    o3vmr = rrtmg_zeros((ncol,nlay))
    co2vmr = rrtmg_zeros((ncol,nlay))
    co2vmr[0,:] = pco2
    ch4vmr = rrtmg_zeros((ncol,nlay))
    ch4vmr[0,:] = pch4
    n2ovmr = rrtmg_zeros((ncol,nlay))
    o2vmr = rrtmg_zeros((ncol,nlay))
    o2vmr[0,:] = 0.2
    asdir = rrtmg_zeros((ncol,))
    asdir[:] = 0.1
    aldir = rrtmg_zeros((ncol,))
    aldir[:] = 0.1
    asdif = rrtmg_zeros((ncol,))
    asdif[:] = 0.1
    aldif = rrtmg_zeros((ncol,))
    aldif[:] = 0.1
    dyofyr = 0
    adjes = 1.0
    coszen = rrtmg_zeros((ncol,))
    coszen[:] = 0.667
    scon = 510.0
    isolvar = 0
    indsolvar = rrtmg_zeros((2,))
    indsolvar[:] = 1.0
    bndsolvar = rrtmg_zeros((NUM_BANDS,))
    bndsolvar[:] = 1.0
    solcycfrac = 0.0



    inflgsw = 2     # provide cloud water paths and droplet radii,
                    # treat water and ice clouds separately
    iceflgsw = 1    # Ebert & Curry 1992
    liqflgsw = 1    # radius-dependent optical properties
    cldfr = rrtmg_zeros((ncol,nlay))
    taucld = rrtmg_zeros((NUM_BANDS,ncol,nlay))
    ssacld = rrtmg_zeros((NUM_BANDS,ncol,nlay))
    asmcld = rrtmg_zeros((NUM_BANDS,ncol,nlay))
    fsfcld = rrtmg_zeros((NUM_BANDS,ncol,nlay))
    cicewp = rrtmg_zeros((ncol,nlay))
    cliqwp = rrtmg_zeros((ncol,nlay))
    reice = rrtmg_zeros((ncol,nlay))
    reliq = rrtmg_zeros((ncol,nlay))
    tauaer = rrtmg_zeros((ncol,nlay,NUM_BANDS))
    ssaaer = rrtmg_zeros((ncol,nlay,NUM_BANDS))
    asmaer = rrtmg_zeros((ncol,nlay,NUM_BANDS))
    ecaer = rrtmg_zeros((ncol,nlay,NUM_AER))
    
    # Output arrays
    swuflx = rrtmg_zeros((ncol,nlay+1))
    swdflx = rrtmg_zeros((ncol,nlay+1))
    swhr = rrtmg_zeros((ncol,nlay))
    swuflxc = rrtmg_zeros((ncol,nlay+1))
    swdflxc = rrtmg_zeros((ncol,nlay+1))
    swhrc = rrtmg_zeros((ncol,nlay))
    wavenumber_range = rrtmg_zeros((2,))

    # Run RRTM
    SWNS = np.zeros((NUM_BANDS,))
    SWDT = np.zeros((NUM_BANDS,))
    nu_low = np.zeros((NUM_BANDS,))
    nu_high = np.zeros((NUM_BANDS,))
    for iband in range(16,16+NUM_BANDS):
        sw.rrtmg_sw_wrapper(
            ncol, nlay, icld, iaer,             
            play, plev, tlay, tlev, tsfc,                   
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,   
            asdir, asdif, aldir, aldif,                     
            coszen, adjes, dyofyr, scon, isolvar,           
            inflgsw, iceflgsw, liqflgsw, cldfr,             
            taucld, ssacld, asmcld, fsfcld,                 
            cicewp, cliqwp, reice, reliq,                   
            tauaer, ssaaer, asmaer, ecaer,                  
            iband, iband,                           
            swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc,  
            wavenumber_range,                               
            bndsolvar, indsolvar, solcycfrac )
        SWNS[iband-16] = swdflx[0,0] - swuflx[0,0]
        SWDT[iband-16] = swdflx[0,-1]
        nu_low[iband-16] = wavenumber_range[0]
        nu_high[iband-16] = wavenumber_range[1]

    # Return results
    return SWNS, SWDT, nu_low, nu_high

# Initialize
sw.rrtmg_sw_ini_wrapper(cp)

# Compute spectrally-resolved surface net flux
# Use reference state from UChicago RRTM page
# http://climatemodels.uchicago.edu/rrtm/
SWNS, SWDT, nu_low, nu_high = run_rrtmg_sw(
    284.42, 6e-3, 15e3, 400e-6, 1.7e-6, 0.8)
isort = np.argsort(nu_low)
SWNS = SWNS[isort]
SWDT = SWDT[isort]
nu_low = nu_low[isort]
nu_high = nu_high[isort]

print(SWNS)
print(nu_low)
print(nu_high)

# Create arrays for boxy plots
SWNS_nu = SWNS/(nu_high - nu_low)
SWDT_nu = SWDT/(nu_high - nu_low)
SWNS_nu_box = np.stack((SWNS_nu, SWNS_nu), axis = -1).ravel()
SWDT_nu_box = np.stack((SWDT_nu, SWDT_nu), axis = -1).ravel()
nu_box = np.stack((nu_low, nu_high), axis = -1).ravel()

# Plot
plt.rc('font', size = 8)
fig, axes = plt.subplots(figsize = (3.25, 3), nrows = 1, ncols = 1,
        constrained_layout = True, dpi = 200)
axes = np.array(axes)[np.newaxis,np.newaxis]
axes[0,0].plot(nu_box, SWNS_nu_box, 'k-', 
        label = 'surface (net)')
axes[0,0].plot(nu_box, SWDT_nu_box, '-', 
        color = 'gray', label = 'TOA (down)')
axes[0,0].set_title(
        'Surface heating: %.1f W m$^{-2}$' % (np.sum(SWNS)), fontsize = 8)
axes[0,0].set_xlabel(r'$\nu$ (cm $^{-1}$)')
axes[0,0].set_ylabel('Flux (W m$^{-2}$ cm, positive downward)')
axes[0,0].set_xlim([nu_box[0], nu_box[-1]])
axes[0,0].legend(loc = 'upper right', frameon = False)
axes[0,0].set_ylim([0, axes[0,0].get_ylim()[1]])

plt.savefig('test.pdf')
plt.show()
