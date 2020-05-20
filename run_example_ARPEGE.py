# Running the wake parameterization over the AMMA case, as simulated with ARPEGE-Climat

import sys, os
sys.path.append('wakelib')

import numpy as np

import netCDF4

import matplotlib.pyplot as plt

import wakelib

###################
# if lplot is True, some plot of input/output variables will be performed and saved in the directories below.
lplot = True
if lplot:
    rep_images_in = './images/input'
    if not(os.path.exists(rep_images_in)):
        os.makedirs(rep_images_in)
    rep_images_out = './images/output'
    if not(os.path.exists(rep_images_out)):
        os.makedirs(rep_images_out)

###################
# Defining the wake parameterization parameters
parameters = {'timestep': 300.,\
              'wapecut': 1.,\
              'sigmad': 0.02,\
              'hwmin': 10.,\
              'sigmaw_max': 0.4,\
              'dens_rate': 0.1,\
              'wdensmin': 1.e-14,\
              'stark': 0.33,\
              'alpk':0.25,\
              'wdens_ref': (8.e-12,8.e-12),\
              'coefgw': 4.,\
              'tau_cv': 4000.,\
              'crep_upper': 0.9,\
              'crep_sol': 1.0,\
              'delta_t_min': 0.2,\
              'rzero':5000.}

# Initialize the wake parameters
wakelib.init(parameters)

###################
# Preparing the input variables for the wake parameterization

nlev = 91 # Number of full levels
nlevp1 = nlev + 1 # Number of half levels

# time axis
time = np.arange(6*3600,24*3600,parameters['timestep'])/3600.
nstep, = time.shape

fin = 'files/Out_arp631.wake2.sigconv_CMIP6.wakeThermo.LPBLE_L91_300s_AMMA.nc'

var2read = ['znatsurf','p','ph','pi','te','qe','omgbe','dtdwn','dqdwn','amdwn','amup','dta','dqa','wgen','sigd0','cin','deltatw0','deltaqw0','sigmaw0','awdens0','wdens0']

datain = {} # contains input data
f = netCDF4.Dataset(fin,'r')
for var in var2read:
    #print var
    datain[var] = f['wake_'+var][:]

    if lplot: # plotting input data
        if len(datain[var].shape) == 2:
            nt,nlev_loc = datain[var].shape
            if nlev_loc == nlev:
                lev = f['pf'][:,::-1]/100.
                datain['pfax'] = f['pf'][:,::-1]/100.
            elif nlev_loc == nlevp1:
                lev = f['ph'][:,::-1]/100.
                datain['phax'] = f['ph'][:,::-1]/100.
            else:
                print 'nlev_loc unexpected:', nlev_loc
                sys.exit()
            l3d = True
        else:
            nt, = datain[var].shape
            l3d = False

        if l3d:
            time_loc = np.tile(time,(nlev_loc,1))
            cs = plt.pcolormesh(time_loc,np.transpose(lev),np.transpose(datain[var]))
            plt.colorbar(cs)
            plt.ylim(1000.,0.)
            plt.ylabel('Pressure (hPa)')
        else:
            cs = plt.plot(time,datain[var])

        plt.xlabel('Hour (UTC)')
        plt.title('{0}'.format(var))
 
        plt.savefig('{0}/{1}.png'.format(rep_images_in,var))
        plt.close()

f.close()

# Need to rename inout variables
datain['deltatw'] = datain['deltatw0']
datain['deltaqw'] = datain['deltaqw0']
datain['sigmaw']  = datain['sigmaw0']
datain['awdens']  = datain['awdens0']
datain['wdens']   = datain['wdens0']

# Add a few more variables
datain['timestep'] = parameters['timestep']

###################
# Executing the wake parameterization over nstep
dataout = wakelib.itere(datain,nstep=nstep,update=True)

###################
# Plotting output, if required
if lplot:
    dt = parameters['timestep']/3600.
    for var in dataout.keys():
        l3d = len(dataout[var].shape) == 2
        if l3d:
            nt_loc,nlev_loc = dataout[var].shape
            time_loc = np.tile(time,(nlev_loc,1))
            if nlev_loc == nlev:
                lev = datain['pfax']
            elif nlev_loc == nlevp1:
                lev = datain['phax']
            else:
                print 'nlev0 unexpected:',nlev, 'for var:', var
                sys.exit()

            cs = plt.pcolormesh(time_loc+dt,np.transpose(lev),np.transpose(dataout[var]))
            plt.colorbar(cs)
            plt.ylim(1000.,0.)
            plt.ylabel('Pressure (hPa)')
        else:
            plt.plot(time+dt,dataout[var],'b-',label='offline')
            try:
                plt.plot(time,datain[var],'k-',label='original')
            except:
                pass
 
            plt.legend()

        plt.xlabel('Hour (UTC)')
        plt.title('{0}'.format(var))

        plt.savefig('{0}/{1}.png'.format(rep_images_out,var))
        plt.close()


###################
# Some final cleaning
wakelib.clean()
