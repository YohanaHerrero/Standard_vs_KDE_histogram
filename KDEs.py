from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits
import math
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import seaborn as sns


#I read the fits file to obtain the wanted data
event_filename = get_pkg_data_filename('60fields.fits')
events = Table.read(event_filename, hdu=1)
hdul = fits.open(event_filename)
data = hdul[1].data  

#extract the colums from the table
dec = data['DEC']
ra = data['RA']
redshift = data['Z']

#some specific selection
dec_sel = dec[:1593]
ra_sel = ra[:1593]
redshift_sel = redshift[:1593]

#redshift selection from 3 to 6
select = (redshift_sel >= 3 ) & (redshift_sel <= 6.)
zf = redshift_sel[select]
decf = dec_sel[select] 
raf = ra_sel[select] 


#plot data with color bar
cm = plt.cm.get_cmap('jet') 
fig = plt.figure().add_subplot(111)
plt.scatter(raf,decf, s=10, c=zf, marker='o', cmap=cm)
plt.gca().invert_xaxis()
colorbar=plt.colorbar()
colorbar.set_label('z')
colorbar.ax.tick_params( direction='in')
plt.clim(3., 6.)  
plt.gca().invert_xaxis()
fig.xaxis.set_ticks_position('both')
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_tick_params(direction='in', which='both')
fig.yaxis.set_tick_params(direction='in', which='both')
plt.xlabel("RA", fontsize=14)
plt.ylabel("Dec", fontsize=14)
plt.grid(False)
plt.tight_layout()
#plt.savefig('scatter with color bar data', dpi=500)
plt.show()


#Traditional histogram with 0.03 as bin width
fig = plt.figure().add_subplot(111)
plt.hist(zf, bins=np.arange(3., 6., 0.03), edgecolor = 'k') 
plt.xlabel('z', fontsize = 14)
plt.ylabel(r'N$_{\rm galaxies}$', fontsize = 14)
plt.grid(False)
fig.xaxis.set_ticks_position('both')
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_tick_params(direction='in', which='both')
fig.yaxis.set_tick_params(direction='in', which='both')
plt.tick_params(labelsize = 'large')
plt.grid(False)
#plt.savefig('Traditional histogram', dpi=500)
plt.show()




#Preferentially, it is better to avoid binning -> KDE-histograms 


# This is just a normalised 1D Gaussian, the kernel
def gaussian1d( x, x0, sigma, dx ): 
    #x are redshifts
    #x0 the center of the gaussian
    #sigma is the std deviation of the gaussian 
    #dx the redshift step
    A = 1./(sigma * math.sqrt(2.*math.pi)) * dx
    gaussian = A * np.exp( -0.5 * ((x-x0)/sigma)**2 )
    return gaussian

# The actual KDE
def kdehistogram( redshifts, sigma, agrid ):
    #agrid is the reference grid for the output array (here assumed to be equidistant)
    kdehisto = np.zeros_like(agrid)
    dx = agrid[1]-agrid[0]  
    for zi in redshifts:
        weights = gaussian1d( zi, agrid, sigma, dx )
        kdehisto += weights
    kdehisto *= len(redshifts)/np.sum(kdehisto) #normalise to total number of objects
    return kdehisto

#usage of the KDE
z_step=0.005 
bw=0.005
zgrid = np.arange(3., 6., z_step)  #z_step can be chosen arbitrarily, but should be <= bw for Nyquist sampling.
kdehist = kdehistogram( zf, bw, zgrid )  
zrefbinsize = 2.3548 * bw   #FWHM of KDE Gaussian, by definition FWHM=2.3548*sigma
kdehist *= (zrefbinsize/z_step)  #renormalise to ~ the number of objects within one bin of width = FWHM


#plot the KDE-histogram
fig = plt.figure().add_subplot(111)
plt.xlabel('$z$', fontsize = 14)
plt.ylabel(r'N$_{\rm LAE}$ per $\Delta z=0.012$', fontsize = 14)
plt.plot(zgrid, kdehist,  label="Redshift distribution")
plt.grid(False)
fig.spines['bottom'].set_color('0')
fig.spines['top'].set_color('0')
fig.spines['right'].set_color('0')
fig.spines['left'].set_color('0')
plt.fill_between(zgrid,kdehist,color='lightblue')
fig.xaxis.set_ticks_position('both')
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_tick_params(direction='in', which='both')
fig.yaxis.set_tick_params(direction='in', which='both')
plt.tick_params(labelsize = 'large')
plt.legend(loc = 'upper right').get_frame().set_edgecolor('black')
plt.grid(False)
#plt.savefig('KDE histogram by hand', dpi=500)
plt.show()




#compare KDEs, that by hand and the already implemented in seaborn
fig = plt.figure().add_subplot(111)
plt.xlabel('$z$', fontsize = 14)
plt.ylabel(r'N$_{\rm LAE}$', fontsize = 14)
plt.plot(zgrid, kdehist,  label="By hand")
sns.kdeplot(np.array(zf), bw=0.03, shade=True, kernel='gau')
sns.distplot(zf, bins=np.arange(3., 6., 0.03), rug=False, kde=True, norm_hist=True, hist_kws={"density":False}, kde_kws={"bw": 0.03, "shade": True, "cut":3, "kernel": "gau"},label='Seaborn')
plt.grid(False)
fig.spines['bottom'].set_color('0')
fig.spines['top'].set_color('0')
fig.spines['right'].set_color('0')
fig.spines['left'].set_color('0')
fig.xaxis.set_ticks_position('both')
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_tick_params(direction='in', which='both')
fig.yaxis.set_tick_params(direction='in', which='both')
plt.tick_params(labelsize = 'large')
plt.legend(loc = 'upper right').get_frame().set_edgecolor('black')
plt.grid(False)
#plt.savefig('KDE from seaborn vs by hand',dpi=500)
plt.show()

