##############################################################################

print('')
print('----------------------------')
print('summa: MaNGA Galaxy Stacking')
print('Version: 2')
print('Using: MPL-7, v2_4_3')
print('----------------------------')

##############################################################################

# Importing necessary modules
import os                                 # Operating system tools
import glob                               # File name tools
import warnings                           # Python warning trigger
import time                               # Filename tool
import multiprocessing as mp              # Parallel processing for stacking
import copy                               # Module for making deepcopies.
from astropy.io import fits               # Loading .fits files
import numpy as np                        # Numerical analysis
warnings.filterwarnings("ignore")         # Stop all the warnings.

##############################################################################

#============# DEFINITIONS #============#

# NAME: 
#   gal_list
#
# TASK:
#   Identify galaxies with characteristics within given user
#   color-magnitude diagram paremeters. Galaxy characteristics are from the 
#   NSA catalogue, extracting from the most current DRPall file. When used 
#   with summa, gal_list will only produce the DRPall index for each 
#   identified galaxy. Alternatively, one can set a keyword to have gal_list 
#   save the list as a text file.
#
#   Note: One should declare the desired bands for galaxy identification,
#   hard coded into gal_list. The default result in FUV-i color vs. i-band
#   magnitudes. 
#
# PAREMETERS:
#   dim_mag, bright_mag:
#       Desired upper and lower magnitude limits. To be thought of as the 
#       x-range.
#   blue_lim, red_lim:
#       Desired upper and lower color limits. To be thought of as the y-range.
#   print_gal:
#       If 'y', gal_list will create a text file containing the name of each
#       identified galaxy in [plate]-[ifudsgn] format. If 'n', gal_list will
#       only produce an array with the galaxy's DRPall indices. 
#       Default is 'n'.

def gal_list(dim_mag,bright_mag,blue_lim,red_lim, print_gal='n'):
    # High quality galaxies (see notebook 3.3.1)
    gal_idx = np.load(path3+'high_qual_gal_idx.npy')
    # Load desired magnitude bands for both color-mag diagram axes.
    magni = absmag[:,mag_band]
    color = absmag[:,color_blue] - absmag[:,color_red]
    # Determine indicies
    tar_idx1 = np.argwhere((magni <= dim_mag)&(magni >= bright_mag)&
        (color <= red_lim)&(color >= blue_lim))
    tar_idx = np.intersect1d(tar_idx1,gal_idx)
    # Printing galaxy names to stacked_galaxy_list.txt
    if print_gal == 'y':
        gal_names = np.empty(tar_idx.shape[0],dtype=object)
        for i in np.arange(tar_idx.shape[0]):
            gal_names[i] = plateifu[tar_idx][i]
        np.savetxt(path2+'Stacks/'+rgn_name+'/stacked_galaxy_list.txt',
            gal_names,fmt='%s')
    # Done, returning index list! 
    return tar_idx


# NAME: 
#   decode2bit
#
# TASK:
#   Decode quality bitmasks. Output is an array of all the bits
#   encoded within the given quality value. 
#
# PARAMETERS:
#   value:
#       Given DRP3QUAL value. 

def decode2bit(value):
    flags = np.array([])
    while value != 0:
        bit = np.floor(np.log10(value)/np.log10(2))
        value -= 2**bit
        flags = np.concatenate((flags,[bit]))
    return flags


# NAME:
#   regions
#
# TASK:
#   Given a value and an array representing bin borders, determine the
#   correct bin a value would be a part of. For example, given a value
#   of 3 and bin borders [0,2,4,6], the value would be in the second bin.
#   In this example, there are only 3 possible bins.
#
# PARAMETERS:
#   value:
#       Given value. To be used by summa as the b/a and radial distance.
#   borders:
#       User designated border values.
#   prefix:
#       Region name prefix. For example, when summa uses this to determine
#       inclination bins, the prefix is set to 'i'. 
#

def regions(val,borders,prefix):
    found = 'n'
    while found == 'n':
        for i in np.arange(borders.shape[0]):
            # If the value is larger than the largest bin edge then this
            # part will trigger. 
            if i == len(r_borders)-1:
                value1 = 'X'
                found = 'y'
                break
            if (val >= borders[i])&(val <= borders[i+1]):
                value1 = prefix+str(i)
                found = 'y'
                break
    return value1


# NAME: 
#   combine_galaxy
#
# TASK:
#   Combine all the spectra in a single galaxy in the desired way.
#   This is the main engine behind the stacking routine. There are
#   a lot of definitions that are made before calling this. 
#
# PARAMETERS:
#   tar_idx: integer
#       DRPall index of galaxy to be stacked.   

def combine_galaxy(tar_idx):
    g = np.int(tar_idx)
    # Loading .fits files. summa will move on to the next galaxy if there
    # are problems with the quality bitmasks.
    # Loading DRP
#============================================================================#
    # Make sure that this is the correct path to your DRP LOGCUBES!
    #drp_gal = glob.glob(path+plateifu[g]+'/*LOGCUBE*7.fits*')
    drp_gal = glob.glob(path+'/drp/'+np.str(plate[g])+
        '/stack/*'+plateifu[g]+'*LOGCUBE.fit*')
#============================================================================#
    # Skips galaxy and continues with stack if data is missing
    if len(drp_gal) == 0:
        print('DRP for '+plateifu[g]+' is missing.')
        return
    # Loads data and quality flag
    obj = fits.open(drp_gal[0])
    flags = decode2bit(obj[0].header['DRP3QUAL'])
    # Skips galaxy and continues with data if flag suggests poor data.
    # I check this before coming up with the galaxy indices, so this should
    # not be necessary unless the user is using a specific galaxy list. 
    if np.setdiff1d(flags,okay_bit).shape[0] >= 1:
        if verbose == 'y':
            print('Skipping galaxy due to DRP quality issues: '+
                plateifu[g])
        return
    # Loading DAP
#============================================================================#
# Make sure that this is the correct path to your DAP MAPS!
    #dap_gal = glob.glob(path+plateifu[g]+'/*MAPS*HYB10*')
    dap_gal = glob.glob(path+'/dap/HYB10-GAU-MILESHC/'+np.str(plate[g])+'/'+
        ifudsgn[g]+'/*MAPS*.fit*')
#============================================================================#
    # Skips galaxy and continues with stack if data is missing
    if len(dap_gal) == 0:
        print('DAP for '+plateifu[g]+' is missing.')
        return
    # Loads data and quality flag
    dap = fits.open(dap_gal[0])
    flags = decode2bit(dap[0].header['DAPQUAL'])
    # Skips galaxy and continues with data if flag suggests poor data.
    # I check this before coming up with the galaxy indices, so this should
    # not be necessary unless the user is using a specific galaxy list. 
    if np.setdiff1d(flags,okay_bit_dap).shape[0] >= 1:
        if verbose == 'y':
            print('Skipping galaxy due to DAP quality issues: '+
                plateifu[g])
        return

    # Print progress..
    if verbose == 'y':
        print('Stacking galaxy: ',plateifu[g])

    # Loading DRP extensions
    flux = obj[1].data
    ivar = obj[2].data
    mask = obj[3].data
    disp = obj[4].data
    logwave = np.log10(obj[6].data)
    obj.close()

    # Loading DAP extensions
    coord = dap[2].data[1,:,:]     # Coordinates for each spaxel in R/R_{eff}
    if (half_light == 'y'): 
        if (hl_units == 'a'):
            coord *= petro_hlr_r[g]    # If we want this in arcsecs instead...
    snr_g = dap[5].data            # mean g-band weighted SNR (weighted?)
    if vel_regi == 'stellar':
        # LOS stellar velocity relative to cz in km/s
        stevel = dap[15].data           
        stevel_mask = dap[17].data      # Quality mask for above.
    #===============================================================#    
        # LOS stellar velocity dispersion in km/s, corrected. 
        #stesig = np.sqrt((np.square(dap[18].data)-np.square(dap[21].data)))       
        #stesig_mask = dap[20].data      # Quality mask for above.
    #===============================================================#

    else:
        # LOS emission-line vel. disp. as per above.
        emvel = dap[36].data[18,:,:]   
        emvel_mask = dap[38].data[18,:,:]    # Quality mask for above. 
    #===============================================================#
        #emsig = np.sqrt((np.square(dap[39].data[18,:,:])-
        #                 np.square(dap[42].data[18,:,:])))   # Sigma values
        #emsig_mask = dap[41].data[18,:,:]                    # Mask for above
    #===============================================================#
    # Closing .fits file. 
    dap.close()

#===============================================================#
    # Need to add some code here to deal with keeing DAP sigmas
    # for viewing after stack and for comparing to sigma measured
    # from the final spectrum.
#===============================================================#
     # If the inclination slicing toggle was set to 'y', this labels all
    # spectra from this galaxy appropriately.
    if incli == 'y':
        incli1 = regions(petro_ba[g],i_borders,'i')
    else:
        incli1 = 'n'

    # Determining velocity registration shift.
    velcorr = np.log10((stevel/c)+1+redshift[g])

    # Determining spaxels with good velocity measurements.
    gidx = np.argwhere((snr_g >= sncut)&(stevel_mask == 0))

    # Creating dictionaries
    # Copying empty dictionary onto other dictionaries...
    nspec = copy.deepcopy(base)
    # Weighted mean arrays
    wmnnum = copy.deepcopy(base)
    wmnden = copy.deepcopy(base)
    wmn = copy.deepcopy(base)
    wmnerr = copy.deepcopy(base)
    wmnerrnum = copy.deepcopy(base)
    wmndispnum = copy.deepcopy(base)
    wmndispden = copy.deepcopy(base)
    wmndisp = copy.deepcopy(base)
    # As above, but for the raw mean. 
    rmnnum = copy.deepcopy(base)
    rmnden = copy.deepcopy(base)
    rmn = copy.deepcopy(base)
    rmndispnum = copy.deepcopy(base)
    rmndispden = copy.deepcopy(base)
    rmndisp = copy.deepcopy(base)
    # Other values for each spectrum, currently not being used...
    spec_vals = copy.deepcopy(base)

    for i in gidx:
        xval,yval = i[0],i[1]
        corflux = flux[:,xval,yval]
        corivar = ivar[:,xval,yval]
        cordisp = disp[:,xval,yval]
        
        # Average of the inverse variance to correct error vector later..
        corivar_avg = np.average(corivar[np.isfinite(corivar)])
        
        if subpix == 'y':
            # Subpixeling onto a finer grid
            corflux = np.interp(fwave,logwave,corflux)
            corivar = np.interp(fwave,logwave,corivar)
            cordisp = np.interp(fwave,logwave,cordisp)

        # Number of spectra in each channel (and weights for raw mean)
        numspec = np.ones(len(corflux))
        numspec[np.argwhere(mask[:,xval,yval] != 0)] = 0

        # Creating corrected wavelength array.
        corwave = fwave - velcorr[xval,yval]

        # Creating flux normalization constant and correcting flux
        nwave = np.argwhere((corwave > nregn_low)&
            (corwave < nregn_hi)).flatten()
        nflux = corflux[nwave]
        norm = nflux[(nflux!=0)&(np.isfinite(nflux))]
        if len(norm) == 0:
            # This only happens when every value inside the desired
            # normalization range is indefinite or equal to 0.
            # Skipping spaxel.
            continue
        norm = np.average(norm)
        corflux /= norm
        corivar *= (norm**2)
        
        # Determines where the new wavelength array fits into the master
        # from the left.
        base_shift = np.modf((corwave[0]-wave_ma[0])/del_wave)
        # Creates modifiers
        s_shift = base_shift[1]
        fmodL = base_shift[0]
        fmodR = 1 - fmodL

        # Fixing arrays to be per pixel. This is effectively a weight.
        # in which left pixels are weighted less than right pixels because of
        # the logarithmic scaling. 
        corflux *= dens_fix
        corivar /= dens_fix**2
        numspec *= dens_fix
        cordisp *= dens_fix
        
        # Transforms flux/ivar into arrays of same size as master, making
        # sure to conserve array values.
        # Flux
        beta = ( ((fmodL*np.concatenate(([0],corflux))) +
            (fmodR*np.concatenate((corflux,[0])))) / 
                ((fmodL*np.concatenate(([0],dens_fix))) +
            (fmodR*np.concatenate((dens_fix,[0])))) )
        # Inverse variance for weighting
        alpha = ((((fmodL*np.concatenate(([0],dens_fix))) + 
                  (fmodR*np.concatenate((dens_fix,[0]))))**2) / ( 
                ((fmodL**2)*np.concatenate(([0],1/corivar))) + 
                ((fmodR**2)*np.concatenate((1/corivar,[0])))) )
        # Number of spectra
        aura = ( ((fmodL*np.concatenate(([0],numspec))) +
            (fmodR*np.concatenate((numspec,[0])))) / 
                ((fmodL*np.concatenate(([0],dens_fix))) +
            (fmodR*np.concatenate((dens_fix,[0])))) )
        # Dispersion
        gamma = np.sqrt(( ((fmodL**2)*np.concatenate(([0],cordisp**2))) + 
                ((fmodR**2)*np.concatenate((cordisp**2,[0]))) ) / (
                ((fmodL**2)*np.concatenate(([0],dens_fix**2))) + 
                ((fmodR**2)*np.concatenate((dens_fix**2,[0]))) ))
        
        # Determines how many trailing zeros need to be concatenated to
        # the new array before it fits with the master.
        endfix = np.int(wave_ma.shape[0] - np.int(s_shift) - 
            beta.shape[0])
        
        # Creating final flux/ivar arrays. Note that this whole process
        # smoothes each array slightly, failing to conserve array values
        # by ~1%. 
        mf = np.concatenate((np.zeros(np.int(s_shift)),beta,
            np.zeros(endfix)))
        mi = np.concatenate((np.zeros(np.int(s_shift)),alpha,
            np.zeros(endfix)))
        # Variance
        err_corr = corivar_avg / np.average(alpha[alpha!=0])  # Normalization
        ei = np.concatenate((np.zeros(np.int(s_shift)),1/(err_corr*alpha),
            np.zeros(endfix)))
        ei[np.argwhere(np.isinf(ei) == True)] = 0    # Removing infinity.
        # Normalizing the weight
        mi /= np.average(mi[mi!=0])
        # Array to keep track of how many spectra are in in each channel
        nspec1 = np.concatenate((np.zeros(np.int(s_shift)),aura,
             np.zeros(endfix)))
        # Keeping track of the dispersion..
        cor_lsf = np.concatenate((np.zeros(np.int(s_shift)),gamma,
             np.zeros(endfix)))

         # If slicing by half light radius, categorize the current spaxel
            # into the correct region.
        if half_light == 'y':
            radi1 = regions(coord[xval,yval],r_borders,'r')
            if radi1 == 'X':
                continue
        else:
            radi1 = 'n'

        # Final vectors are: mf,mi,nspec1,cor_lsf
        # Stacking in each desired region..
        # Number of spectra
        nspec[radi1+'-'+incli1] += nspec1
        # Raw Mean
        rmnnum[radi1+'-'+incli1] += (mf*nspec1)
        rmnden[radi1+'-'+incli1] += (nspec1)
        rmndispnum[radi1+'-'+incli1] += (cor_lsf*nspec1)**2
        rmndispden[radi1+'-'+incli1] += (nspec1)**2
        # Weighted Mean
        wmnnum[radi1+'-'+incli1] += (mf*mi*nspec1)
        wmnden[radi1+'-'+incli1] += (mi*nspec1)
        wmndispnum[radi1+'-'+incli1] += (cor_lsf*mi*nspec1)**2
        wmndispden[radi1+'-'+incli1] += (mi*nspec1)**2
        wmnerrnum[radi1+'-'+incli1] += (ei*(nspec1*mi)**2)

    #===============================================================#
        # Need to add some code here to deal with keeing DAP sigmas
        # for viewing after stack and for comparing to sigma measured
        # from the final spectrum.
    #===============================================================#

    if save_every_galaxy == 'y':
        spec_vals[plateifu[g]] = np.zeros(len(nspec1))

    return (nspec,rmnnum,rmnden,rmndispnum,rmndispden,
        wmnnum,wmnden,wmndispnum,wmndispden,wmnerrnum,spec_vals)


##############################################################################

#Starting the clock!
start_time = time.time()

# Path to input files
#path3 = path
path3 = '/usr/data/dizzy/machuca/summa_inputs/'

# Read user parameters.
inp = np.genfromtxt(path3+'input.txt',dtype='U')
verbose = inp[0]                                    # Verbose toggle
sncut = np.int(inp[1])                              # S/N Cutoff
agn = inp[2]                                        # AGN toggle
vel_regi = inp[3]                                   # Vel. type
half_light = inp[4][0]                              # Radial slicing toggle
if half_light == 'y':
    hl_units = inp[4][2:]
    r_borders = np.fromstring(inp[5],dtype='float64',sep=',')
else:
    r_borders = np.array([0])
incli = inp[6]                                      # b/a slicing toggle         
if incli == 'y':
    i_borders = np.fromstring(inp[7],dtype='float64',sep=',')
else:
    i_borders = np.array([0])
subpix = inp[8][0]                                  # Subpixeling toggle
rgn_name = inp[9]                                   # Region name
gal_list_tog = inp[10][0]                           # List toggle.
if gal_list_tog == 'y':
    bin_num = np.int(inp[10][2:])
else:
    raw_list = inp[11]                              # List name
    print('Stacking from: '+raw_list)
    usr_list = np.array([np.genfromtxt(path3+raw_list,dtype='U')]).flatten()

#=============================================================#
# Changeable parameters otherwise hard-coded into summa:
# Speed of light..
c = 2.998e5

# Ignorable DRP3QUAL bits. Include any bits describing problematic DRP-3d 
# products you want to include in the stacking the stacking. The bitmask
# description can be found at:
# https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-7/metadata
# If you wish to only use the best quality data, set the array to only 
# include the value 0, like okay_bit = np.array([0])
okay_bit = np.array([0])

# Ignorable DAPQUAL bits. Same as above. If you want to only use the best
# quality data, set the array to an only include the value zero.
okay_bit_dap = np.array([0])

# If using gas velocities, this keyword determines which emission line
# velocity summa will use for the velocity registration. The velocities are
# identical for every line in MPL-7, but just in case that changes in the
# future...
em_line = 18

# During the stack, summa normalizes each given spectrum by the median
# flux value within a given wavelength range. The user can change this
# wavelength range below. The values are to be in angstroms. 
nregn_low = np.log10(5450)
nregn_hi = np.log10(5500)

# Subpixeling integer. Interpolates all the arrays to fit a wavelength vector
# in which the delta(lambda) is 1e-4 / (this_integer).
subpixint = 4

# Path for the DRP and DAP data
#path = '/home/machuca/Wisconsin/Research/fits/'
path = '/usr/data/dizzy2/mab/manga/mpl7_data/'
# Path to place created files into. Will create the folders in this directory.
#path2 = '/home/machuca/Wisconsin/Research/'
path2 = '/usr/data/dizzy/machuca/'

# Sigma tracking is an option for if you wish summa to maintain an
# array of all the Ha sigmas from the DAP for all the spectra that are
# stacked. If yes, a bunch of other extra files will be saved.
# Specifically: the Ha sigma of all the encompasing spectra, the
# g-band SNR for those spectra, and the average inverse variance weight
# around Ha for each spectra. 
sigma_tracking = 'n'

# Save each individual galaxy stack? This is mostly to check for any crazy
# errors that might have come up. Recommend to leave as 'n'.
save_every_galaxy = 'n'

# Making directory
os.makedirs(path2+'Stacks/'+rgn_name,exist_ok=True)

# Loading necessary values from the NSA catalogue..
drp = fits.open(path3+'drpall-v2_4_3.fits')
plate = drp[1].data['plate']
ifudsgn = drp[1].data['ifudsgn']
plateifu = drp[1].data['plateifu']
redshift = drp[1].data['nsa_z']
if (half_light == 'y'):
    petro_hlr_r = drp[1].data['nsa_elpetro_th50_r']
if (incli == 'y'):
    petro_ba = drp[1].data['nsa_elpetro_ba']
if gal_list_tog == 'y':
    # Declare bands to be used for galaxy identification in color-magnitude
    # space. Available broadbands are FNugriz, in order from 0-6. 
    mag_band = 5
    color_blue = 0
    color_red = 5
    # Loading magnitudes
    absmag = drp[1].data['nsa_elpetro_absmag']
drp.close()

# Creates galaxy index list from CMD values
if gal_list_tog == 'y':
    usr_rgn = np.loadtxt(path3+'bin_info.txt')
    # Printing stacked region (CMD)
    print('Stacking region '+np.str(bin_num)+': '+np.str([usr_rgn[bin_num,0],
        usr_rgn[bin_num,1],usr_rgn[bin_num,2],usr_rgn[bin_num,3]]))
    # Making list
    tar_idx = gal_list(usr_rgn[bin_num,0],usr_rgn[bin_num,1],
        usr_rgn[bin_num,2],usr_rgn[bin_num,3])
# Creates galaxy index list using a user-supplied text file in which the
# galaxy names are in plate-ifudsgn format. 
else:
    tar_idx = np.zeros(usr_list.shape[0])
    for i in np.arange(usr_list.shape[0]):
        tar_idx[i] = np.argwhere(plateifu == usr_list[i])

#Creating master arrays for final stacked spectra
del_wave = 1e-4
if subpix == 'y':
    del_wave /= subpixint
wave_ma = np.arange(np.log10(3125),np.log10(10400),del_wave)

# Empty dictionary
# Using this as a template to then make the rest of the dictionaries.
base = {}

# Creating every possible dictionary key based on the given bounderies..
if (incli=='n')&(half_light=='n'):
    base['n-n'] = np.zeros(len(wave_ma))
elif (incli=='n'):
    for i in np.arange(0,len(r_borders)-1,1):
        base['r'+str(i)+'-n'] = np.zeros(len(wave_ma))
elif (half_light=='n'):
    for j in np.arange(0,len(i_borders)-1,1):
        base['n-i'+str(j)] = np.zeros(len(wave_ma))
else:
    for i in np.arange(0,len(r_borders)-1,1):
        for j in np.arange(0,len(i_borders)-1,1):
            base['r'+str(i)+'-i'+str(j)] = np.zeros(len(wave_ma))

# Determining necessary values in order to correct flux from 
# erg/s/cm^2/ang/spaxel to erg/s/cm^2/pixel/spaxel before stacking.
# This is important because summa will seek to velocity register
# every spectrum before stacking them together. The final product
# will return to erg/s/cm^2/ang/spaxel.
# fwave here is identical to every galaxy wavelength vector. 
fwave = np.arange(3.5589,4.0151,del_wave)
dens_fix = 10**(fwave+(del_wave/2)) - 10**(fwave-(del_wave/2))

if verbose == 'y':
    print('Stacking '+np.str(len(tar_idx))+' galaxies..')

# Parallel Processing Stacking
if __name__ == "__main__":
    pool = mp.Pool()
    results = pool.map_async(combine_galaxy, tar_idx)
    pool.close()
    pool.join()

# In order to check each individual galaxy stack..
if save_every_galaxy == 'y':
    np.save(path2+'Stacks/'+rgn_name+'/every_gal_idx',tar_idx)
    np.save(path2+'Stacks/'+rgn_name+'/every_gal',results.get())

# Copying base dictionary onto other dictionaries...
nspec = copy.deepcopy(base)
# Weighted mean arrays
wmnnum = copy.deepcopy(base)
wmnden = copy.deepcopy(base)
wmn = copy.deepcopy(base)
wmnerr = copy.deepcopy(base)
wmnerrnum = copy.deepcopy(base)
wmndispnum = copy.deepcopy(base)
wmndispden = copy.deepcopy(base)
wmndisp = copy.deepcopy(base)
# As above, but for the raw mean. 
rmnnum = copy.deepcopy(base)
rmnden = copy.deepcopy(base)
rmn = copy.deepcopy(base)
rmndispnum = copy.deepcopy(base)
rmndispden = copy.deepcopy(base)
rmndisp = copy.deepcopy(base)
# Other values for each spectrum, currently not being used...
spec_vals = copy.deepcopy(base)

# Galaxy number, array number, dictionary key.
# Galaxy number is 0 through len(results.get()) - 1.
# Array numbers: nspec:0, rmnnum:1, rumnden:2, rmndispnum:3, rmndispden:4, 
# wmnnum:5, wmnden:6, wmndispnum:7, wmndispden:8.
# Dictionary keys are like 'r0-i0' or similar.
for i in np.arange(len(results.get())):
# Getting final results from this stack..
    finarr = (results.get())[i]
    # Putting into final arrays...
    for j in nspec.keys():
        nspec[j] += finarr[0][j]
        rmnnum[j] += finarr[1][j]
        rmnden[j] += finarr[2][j]
        rmndispnum[j] += finarr[3][j]
        rmndispden[j] += finarr[4][j]
        wmnnum[j] += finarr[5][j]
        wmnden[j] += finarr[6][j]
        wmndispnum[j] += finarr[7][j]
        wmndispden[j] += finarr[8][j]
        wmnerrnum[j] += finarr[9][j]

# Finding final means..
for i in nspec.keys():
    # Raw Mean
    rmnden[i][rmnden[i] == 0] = 1
    rmndispden[i][rmndispden[i] ==0] = 1
    rmn[i] = rmnnum[i]/rmnden[i]
    rmndisp[i] = np.sqrt(rmndispnum[i]/rmndispden[i])
    # Weighted Mean
    wmnden[i][wmnden[i] == 0] = 1
    wmndispden[i][wmndispden[i] == 0] = 1
    wmn[i] = wmnnum[i]/wmnden[i]
    wmndisp[i] = np.sqrt(wmndispnum[i]/wmndispden[i])
    wmnerr[i] = np.sqrt(wmnerrnum[i]/(wmnden[i]**2))

#===============================================================#
# Need to add some code here to include a header on one of the 
# dictionaries so that it's easy to figure out what went into a
# particular stack. 
#===============================================================#

# Ending the clock!
end_time = time.time()
print('Done! This took '+np.str(np.int(np.ceil(end_time - start_time)))+
    ' seconds!')

# Saving files
np.save(path2+'Stacks/'+rgn_name+'/wave',10**wave_ma)
np.save(path2+'Stacks/'+rgn_name+'/nspec',nspec)
np.save(path2+'Stacks/'+rgn_name+'/wmn',wmn)
np.save(path2+'Stacks/'+rgn_name+'/rmn',rmn)
np.save(path2+'Stacks/'+rgn_name+'/wmndisp',wmndisp)
np.save(path2+'Stacks/'+rgn_name+'/rmndisp',rmndisp)
np.save(path2+'Stacks/'+rgn_name+'/wmnerr',wmnerr)

# Printing some stuff for the end..
if verbose == 'y':
    print('----------------------------')
    print('Stack Overview:')
    print(np.str(len(results.get()))+' total combined galaxies.')
    if gal_list_tog == 'y':
    # Printing stacked region (CMD)
        print('Stacking region '+np.str(bin_num)+': '+
            np.str([usr_rgn[bin_num,0],usr_rgn[bin_num,1],usr_rgn[bin_num,2],
            usr_rgn[bin_num,3]]))
    else:
        print('Stacked from user list: '+raw_list)
    print('SNR cut-off: ',np.str(sncut))
    if half_light == 'y':
        print('Radial bin edges: ',r_borders)
    else:
        print('No radial binning.')
    if incli == 'y':
        print('b/a bin edges: ',i_borders)
    else:
        print('No inclination-based binning.')
    print('')

# Done!

##############################################################################
#
# Version 2: 04/2019
# Included subpixeling onto a finer grid and a couple of other options. 
#
# Version 1: 02/2019
# Essentially completely re-done.
# Multiprocessing included, major fixes to stacking methods.
#
# Version 0: 
# Old summa.py from 06/2018
#
##############################################################################