##############################################################################
#
# NAME: summa
#
# UPDATES: 
#   https://github.com/c-machuca/summa
#
# PURPOSE:
#   summa is a flexible spectral stacking routine for MaNGA data. 
#   The final spectrum can be the result of a variety of averages, 
#   including the raw mean, weighted mean, and median. summa accepts 
#   galaxy input lists given in typical MaNGA convention (plate-ifudsgn). 
#   For help with galaxy list making, see gal_list below.. 
#
#   Key optional features include:
#   1)  Subpixeling data to a finer grid to artificially increase
#       spectral resolution of final spectrum. 
#   2)  Slicing input spectra by half-light radius and host galaxy b/a.
#   3)  Produce 2D histograms of all input spectra for visualization.
#
# CALLING:
#   The desired input parameters should be declared in a separate file
#   titled 'input.txt'. An example detailing the appropriate order and 
#   syntax of the input file can be found below:
#
#   y           # Verbose toggle. Print progress?
#   wmn         # Type of average?
#   6           # Signal-to-noise cutoff? Minimum of 3. 
#   y           # Discard AGN candidates?
#   stellar	    # Use stellar or gas velocities?
#   y,a	        # Two toggles. Slice by half-light radius? Re or arcsecs?
#   0,.5,1      # If 'y' above, borders of desired radius regions.
#   y           # Slice by galaxy b/a?
#   0,5,10,15   # If 'y' above, borders of desired b/a (inclination) regions. 
#   y,2	        # Subpixel? If y, give an integer to describe by how much.
#   n,400,500   # Produce 2D histograms? If y, give desired wavelength range.
#   full        # Region name?
#   y,7	        # Use pre-discriminated ranges for galaxy list making?
#   gal.txt	    # If 'n' above, give list of galaxies to be stacked.
#
#   More information on input parameters:
#   Line 1: This is a simple verbose toggle. If 'y', summa will print 
#           when it begins a new galaxy and when it skips a galaxy due
#           to problematic flags in its error mask. 
#   Line 2: Determines the type of average summa will attempt. The current
#           choices are raw mean (rmn), weighted mean (wmn), and median (med).
#           Sigma-clipping is a planned feature in the near future.
#   Line 3: Value used to determine at what point summa will begin to ignore
#           noisy spectra. summa will use the signal-to-noise in the g-band
#           to compare to the given value and determine progression.
#           It's recommended to use a minimum value of 3 here or otherwise
#           problems have been observed. Preferably, use 5. 
#   Line 4: Toggle to discard any AGN candidates in the galaxy stacking list.
#           AGN candidates are from the work by Wylezalek et al, 2018. 
#   Line 5: Stellar vs. gas velocity for use in velocity registration of
#           spectra. Velocities are from the MaNGA DAP products. For gas
#           velocities, summa uses H-alpha measurements.
#   Line 6: If first toggle is yes, summa will determine the half-light radius
#           position of the current spaxel and assign it a region (from a 
#           predetermined list of possibilities). The half-light radius used 
#           will come from the NSA elliptical petrosian in the r-band. If you 
#           would like the the radial regions to be in half-light radii, the 
#           second toggle should be 'r'. If you want them to be in arcseconds,
#           it should be 'a'. 
#   Line 7: If the above line is yes, these are the defined border regions for
#           the user's desired radial bins. For example, if you want two bins 
#           - one between 0 and 3", and one between 3 and 10", you would write
#           0,3,10 for this line. 
#   Line 8: If yes, summa will use the galaxy's b/a (from the NSA catalog) to
#           determine to which set a spectrum should be stacked to.
#   Line 9: If the above line is yes, these are the defined border regions for
#           the user's desired b/a bins. 
#   Li. 10: If yes, summa will resample each spectrum to a finer grid before
#           stacking the set. The new grid will be finer by the integer
#           declared by the user, i.e. a '2' would result in a grid with
#           half-sized pixels. 
#   Li. 11: If yes, summa will produce a 2D histogram to help the user
#           visualize the full set of stacked spectra. This task requires
#           a large amount of memory and will therefore increase the total
#           computational time significantly. It is recommended to therefore
#           limit the stack to a smaller wavelength region, given in nm. 
#   Li. 12: A region name to be used in creating the output file. A new folder
#           with this name will be made in the current directory. 
#   Li. 13: The user can declare pre-determined regions in a text file titled
#           'bin_info.txt' and summa will use the gal_list routine below to 
#           create a galaxy list. The necessary integer declares which line 
#           summa will read. For example, an example of how to declare 
#           regions is below:
#
#           # dim_mag   bright_mag  blue_lim   red_line
#           -17.25      -18.25      1.5	        2.5
#           -20,25      -21.25      3.5	        4.5
#
#   Li. 14: If the user already has a list of galaxies to be stacked, the file
#           should be declared here. There should be one galaxy per line.
#
# OUTPUT:
#   Outputs are created by the Numpy save module (.npy) and can be restored
#   via the Numpy load module.  Most arrays are held in the dictionary data 
#   structure, each accesible by a set of keys (strings) with the naming
#   convention '(Half light radial region)-(inclination region)'. If your stack
#   does not slice by half light radius or by inclination, they will be
#   accesible in wmn['n-n']. All other letter possibilities are described
#   below. 
#
#       Letter guides:
#       i#   :   Inclination regions, 0 represents the first of the user bins.
#                If the user indicates 4 regions, i3 will be the last.
#       r#   :   Radial regions. As above. 
#       n    :   Not sliced in this way.
#
#   wmn/rmn: dictionary
#       Weighted or raw mean arrays. The units are arbitrary, as the original
#       spectra are normalized during the stacking.
#
#   wmnerr/rmnerr: dictionary
#       Weighted or raw mean error arrays. Keep in mind that these are rough
#       approximations. summa does not attempt to correctly account for
#       covariance in order to do the error propagation correctly. These
#       numbers should only be trusted to within a factor of 2 or 3. The
#       weighted mean errors are of better quality than the raw mean errors.
#
#   nspec: dictionary
#       Number of spectra that make up each wavelength pixel.
#
#       Note:
#          This dictionary will also be populated by various other keys that
#          correspond to extra information like the signal-to-noise cut used,
#          or the inclination/radial borders defined. Think of this as a
#          'header', like headers in .fits files. 
#
#   wave: array
#       Wavelength array of produced spectra. Units are in angstroms. 
#
#   resolution: dictionary
#       1-sigma resolution of each stacked spectrum. Units in angstroms.
#
# HARD-CODED PARAMETERS:
#   These are parameters that are hard-coded into summa. These can be changed,
#   but at the user's discretion. They can be found below. 
#
# UPDATE HISTORY:
#   V1.0: Created by Camilo Machuca. UW - Madison, March 2018. 
#   V1.1: Enabled raw mean, June 5th.
#
##############################################################################

print('')
print('summa: MaNGA Galaxy Stacking')
print('----------------------------')

# Load necessary packages.
import numpy as np
from astropy.io import fits
import time
import glob
import os
import warnings

#=============================================================#

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
#   Note: One should declare the desired bands for galaxy identification, hard
#   coded into gal_list. The default result in NUV-i color vs. i-band
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
#       only produce an array with the galaxy's DRPall indices. Default is 'n'.

def gal_list(dim_mag,bright_mag,blue_lim,red_lim, print_gal='n'):
    # Declare bands to be used for galaxy identification in color-magnitude
    # space. Available broadbands are FNugriz, in order from 0-6. 
    mag_band = 5
    color_blue = 1
    color_red = 5
    # Load desired DRPall extensions
    drp = fits.open(path+'drpall-v2_3_1.fits')
    if print_gal == 'y':
        plate = drp[1].data['plate']
        ifudsgn = drp[1].data['ifudsgn']
    absmag = drp[1].data['nsa_elpetro_absmag']
    drp.close()
    # Load desired magnitude bands for both color-mag diagram axes.
    magni = absmag[:,mag_band]
    color = absmag[:,color_blue] - absmag[:,color_red]
    # Determine indicies
    idx1 = np.argwhere((magni <= dim_mag)&(magni >= bright_mag))
    idx2 = np.argwhere((color <= red_lim)&(color >= blue_lim))
    tar_idx = np.intersect1d(idx1,idx2)
    if print_gal == 'y':
        gal_names = np.empty(tar_idx.shape[0],dtype=object)
        for i in np.arange(tar_idx.shape[0]):
            gal_names[i] = np.str(plate[tar_idx][i])+'-'+ifudsgn[tar_idx][i]
        np.savetxt('galaxy_list.txt',gal_names,fmt='%s')
    return tar_idx

# NAME: 
#   decode2bit
#
# TASK:
#   Decode DRP3QUAL quality bitmask. Output is an array of all the bits
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
#   header
#
# TASK:
#   Fill dictionary with important stack information, like the radial and
#   inclination bin borders, the signal-to-noise cut-off value, etc, etc.
#   Consider it like a header for .fits files.  
#
# PARAMETERS:
#   flux_array:
#       Dictionary to be filled.  

def header(flux_array):
    flux_array['sncut'] = sncut
    flux_array['ngal'] = np.array([n,prob_gal.shape[0]])
    if half_light == 'y':
        flux_array['r_borders'] = r_borders
    if incli == 'y':
        flux_array['i_borders'] = i_borders
    if gal_list_tog == 'y':
        flux_array['bin_rgn'] = [usr_rgn[bin_num,0],usr_rgn[bin_num,1],
        usr_rgn[bin_num,2],usr_rgn[bin_num,3]]
    else:
        flux_array['gal_list'] = usr_list
    if subpix == 'y':
        flux_array['subpix_int'] = subpixint
    if histogram == 'y':
        flux_array['wave_range'] = np.array([rgn_low,rgn_hi])
    return flux_array

# NAME: 
#   estimate
#
# TASK:
#   Predict how long it will take summa to finish. 
#
# PARAMETERS:
#   tar_idx: array
#       DRPall indices of galaxies to be stacked.   

def estimate(tar_idx):
    ifu_size = np.array(ifudsgn[tar_idx],dtype='float')
    n127 = len(np.argwhere(ifu_size > 12700))
    n19 = len(np.argwhere(ifu_size < 3700))
    n37 = len(np.argwhere((ifu_size > 3700)&(ifu_size < 6100)))
    n61 = len(np.argwhere((ifu_size > 6100)&(ifu_size < 9100)))
    n91 = len(np.argwhere((ifu_size > 9100)&(ifu_size < 12700)))
    tpredict = np.int(((n127*10.4)+(n19*2.4)+(n37*3.7)+(n61*5.6)+(n91*7.8)) 
        * 1.05 / 60)
    return tpredict


#=============================================================#

# Starting the clock!
start_time = time.time()

# Read user parameters.
inp = np.genfromtxt('input.txt',dtype='U')
verbose = inp[0]
avgtype = inp[1]
sncut = np.int(inp[2])
agn = inp[3]
vel_regi = inp[4]
half_light = inp[5][0]
if half_light == 'y':
    hl_units = inp[9][2:]
    r_borders = np.fromstring(inp[6],dtype='float64',sep=',')
else:
    r_borders = np.array([0])
incli = inp[7]
if incli == 'y':
    i_borders = np.fromstring(inp[8],dtype='float64',sep=',')
else:
    i_borders = np.array([0])
subpix = inp[9][0]
if subpix == 'y':
    subpixint = np.int(inp[9][2:])
histogram = inp[10][0]
if histogram == 'y':
    rgn_low = np.int(inp[10][2:inp[10][2:].find(',')+2])
    rgn_hi = np.int(inp[10][inp[10][2:].find(',')+3:])
rgn_name = inp[11]
gal_list_tog = inp[12][0]
if gal_list_tog == 'y':
    bin_num = np.int(inp[12][2:])
else:
    raw_list = inp[13]
    print('Stacking from: '+raw_list)
    usr_list = np.array([np.genfromtxt(raw_list,dtype='U')])

#=============================================================#
# Changeable parameters otherwise hard-coded into summa:

# Ignorable DRP3QUAL bits. Include any bits describing problematic DRP-3d 
# products you want to include in the stacking the stacking. The bitmask
# description can be found at:
# https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-6/metadata
# If you wish to only use the best quality data, set the array to only 
# include the value 0, like okay_bit = np.array([0])
okay_bit = np.array([0])

# Ignorable DAPQUAL bits. Same as above. If you want to only use the best
# quality data, set the array to an only include the value zero.
okay_bit_dap = np.array([0])

# If using gas velocities, this keyword determines which emission line
# velocity summa will use for the velocity registration. The velocities are
# identical for every line in MPL-6, but just in case that changes in the
# future...
em_line = 18

# During the stack, summa normalizes each given spectrum by the median
# flux value within a given wavelength range. The user can change this
# wavelength range below. The values are to be in angstroms. 
nregn_low = 6750
nregn_hi = 6900

# Path for DRP and DAP data
#path = '/usr/data/dizzy2/mab/manga/mpl6_data/'
path = '/Users/machuca/Wisconsin/Manga/testing_data/'

# Debug mode for when the produced stacks are less than ideal. Plots
# a spectrum of the stack after each new galaxy. 
debug = 'n'
if debug == 'y':
    import matplotlib.pyplot as plt

# Sigma tracking is an option for if you wish summa to maintain an
# array of all the Ha sigmas from the DAP for all the spectra that are
# stacked. If yes, a bunch of other extra files will be saved.
# Specifically: the Ha sigma of all the encompasing spectra, the
# g-band SNR for those spectra, and the average inverse variance weight
# around Ha for each spectra. 
sigma_tracking = 'y'

#=============================================================#

# Declare constants
c = 2.998*10**5

# Creates directory for stack.
fol_name = (half_light+np.str(r_borders.shape[0]-1)+incli+
    np.str(i_borders.shape[0]-1))
os.makedirs('stacks/'+rgn_name,exist_ok=True)
os.makedirs('stacks/'+rgn_name+'/'+fol_name,exist_ok=True)

# Creates galaxy indexlist.
if gal_list_tog == 'y':
    usr_rgn = np.loadtxt('bin_info.txt')
    print('Stacking region '+np.str(bin_num)+': '+np.str([usr_rgn[bin_num,0],
        usr_rgn[bin_num,1],usr_rgn[bin_num,2],usr_rgn[bin_num,3]]))
    tar_idx = gal_list(usr_rgn[bin_num,0],usr_rgn[bin_num,1],
        usr_rgn[bin_num,2],usr_rgn[bin_num,3])
else:
    drp = fits.open(path+'drpall-v2_3_1.fits')
    plate = drp[1].data['plate']
    ifudsgn = drp[1].data['ifudsgn']
    drp.close()
    tar_idx = np.zeros(usr_list.shape[0])
    for i in np.arange(usr_list.shape[0]):
        tar_idx[i] = np.intersect1d(np.argwhere(plate ==
            np.int(usr_list[i][0:usr_list[i].find('-')])),
            np.argwhere(ifudsgn == usr_list[i][usr_list[i].find('-')+1:]))
ngal = tar_idx.shape[0]

#	Identifies AGN. 
if agn == 'y':
    agn_idx = np.load('agn_idx_v2_3_1.npy')
    tar_idx = np.setdiff1d(tar_idx,agn_idx)
    if verbose == 'y':
        print(np.str(np.int(ngal-tar_idx.shape[0]))+
            ' galaxies found to be AGN candidates.')

# Loads necessary values from the NSA catalogue.
drp = fits.open(path+'drpall-v2_3_1.fits')
plate = drp[1].data['plate']
ifudsgn = drp[1].data['ifudsgn']
redshift = drp[1].data['nsa_z']
if (half_light == 'y'):
    petro_hlr_r = drp[1].data['nsa_elpetro_th50_r']
if (incli == 'y'):
    petro_ba = drp[1].data['nsa_elpetro_ba']
drp.close()

# Print prediction of run time.
tar_idx = np.array(tar_idx,dtype=np.int)
tpredict = estimate(tar_idx)
print('Predicted run time: approx. '+np.str(tpredict)+' minutes.')

# Creating master arrays for final stacked spectra
del_wave = .0001
if subpix == 'y':
    del_wave = del_wave/subpixint
wave_ma = np.arange(np.log10(3125),np.log10(10400),del_wave)
nspec = {}
reso = {}
total_reso = {}
if avgtype == 'wmn':
    wmnnum = {}
    wmnden = {}
    wmn = {}
    wmnerr = {}
elif avgtype == 'rmn':
    rmnnum = {}
    rmnden = {}
    rmnerr = {}
    rmn = {}

if debug == 'y':
    r = np.arange(len(wave_ma))[(10**wave_ma >= 5000)&(10**wave_ma <= 6000)]
if sigma_tracking == 'y':
    ha_sig = {}
    sig_noi = {}
    sig_weight = {}

# Loading each individual galaxy. 
n = 1
prob_gal = np.array([])
for g in tar_idx:
    g = np.int(g)
    # Verbose toggle to mark progress.
    if verbose == 'y':
        cur_time = time.time()
        print('Stacking galaxy number '+np.str(n)+' out of '
            +np.str(tar_idx.shape[0])+'. Current run time: '+
            np.str(np.int(cur_time - start_time))+' seconds.')
    # Loading .fits files. summa will move on to the next galaxy if there
    # are problems with the quality bitmasks.
    #drp_gal = glob.glob(path+'drp/'+np.str(plate[g])+'/stack/manga-'+
    #    np.str(plate[g])+'-'+ifudsgn[g]+'-LOGCUBE.fit*')
    drp_gal = glob.glob(path+'/manga-'+
        np.str(plate[g])+'-'+ifudsgn[g]+'-LOGCUBE.fit*')
    if len(drp_gal) == 0:
        print('DRP for '+np.str(plate[g])+'-'+ifudsgn[g]+' is missing.')
        n += 1
        prob_gal = np.concatenate((prob_gal,[g]))
        continue
    obj = fits.open(drp_gal[0])
    flags = decode2bit(obj[0].header['DRP3QUAL'])
    if np.setdiff1d(flags,okay_bit).shape[0] >= 1:
        if verbose == 'y':
            print('Skipping galaxy due to DRP quality issues: '+
                np.str(plate[g])+'-'+ifudsgn[g])
            n += 1
            prob_gal = np.concatenate((prob_gal,[g]))
        continue 
    #dap_gal = glob.glob(path+'dap/SPX-GAU-MILESHC/'+np.str(plate[g])+'/'+
    #    ifudsgn[g]+'/manga-'+np.str(plate[g])+'-'+ifudsgn[g]+
    #    '-MAPS-SPX-GAU-MILESHC.fit*')
    dap_gal = glob.glob(path+'manga-'+np.str(plate[g])+'-'+ifudsgn[g]+
        '-MAPS-SPX-GAU-MILESHC.fit*')
    if len(dap_gal) == 0:
        print('DAP for '+np.str(plate[g])+'-'+ifudsgn[g]+' is missing.')
        n += 1
        prob_gal = np.concatenate((prob_gal,[g]))
        continue
    maps = fits.open(dap_gal[0])
    flags = decode2bit(maps[0].header['DAPQUAL'])
    if np.setdiff1d(flags,okay_bit_dap).shape[0] >= 1:
        if verbose == 'y':
            print('Skipping galaxy due to DAP quality issues: '+
                np.str(plate[g])+'-'+ifudsgn[g])
            n += 1
            prob_gal = np.concatenate((prob_gal,[g]))
        continue
    
    # Loading DRP extensions
    flux = obj[1].data
    ivar = obj[2].data
    mask = obj[3].data
    disp = obj[4].data
    logwave = np.log10(obj[6].data)
    obj.close()
    
    # Loading DAP extensions
    if vel_regi == 'stellar':
        velo = maps[15].data
        velomask = maps[17].data
    else:
        velo = maps[36].data[em_line,:,:]
        velomask = maps[38].data[em_line,:,:]
    if half_light == 'y':
        # Multiplying by half-light radius to get regions in terms of arcsecs.
        if hl_units == 'a':
            coord = maps[2].data[1,:,:]*petro_hlr_r[g]
        else:
            coord = maps[2].data[1,:,:]
    snr = maps[5].data

    # Sigma tracking..
    if sigma_tracking == 'y':
        ha_sigma = maps[39].data[em_line,:,:]
        ha_sigma_corr = maps[42].data[em_line,:,:]
        warnings.filterwarnings('ignore')
        corr_sigma = np.sqrt(ha_sigma**2 - ha_sigma_corr**2)
        sig_weight_wave = np.ndarray.flatten(np.argwhere(
                (10**wave_ma > 6555)&(10**wave_ma < 6570) == True))
    maps.close()
    
    # If the inclination slicing toggle was set to 'y', this labels all
    # spectra from this galaxy appropriately.
    if incli == 'y':
        incli1 = regions(petro_ba[g],i_borders,'i')
    else:
        incli1 = 'n'
    
    # Determining necessary values in order to correct flux from 
    # erg/s/cm^2/ang/spaxel to erg/s/cm^2/pixel/spaxel before stacking.
    # This is important because summa will seek to velocity register
    # every spectrum before stacking them together. The final product
    # will return to erg/s/cm^2/ang/spaxel. 
    dens_fix = (10**(logwave+(del_wave/2)) - 10**(logwave-(del_wave/2)))
    
    # Determining velocity registration shift.
    zred = redshift[g]
    velcorr = np.log10((velo/c)+1+zred)
    
    # Determining spaxels with good velocity measurements.
    gvel = np.argwhere(velomask == 0)
    
    # Begin for-loop for spectral stacking.
    for i in gvel:
        # Eliminating spaxels with signal-to-noise ratio lower than given cut.
        if snr[i[0],i[1]] >= sncut:
            # Currecting flux and the inverse variance
            corflux = flux[:,i[0],i[1]] * dens_fix
            corivar = ivar[:,i[0],i[1]] / (dens_fix**2)
            # Creating array describing the number of spectra
            num_spec = np.ones(len(corflux))
            # Fixing array with mask to take into account bad pixels for
            # weighting purposes. Similarly, fixing inveres variance..
            num_spec[np.argwhere(mask[:,i[0],i[1]] != 0)] = 0
            corivar *= num_spec
            # Dispersion
            cordisp = disp[:,i[0],i[1]]
            # Creating corrected wavelength array.
            corwave = logwave - velcorr[i[0],i[1]]

            # Creating flux normalization constant and correcting flux/ivar.
            nwave = np.ndarray.flatten(np.argwhere((corwave > 
                np.log10(nregn_low))& (corwave < np.log10(nregn_hi)) == True))
            nflux = corflux[nwave]
            norm = nflux[(nflux!=0)&(np.isfinite(nflux))]
            if len(norm) == 0:
                # This only happens when every value inside the desired
                # normalization range is indefinite or equal to 0.
                continue
            norm = np.median(norm)
            corflux /= norm
            corivar *= norm**2

            # Determines where the new wavelength array fits into the master
            # from the left.
            base_shift = np.modf((corwave[0]-wave_ma[0])/del_wave)
            # Creates modifiers
            s_shift = base_shift[1]-1
            fmodL = base_shift[0]
            fmodR = 1 - fmodL

            # Transforms flux/ivar into arrays of same size as master, making
            # sure to conserve array values.  
            beta = ((fmodL*np.concatenate(([0],corflux))) +
                (fmodR*np.concatenate((corflux,[0]))))
            alpha = ((fmodL*np.concatenate(([0],corivar))) +
                (fmodR*np.concatenate((corivar,[0]))))
            # Number of spectra and dispersion..
            aura = ((fmodL*np.concatenate(([0],num_spec))) +
                (fmodR*np.concatenate((num_spec,[0]))))
            gamma = ((fmodL*np.concatenate(([0],cordisp))) +
                (fmodR*np.concatenate((cordisp,[0]))))
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
            # Array to keep track of how many spectra are in in each channel
            nspec1 = np.concatenate((np.zeros(np.int(s_shift)),aura,
                 np.zeros(endfix)))
            # Keeping track of the dispersion..
            cor_lsf = np.concatenate((np.zeros(np.int(s_shift)),gamma,
                 np.zeros(endfix)))

            # If slicing by half light radius, categorize the current spaxel
            # into the correct region.
            if half_light == 'y':
                radi1 = regions(coord[i[0],i[1]],r_borders,'r')
                if radi1 == 'X':
                    continue
            else:
                radi1 = 'n'

            # Create dictionary keys to save data into.
            if radi1+'-'+incli1 not in nspec.keys():
                nspec[radi1+'-'+incli1] = np.zeros(wave_ma.shape[0])
                reso[radi1+'-'+incli1] = np.zeros(wave_ma.shape[0])
                if avgtype == 'wmn':
                    wmnnum[radi1+'-'+incli1] = np.zeros(wave_ma.shape[0]) 
                    wmnden[radi1+'-'+incli1] = np.zeros(wave_ma.shape[0])
                elif avgtype == 'rmn':
                    rmnnum[radi1+'-'+incli1] = np.zeros(wave_ma.shape[0]) 
                    rmnden[radi1+'-'+incli1] = np.zeros(wave_ma.shape[0])
                
                # Sigma tracking stuff..
                if sigma_tracking == 'y':
                    ha_sig[radi1+'-'+incli1] = np.array([])
                    sig_noi[radi1+'-'+incli1] = np.array([])
                    sig_weight[radi1+'-'+incli1] = np.array([])

            # Stacking in regions
            # Number of spectra
            nspec[radi1+'-'+incli1] += nspec1
            # Tracking the resulting spectral resolution
            reso[radi1+'-'+incli1] += (cor_lsf**2)
            if avgtype == 'wmn':
                wmnnum[radi1+'-'+incli1] += (mf*mi) 
                wmnden[radi1+'-'+incli1] += mi
            elif avgtype == 'rmn':
                rmnnum[radi1+'-'+incli1] += (mf*nspec1)
                rmnden[radi1+'-'+incli1] += (1/mi)

            # Sigma tracking stuff..
            if sigma_tracking == 'y':
                ha_sig[radi1+'-'+incli1] = np.concatenate((
                        ha_sig[radi1+'-'+incli1],[corr_sigma[i[0],i[1]]]))
                sig_noi[radi1+'-'+incli1] = np.concatenate((
                        sig_noi[radi1+'-'+incli1],[snr[i[0],i[1]]]))
                sig_weight[radi1+'-'+incli1] = np.concatenate((
                        sig_weight[radi1+'-'+incli1],
                        [np.average(mi[sig_weight_wave])]))

    # Galaxy counter            
    n += 1
    # Debug plotting
    if debug == 'y':
        if (incli1 == 'i2'):
            plt.plot(wave_ma[r],wmnnum['r2-i2'][r]/wmnden['r2-i2'][r],'b')
            plt.plot(wave_ma[r],wmnnum['r2-i1'][r]/wmnden['r2-i1'][r],'r')
            plt.title(np.str(plate[g])+'-'+ifudsgn[g])
            plt.show()

# Next few lines will give divide by zero/NaN warnings. This is expected.
# Ignoring warnings..
warnings.filterwarnings('ignore')

# Calculating appropriate average and correcting fluxes back 
# to erg/s/cm^2/ang/spaxel. The 4.2 multiplicative factor in the errors is
# an attempt to correct for covariance.
dens_fix_ma = (10**(wave_ma+(del_wave/2)) - 10**(wave_ma-(del_wave/2)))
for i in nspec.keys():
    total_reso[i] = np.sqrt(reso[i]/nspec[i])
    if avgtype == 'wmn':
        wmn[i] = (wmnnum[i]/wmnden[i]) / dens_fix_ma
        wmnerr[i] = (np.sqrt(1/wmnden[i]) / dens_fix_ma) * 4.2
    elif avgtype == 'rmn':
        rmn[i] = (rmnnum[i]/nspec[i]) / dens_fix_ma
        rmnerr[i] =  ((np.sqrt(rmnden[i])/nspec[i]) / dens_fix_ma) * 4.2

# Adding a couple of extra pieces of useful information to the dictionaries..
nspec = header(nspec)

# Ending the clock!
end_time = time.time()
print('Done! This took '+np.str(np.int(np.ceil(end_time - start_time)))+
    ' seconds!')

# Saving files
np.save('stacks/'+rgn_name+'/'+fol_name+'/wave',10**wave_ma)
np.save('stacks/'+rgn_name+'/'+fol_name+'/nspec',nspec)
np.save('stacks/'+rgn_name+'/'+fol_name+'/resolution',total_reso)
np.save('stacks/'+rgn_name+'/'+fol_name+'/prob_gal',prob_gal)
if avgtype == 'wmn':
    np.save('stacks/'+rgn_name+'/'+fol_name+'/wmn',wmn)
    np.save('stacks/'+rgn_name+'/'+fol_name+'/wmnerr',wmnerr)
elif avgtype == 'rmn':
    np.save('stacks/'+rgn_name+'/'+fol_name+'/rmn',rmn)
    np.save('stacks/'+rgn_name+'/'+fol_name+'/rmnerr',rmnerr)

# Sigma tracking..
if sigma_tracking == 'y':
    np.save('stacks/'+rgn_name+'/'+fol_name+'/ha_sigma',ha_sig)
    np.save('stacks/'+rgn_name+'/'+fol_name+'/sigmaSNR',sig_noi)
    np.save('stacks/'+rgn_name+'/'+fol_name+'/sigmaWeights',sig_weight)

# Printing some stuff for the end..
if verbose == 'y':
    print('')
    print('Stack Overview:')
    print(np.str(n)+' total galaxies. '+np.str(prob_gal.shape[0])+' skipped.')
    if gal_list_tog == 'y':
        print('Stacking region '+np.str(bin_num)+': '+np.str([usr_rgn[bin_num,0],
            usr_rgn[bin_num,1],usr_rgn[bin_num,2],usr_rgn[bin_num,3]]))
    else:
        print('Stacking region: '+rgn_name)
    print('SNR cut-off: '+np.str(nspec['sncut']))
    if 'r_borders' in list(nspec.keys()):
        print('Radial bin edges:')
        print(np.str(nspec['r_borders']))
    if 'i_borders' in list(nspec.keys()):
        print('Inclination bin edges:')
        print(np.str(nspec['i_borders']))
    regi = np.argwhere((10**wave_ma>5500)&(10**wave_ma<7000))
    print('Key'.center(5)+' : # spectra : SNR at r-band')
    keys = sorted(list(total_reso.keys()))
    for i in keys:
        if avgtype == 'wmn':
            flux = wmn[i][regi]
            sig = wmnerr[i][regi]
        elif avgtype == 'rmn':
            flux = rmn[i][regi]
            sig = rmnerr[i][regi]
        snr = flux/sig
        snr = np.int(np.median(snr[np.isfinite(snr)]))
        print(i.center(5)+' : '+np.str(np.int(np.max(nspec[i]))).center(9)+
            ' : '+np.str(snr))
print('')

# Done!
#
##############################################################################