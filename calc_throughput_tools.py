# throughput measuring tools for KPF cap-1h with example usage at bottom
# need to make telluric spectrum at resolution of data
from scipy import interpolate

import numpy as np
import matplotlib.pylab as plt

from astropy.io import fits
from astropy.table import Table


def mod_abmag(filename ='./mhr4544.csv'):
	"""
	open ab mag csv file and resave as fits file in uniots of  egs/cm2/s/A

	checked if conversion worked by opening fhr4544 file which is in uniots
	of egs/cm2/s/A --> it matches

	"""
	f = np.loadtxt(filename).T
	x,y = f[0],f[1]

	# convert ab mag to flux
	#https://en.wikipedia.org/wiki/AB_magnitude
	y_ergdens_hz  = 10**((y+48.6)/-2.5) #erg s−1 cm−2 Hz−1
	y_ergdens_ang = y_ergdens_hz * 3e18 / x**2

	cols=[]
	cols.append(fits.Column(name='WAVELENGTH', format='D', array=x))
	cols.append(fits.Column(name='FLUX', format='D', array=y_ergdens_ang))
	
	primary_hdu = fits.PrimaryHDU()

	tbhdu = fits.BinTableHDU.from_columns(cols)
	hdu = fits.HDUList([primary_hdu,tbhdu])
	hdu.writeto(filename.strip('csv')+'fits',overwrite=True)


def gaussian(x, shift, sig):
    ' Return normalized gaussian with mean shift and var = sig^2 '
    return np.exp(-.5*((x - shift)/sig)**2)/(sig * np.sqrt(2*np.pi))


def define_lsf(v,res):
    """
    define gaussian in pixel elements to convolve resolved spectrum with to get rightish resolution
    """
    dlam  = np.median(v)/res
    fwhm  = dlam/np.mean(np.diff(v)) # desired lambda spacing over current lambda spacing resolved to give sigma in array elements
    sigma = fwhm/2.634 # FWHM is dl/l but feed sigma    
    x = np.arange(sigma*10)
    gaussian = (1./sigma/np.sqrt(2*np.pi)) * np.exp(-0.5*( (x - 0.5*len(x))/sigma)**2 )

    return gaussian

def degrade_spec(x,y,res):
    """
    given wavelength, flux array, and resolving power R, return  spectrum at that R
    """
    lsf      = define_lsf(x,res=res)
    y_lowres = np.convolve(y,lsf,mode='same')

    return y_lowres


def calc_coupling(seeing=0.7,fiber_size=1.14):
	"""
	make moffat function, integrate in 3 dimensions, 
	sum that integral out to KPF fiber bounds
	checked it matches steve's results
	Returns the overlap into the fiber

	inputs
	------
	seeing: [arcsec]
		seeing for the observation
	fiber_size [arcsec]
		fiber size in arcsec

	outputs:
	-------
	inFiber 
		the fraction of light coupled into the fiber
	"""
	#def moffat(theta,seeing):

	#return I_of_theta
	
	theta = np.arange(0,30,0.001)

	HW = seeing/2.0
	num1 = (2**(1/7) - 1) * (7-1) * 0.8
	den1 = (1+(2**(1/7) - 1)*(theta/HW)**2)**7
	num2 = (2**0.5 - 1)*(2-1)*0.2
	den2 = (1+(2**0.5 -1)*(theta/HW)**2)**2
	coeff= (1/np.pi/HW**2)
	I_of_theta = coeff * (num1/den1 + num2/den2)

	#moffat_image = np.zeros((len(theta),len(theta)))
	#for i, th1 in enumerate(theta):
	#	for j, th2 in enumerate(theta):
	#		theta_temp =  np.sqrt(th1**2  + th2**2)
	#		moffat_image[i,j] = moffat(theta_temp,seeing)
	
	Inew      = np.zeros_like(theta)
	for i,th in enumerate(theta):
		Inew[i] = 2*np.pi*i * I_of_theta[i]

	#KPF_fiber_size = 4.4 #arcsec, edge to edge
	isub = np.where(theta < fiber_size/2.)

	inKPF = np.trapz(Inew[isub], theta[isub])
	allLight = np.trapz(Inew,theta)

	inFiber = inKPF/allLight

	return inFiber


def vac_to_stand(wave_vac):
    """Convert vacuum wavelength (Ang) to standard wavelength in air since we're
    doing ground based stuff.

	https://idlastro.gsfc.nasa.gov/ftp/pro/astro/vactoair.pro
    Equation from Prieto 2011 Apogee technical note
    and equation and parametersfrom Cidor 1996

    wavelength units in Angstroms!"""
    # eqn
    sigma2= (1e4/wave_vac)**2.
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) + \
                            1.67917e-3/( 57.362 - sigma2)

    # return l0 / n which equals lamda
    return wave_vac/fact

def load_standard_spec(filename,exptime=20):
	"""
	load flux standard spectrum
	
	inputs:
	------
	filename - str
		filename of standard star data, fits format matching calspec format
		units of model data: erg s-1 cm-2 A-1
		model wave units: A

	exptime - [seconds]
		exposure time to convert from photons/s/A to photons/A
	
	output:
	------
	xmod_new, ymod_phot_ang - wavelength and flux density
	units of output: 
		ymod_phot_ang: photons per angstrom, 
		xmod_new: nanometers
	"""
	fmod = fits.open(filename)
	#fstis = fits.open('10lac_stis_006.fits')
	
	xmod,ymod = fmod[1].data['WAVELENGTH'], fmod[1].data['FLUX']
	isub = np.where(xmod < 9000)[0]
	xmod,ymod = xmod[isub],ymod[isub]

	# resample onto finer and equal grid
	fint = interpolate.interp1d(xmod,ymod,bounds_error=False,fill_value=0)
	xmod_new = 10*np.arange(350,900,0.001) # angstrom
	ymod_new = fint(xmod_new)

	ymod_photdensity = ymod_new * 5.03e7 * xmod_new # now in phot/cm2/A/s - from harvard page ..times lambda
	
	#exptime=20#s , header has 600s but that is wrong
	keck_area = 75.76 #m2
	cm_per_meter = 100 #cm/m

	# modify ymod to phot/A
	ymod_phot_ang = ymod_photdensity * exptime*keck_area*cm_per_meter**2
	xmod_new/=10 # convert to nanometers
	
	return xmod_new, ymod_phot_ang 

def load_em_data(filename='./data/10Lac/KP.20221110.26825.94_L1.fits',ploton=False):
	"""
	"""
	df_SCI_EM = Table.read(filename, format='fits',hdu='EXPMETER_SCI').to_pandas()

	EM_gain = 1.48424 # e-/ADU

	# Define wavelength arrays and disperion at each wavelength (nm per pixel)
	wav_SCI_str = df_SCI_EM.columns[2:] # string (center) wavelengths of each pixel
	wav_SCI     = df_SCI_EM.columns[2:].astype(float) # float (center) wavelengths of each pixel
	disp_SCI = wav_SCI*0+np.gradient(wav_SCI,1)*-1
	disp_SCI_smooth   = np.polyval(np.polyfit(wav_SCI,disp_SCI, deg=6),wav_SCI)
	wav_SCI_smooth    = wav_SCI[0] + -1*np.cumsum(disp_SCI_smooth)
	
	# define normalized flux array (e- / time)
	df_SCI_EM_norm        = df_SCI_EM[wav_SCI_str] * EM_gain

	# define time arrays
	date_beg = np.array(df_SCI_EM["Date-Beg"], dtype=np.datetime64)
	date_end = np.array(df_SCI_EM["Date-End"], dtype=np.datetime64)
	tdur_sec = (date_end-date_beg).astype(float)/1000. # array exposure durations in sec
	time     = (date_beg-date_beg[0])/1000 # array of times since beginning in sec

	int_SCI_spec         = df_SCI_EM_norm.sum(axis=0) / np.sum(tdur_sec) # flux vs. wavelength per sec (use first five samples)
	int_SCI_flux         = df_SCI_EM.sum(axis=1)                         # flux (ADU) vs. time (per sample)

	header = fits.getheader(filename)
	#exptime=header['EXPTIME']
	exptime=1 # data is already converted to per second
	airmass=header['AIRMASS']
	note   =header['OBJECT']

	# reshape to match structure of kpf main data
	return  wav_SCI_smooth.reshape(1,len(wav_SCI)), int_SCI_spec.values.reshape(1,len(wav_SCI)), exptime,airmass, note

def load_kpf_data(filename,ploton=False):
	"""
	"""
	f = fits.open(filename)
	exptime=f[0].header['EXPTIME']
	airmass=f[0].header['AIRMASS']
	XPO = f[0].header['AUTXCENT']
	YPO = f[0].header['AUTYCENT']
	# not correct to sum them if wavelength sol is off for any of the traces but doing it bc wvl sols seems close enough
	greenflux = f['GREEN_SCI_FLUX1'].data + f['GREEN_SCI_FLUX2'].data+f['GREEN_SCI_FLUX3'].data
	redflux   = f['RED_SCI_FLUX1'].data + f['RED_SCI_FLUX2'].data+f['RED_SCI_FLUX3'].data
	redwave   = f['RED_SCI_WAVE2'].data/10.
	greenwave = f['GREEN_SCI_WAVE2'].data/10

	allflux = np.concatenate((greenflux,redflux))
	allwave = np.concatenate((greenwave,redwave))
	if ploton:
		plt.figure()
		for order in np.arange(np.shape(allwave)[0]):
			plt.plot(allwave[order],allflux[order])

	return allwave,allflux,exptime,airmass,XPO,YPO

def load_hk_data(filename,ploton=False):
	"""
	"""
	f = fits.open(filename)
	exptime=f[0].header['EXPTIME']
	airmass=f[0].header['AIRMASS']
	note   =f[0].header['OBJECT']
	# not correct to sum them if wavelength sol is off for any of the traces but doing it bc wvl sols seems close enough
	sci = f['CA_HK_SCI'].data
	sky = f['CA_HK_SKY'].data
	wave = f['CA_HK_SCI_WAVE'].data

	return wave,sky,exptime,airmass,note

def load_telluric(datapath='./inputs/telluric/',l0=380,l1=900,airmass=1,pwv=1):
	"""
	filename is hard coded
	Open psg model for all molecules
	returns wavelength, h2o, co2, ch4, co for l0-l1 range specified

	no o3 here .. make this more flexible
	--------
	"""
	filename = datapath +  'psg_out_2015.06.17_l0_380nm_l1_900nm_res_0.002nm_lon_204.53_lat_19.82_pres_0.5826.fits'
	f        = fits.getdata(filename)
	pwv0     = fits.getheader(filename)['PWV']
	airmass0 = fits.getheader(filename)['airmass']

	x = f['Wave/freq']

	h2o = f['H2O']
	co2 = f['CO2']
	ch4 = f['CH4']
	co  = f['CO']
	o3  = f['O3'] # O3 messes up in PSG at lam<550nm and high resolution bc computationally expensive, so don't use it if l1<550
	n2o = f['N2O']
	o2  = f['O2']
	#ray = f['Rayleigh']

	# load better rayleigh and o3 and interpolate onto other grid
	f = np.loadtxt(datapath +  'maunakea_extinction_2013.txt').T
	tempwl,tempext,std = f[0],f[1],f[2] # wavelength (nm)
	fint_ray = interpolate.interp1d(tempwl,10**((tempext+std)/-2.5),bounds_error=False,fill_value=0)
	ray = fint_ray(x)

	f = np.loadtxt(datapath +  'psg_trn.txt').T
	tempwl,tempo3,tempray = f[0],f[4],f[10]
	fint_o3 = interpolate.interp1d(tempwl,tempo3,bounds_error=False,fill_value=0)
	fint_ray = interpolate.interp1d(tempwl,tempray,bounds_error=False,fill_value=0)
	fac = 1.5 # use this to match PSG rayleigh data to measured mauna kea data
	o3  = fint_o3(x)**fac
	ray = fint_ray(x)**fac

	idelete = np.where(np.diff(x) < 0.0001)[0]  # delete non unique points - though I fixed it in code but seems to pop up still at very high resolutions
	x, h2o, co2, ch4, co, o3, n2o, o2, ray= np.delete(x,idelete),np.delete(h2o,idelete), np.delete(co2,idelete),np.delete(ch4,idelete),np.delete(co,idelete),np.delete(o3,idelete),np.delete(n2o,idelete),np.delete(o2,idelete),np.delete(ray,idelete)
	isub = np.where((x > l0) & (x < l1))[0]

	return x[isub], h2o[isub]**((pwv/pwv0) * (airmass/airmass0)) *  (co2[isub] * ch4[isub]* co[isub]* o3[isub]* n2o[isub]* o2[isub]* ray[isub])**(airmass/airmass0)

def compute_throughput(kpf_file,standard_spec_file,mode='expmeter',seeing=1,expected_throughput_file='../transmission_kpf.txt'):
	"""
	Loads files and computes the throughput for the observation 

	inputs:
	-------
	kpf_file (str)
		filename of KPF file to consider, EM must be L0 file, otherwise give it L1 file
	standard_spec_file (str)
		file  

	per wavelength bin
	fac serves to modulate expected throughput in case want to see the factor needed to match teh calculated
	"""
	if mode=='kpf': xdat,ydat,exptime,airmass,xpo,ypo   = load_kpf_data(kpf_file); R=100000
	if mode=='em':  xdat,ydat,exptime,airmass,note      = load_em_data(kpf_file); R=100
	if mode=='hk':  xdat,ydat,exptime,airmass,note      = load_hk_data(kpf_file); R=15000

	Norders = np.shape(xdat)[0]

	x_star,y_star         = load_standard_spec(standard_spec_file,exptime=exptime)
	y_star_lowres         = degrade_spec(x_star,y_star,R)
	x_telluric,y_telluric = load_telluric(airmass=airmass,pwv=1.5)
	y_telluric_lowres     = degrade_spec(x_telluric,y_telluric,R)
	try:
		kpf_throughput_x, kpf_throughput_y = np.loadtxt(expected_throughput_file,delimiter=',').T
	except:
		kpf_throughput_x, kpf_throughput_y = np.loadtxt(expected_throughput_file).T		

	# make interpolation functions for model components
	star_interp     = interpolate.interp1d(x_star,y_star_lowres,bounds_error=False,fill_value=0)
	telluric_interp = interpolate.interp1d(x_telluric,y_telluric_lowres,bounds_error=False,fill_value=0)
	kpf_throughput_interp  = interpolate.interp1d(kpf_throughput_x, kpf_throughput_y,bounds_error=False,fill_value=0)

	# calc seeing factor to adjust seeing conditions
	seeing_ratio = calc_coupling(seeing=seeing)/calc_coupling(seeing=0.7) # steve assumes 0.7 in throughput model

	# loop through orders to get model photons using flux density
	# calculation is F_ang [phot/ang] * dx [ang] where dx is wavelength element of KPF pixel in real data array
	measured_throughput = np.zeros_like(xdat)
	expected_throughput = np.zeros_like(xdat)
	for order in np.arange(Norders):
		# reinterp model on data to divide to get throughput
		ymod_temp = np.zeros_like(ydat[order])
		for i,xi in enumerate(xdat[order]):
			if xi==xdat[order][-1]: 
				ymod_temp[i] = ymod_temp[i-1]
				continue
			model_photons_per_ang = star_interp(xi)
			ymod_temp[i] = model_photons_per_ang * 10*np.abs((xdat[order][i+1] - xdat[order][i]))
		# store measured throughput as data / model
		measured_throughput[order] = ydat[order]/ymod_temp
		# store expected throughput 
		expected_throughput[order] = seeing_ratio * kpf_throughput_interp(xdat[order])*telluric_interp(xdat[order]) #throughput  model times telluric spec

	return xdat, measured_throughput, expected_throughput

def plot_throughput(xdat, measured_throughput, expected_throughput,fac=1,title=''):
	"""
	"""
	Norders = np.shape(xdat)[0]
	plt.figure(figsize=(8,4))
	for order in np.arange(Norders):
		if order==0:
			plt.plot(xdat[order],100*fac*expected_throughput[order],'k--',label='Expected')
			plt.plot(xdat[order],100*measured_throughput[order],label='Observed')
		else:
			plt.plot(xdat[order],100*fac*expected_throughput[order],'k--')
			plt.plot(xdat[order],100*measured_throughput[order])
		plt.ylabel('Throughput (%)')
		plt.xlabel('Wavelength (nm)')
		plt.title(title)
		plt.legend()
	#plt.savefig('ThroughputModel')


def plot_throughput_peaks(xdat, measured_throughput, expected_throughput,fac=1,title=''):
	"""
	"""
	Norders = np.shape(xdat)[0]
	xmeasured, xexpected, measured, expected = [], [], [], []
	for order in np.arange(Norders):
		xmeasured.append(xdat[order][np.argmax(measured_throughput[order])])
		measured.append(np.max(measured_throughput[order]))
		xexpected.append(xdat[order][np.argmax(expected_throughput[order])])
		expected.append(np.max(expected_throughput[order]))

	plt.figure(figsize=(8,4))
	plt.plot(xmeasured,measured,lw=2,label='Measured')
	plt.plot(xexpected,expected,'k--',label='Expected')
	plt.ylabel('Throughput (%)')
	plt.xlabel('Wavelength (nm)')
	plt.title(title)
	plt.grid()
	plt.legend()
	#plt.savefig('ThroughputModelPeaks')



if __name__=='__main__':
	#mod_abmag(filename ='./mhr4554.csv') # if need to make fits version of model

	#load inputs
	expected_throughput_file='inputs/throughput_models/transmission_kpf.txt'
	kpf_file             = 'inputs/HR4554/KP.20230109.55357.76_L1.fits';note='76'
	standard_spec_file   = 'inputs/HR4554/mhr4554.fits' # this might have to be made depending on source of model 
	seeing               = 0.58 # take from log

	xdat, measured_throughput, expected_throughput = compute_throughput(kpf_file,standard_spec_file,seeing=seeing,expected_throughput_file=expected_throughput_file)
	
	title = 'HR4554 (seeing=%s) \n(file %s)'%(seeing,kpf_file.strip('_L1.fits')[1:])
	plot_throughput(xdat, measured_throughput, expected_throughput,fac=1,title=title)
	plot_throughput_peaks(xdat, measured_throughput, expected_throughput,fac=1,title=title)


