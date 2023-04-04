import matplotlib.pylab as plt

from calc_throughput_tools import compute_throughput, plot_throughput, plot_throughput_peaks

if __name__=='__main__':
	#load inputs
	starname = '10Lac'# notes say tip tilt not working, so dont expect best results
	mode='kpf'
	expected_throughput_file= 'inputs/throughput_models/transmission_kpf.txt'
	kpf_file                = 'inputs/10Lac/KP.20221114.16441.86_L1.fits'
	standard_spec_file      = 'inputs/10Lac/10lac_mod_004.fits' # this might have to be made depending on source of model 
	seeing                  = 1 # take from log

	xdat, measured_throughput, expected_throughput = compute_throughput(kpf_file,standard_spec_file,mode=mode,seeing=seeing,expected_throughput_file=expected_throughput_file)
	
	title = starname + '(seeing=%s) \n(file %s)'%(seeing,kpf_file.split('/')[-1].strip('_L1.fits'))
	plot_throughput(xdat, measured_throughput, expected_throughput,fac=1,title=title)
	plt.savefig('outputs/%s_throughput.png'%starname) # maybe change this to label it based on file used
	plot_throughput_peaks(xdat, measured_throughput, expected_throughput,fac=1,title=title)
	plt.savefig('outputs/%s_throughput_peaks.png'%starname)


