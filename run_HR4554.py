import matplotlib.pylab as plt

from calc_throughput_tools import compute_throughput, plot_throughput, plot_throughput_peaks

if __name__=='__main__':
	#load inputs
	expected_throughput_file='inputs/throughput_models/transmission_kpf.txt'
	kpf_file                = 'inputs/HR4554/KP.20230109.55357.76_L1.fits';note='76'
	standard_spec_file      = 'inputs/HR4554/mhr4554.fits' # this might have to be made depending on source of model 
	seeing                  = 0.58 # take from log

	xdat, measured_throughput, expected_throughput = compute_throughput(kpf_file,standard_spec_file,seeing=seeing,expected_throughput_file=expected_throughput_file)
	
	title = 'HR4554 (seeing=%s) \n(file %s)'%(seeing,kpf_file.strip('_L1.fits')[1:])
	plot_throughput(xdat, measured_throughput, expected_throughput,fac=1,title=title)
	plt.savefig('outputs/HR4554_throughput.png')
	plot_throughput_peaks(xdat, measured_throughput, expected_throughput,fac=1,title=title)
	plt.savefig('outputs/HR4554_throughput_peaks.png')


