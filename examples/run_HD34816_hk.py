import matplotlib.pylab as plt

from calc_throughput_tools import compute_throughput, plot_throughput, plot_throughput_peaks

if __name__=='__main__':
	#load inputs
	mode='hk'
	starname = 'HD34816'
	expected_throughput_file= 'inputs/throughput_models/hk_transmission_total_042921.txt'
	kpf_file                = 'inputs/HD34816/KP.20230114.28722.45_L1.fits'
	standard_spec_file      = 'inputs/HD34816/lamlep_mod_004.fits' # this might have to be made depending on source of model 
	seeing                  = 0.55 # take from log

	xdat, measured_throughput, expected_throughput = compute_throughput(kpf_file,standard_spec_file,mode=mode,seeing=seeing,expected_throughput_file=expected_throughput_file)
	
	title = starname + '(seeing=%s) \n(file %s)'%(seeing,kpf_file.split('/')[-1].strip('_L1.fits'))
	plot_throughput(xdat, measured_throughput, expected_throughput,fac=1,title=title)
	plt.savefig('outputs/hk_%s_throughput.png'%starname) # maybe change this to label it based on file used


