#!/usr/bin/env python
from pylab import *
from analysis_library import *
import glob

def single_obs_view(fnm,chain_fnm, nwalkers=100, n_iterations = 1000, n_iteration_filter = 600, true_delay = 1.5, true_delta_mag = 0.3, true_mag = 19., true_sigma = 0.07, true_tau=121.):

	mag1, mag2, time_array, mag1_dat, mag2_dat = read_sim_hubble_light_curve(fnm)

	figure()
	n1 = len(time_array)/2-12
	n2 = len(time_array)/2+12
	subplot(311)
	errorbar(time_array[n1:n2], mag1_dat[n1:n2],fmt='b+', yerr=0.02*ones(len(time_array[n1:n2])), mec='b', ms=3)
	errorbar(time_array[n1:n2], mag2_dat[n1:n2] ,fmt='rx', yerr=0.02*ones(len(time_array[n1:n2])), mec='r', ms=3 )
	ylim(ylim()[1],ylim()[0])
        ma = ylim()[0]
	mb = ylim()[1]
	for k in range(n1,n2):
		plot([time_array[k],time_array[k]] , [ma,mb],'k--')
	errorbar(time_array[n1:n2], mag1_dat[n1:n2], fmt='b+', yerr=0.02*ones(len(time_array[n1:n2])), mec='b', ms=3)
	errorbar(time_array[n1:n2], mag2_dat[n1:n2] ,fmt='rx', yerr=0.02*ones(len(time_array[n1:n2])), mec='r', ms=3 )
	xlabel('time, days')
	ylabel('magnitude')

	subplot(312)
	errorbar(time_array, mag1_dat,fmt='b+', yerr=0.02*ones(len(time_array)), mec='b', ms=3)
	errorbar(time_array, mag2_dat ,fmt='rx', yerr=0.02*ones(len(time_array)), mec='r', ms=3 )
	plot([time_array[n1],time_array[n2]],[ma,ma], 'k-')
	plot([time_array[n1],time_array[n2]],[mb,mb], 'k-')
	plot([time_array[n1],time_array[n1]],[ma,mb], 'k-')
	plot([time_array[n2],time_array[n2]],[ma,mb], 'k-')
	ylim(ylim()[1],ylim()[0])
	xlabel('time, days')
	ylabel('magnitude')

	subplot(313)
	errorbar(time_array, mag1_dat,fmt='b+', yerr=0.02*ones(len(time_array)), mec='b', ms=3)
	errorbar(time_array-true_delay, mag2_dat-true_delta_mag ,fmt='rx', yerr=0.02*ones(len(time_array)), mec='r', ms=3 )
	ylim(ylim()[1],ylim()[0])
	xlabel('time, days')
	ylabel('magnitude')
	subplots_adjust(hspace=0.25, bottom=0.1, top=0.9, left=0.1, right=0.95)
	show()


	samples = get_chain(chain_fnm)
	array_index = range(0,len(samples[:,[0]]))
	delay, d_mag, sigma, tau, avg_mag = get_parameter_samples(samples)

	delay     = filter_samples(delay,   n_iterations, n_iteration_filter)
	d_mag     = filter_samples(d_mag,   n_iterations, n_iteration_filter)
	sigma     = filter_samples(sigma,   n_iterations, n_iteration_filter)
	log10_tau = filter_samples(tau,   n_iterations, n_iteration_filter)
	avg_mag   = filter_samples(avg_mag,   n_iterations, n_iteration_filter)

	figure(figsize=(24,12))
	subplot(521)
	hist(delay, bins=100)
	xlabel('delay, days')
	subplot(522)
	plot(delay, ',')
	ylabel('delay, days')
	xlabel('sample')

	subplot(523)
	hist(d_mag, bins=100)
	xlabel('$\Delta$ mag')
	subplot(524)
	plot(d_mag, ',')
	ylabel('$\Delta$ mag')
	xlabel('sample')

	subplot(525)
	hist(sigma, bins=100)
	xlabel('sigma')
	subplot(526)
	plot(sigma, ',')
	ylabel('sigma')
	xlabel('sample')

	subplot(527)
	hist(log10_tau, bins=100)
	xlabel('log10(tau)')
	subplot(528)
	plot(log10_tau, ',')
	ylabel('log10(tau)')
	xlabel('sample')

	subplot(529)
	hist(avg_mag, bins=100)
	xlabel('avg. mag')
	subplot(5,2,10)
	plot(avg_mag, ',')
	ylabel('avg. mag')
	xlabel('sample')

	dt_min = min(delay)
	dt_max = max(delay)
	delta_dt = (dt_max-dt_min)/1000
	p_dt,  x_dt  = hist_param(delay, arange(dt_min,dt_max,delta_dt))
	dm_min = min(d_mag)
	dm_max = max(d_mag)
	delta_dm = 0.001
	p_dm,  x_dm  = hist_param(d_mag, arange(dm_min,dm_max,delta_dm))

	
	filtered_samples = samples[mod(array_index,n_iterations)>n_iteration_filter]
	fig= triangle.corner(filtered_samples[:,[0,1,4,2,3]], labels=["delay", "$\Delta$m", r"$\langle m \rangle$", "$\sigma$", r"$log_{10}(\tau)$"], truths=[true_delay, true_delta_mag, true_mag, true_sigma, log10(true_tau)], truth_color='red')

	a1 = axes([.55, .83, .4, .14], axisbg='none')
	errorbar(time_array, mag1_dat,fmt='b+', yerr=0.02*ones(len(time_array)), mec='b', ms=3)
	errorbar(time_array, mag2_dat ,fmt='rx', yerr=0.02*ones(len(time_array)), mec='r', ms=3 )
	title('Light Curves')
	#xlabel('time, days')
	ylabel('Magnitude')
	a1.set_xticklabels('')
	mini = min( float(int(min(mag1_dat/0.1)))*0.1, float(int(min(mag2_dat/0.1)))*0.1) 
	maxi = max( float(int(max(mag1_dat/0.1)))*0.1, float(int(max(mag2_dat/0.1)))*0.1)
	print 'mini, max', mini,maxi
	setp(a1, yticks=arange(mini,maxi+0.11, 0.2))
	ylim(maxi+0.1,mini-0.1)

	a2 = axes([.55, .65, .4, .14], axisbg='none')
	dt = x_dt[argmax(p_dt)]
	dm = x_dm[argmax(p_dm)]
	errorbar(time_array, mag1_dat,fmt='b+', yerr=0.02*ones(len(time_array)), mec='b', ms=3, label='light curve 1')
	errorbar(time_array-dt, mag2_dat-dm ,fmt='rx', yerr=0.02*ones(len(time_array)), mec='r', ms=3, label='light curve 2' )
	legend(loc=1)
	xlabel('time, days')
	ylabel('Magnitude')
	mini = min( float(int(min(mag1_dat/0.1)))*0.1, float(int(min((mag2_dat-dm)/0.1)))*0.1 )  
	maxi = max( float(int(max(mag1_dat/0.1)))*0.1, float(int(max((mag2_dat-dm)/0.1)))*0.1 )  
	setp(a2, yticks=arange(mini ,maxi+0.11, 0.1))
	ylim(maxi+0.1,mini-0.05)


if __name__ == "__main__":
    
    import argparse
    import os

    parser=argparse.ArgumentParser(description='view_run routine to view the output of fot_hubble_sim')
    
    parser.add_argument("-fs",  "--data_file_name",        help="file name for the simulated data",     type=str                   )
    parser.add_argument("-fc",  "--chain_file_name",       help="filename for the emcee output chains", type=str                   )
    parser.add_argument("-td",  "--true_delay",            help="true delay",                           type=float, default = 1.5  )
    parser.add_argument("-tdm", "--true_delta_mag",        help="true delta magnitude",                 type=float, default=0.3    )
    parser.add_argument("-tm",  "--true_mag",              help="true magnitude",                       type=float, default = 19.  )
    parser.add_argument("-ts",  "--true_sigma",            help="true sigma",                           type=float, default = 0.07 )
    parser.add_argument("-tt",  "--true_tau",              help="true tau",                             type=float, default = 121. )
    parser.add_argument("-nw",  "--num_walkers",           help="number of walkers",                    type=int,   default = 100  )
    parser.add_argument("-ni",  "--num_iterations",        help="number of iterations",                 type=int,   default = 1000 )
    parser.add_argument("-nif", "--num_iteration_filter",  help="number of iterations filter",          type=int,   default = 600  )
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    #PAESE ARGUMENTS
    args=parser.parse_args()
 
    print args.data_file_name
    print args.chain_file_name
    single_obs_view( args.data_file_name,                 
		     args.chain_file_name,  
		     nwalkers=args.num_walkers, 
		     n_iterations = args.num_iterations, 
		     n_iteration_filter = args.num_iteration_filter,              
	             true_delay = args.true_delay,        
	             true_delta_mag = args.true_delta_mag,
	             true_mag = args.true_mag,            
		     true_sigma = args.true_sigma,        
		     true_tau=args.true_tau) 
    show()  




