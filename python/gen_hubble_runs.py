#!/usr/bin/env python
'''
2014 April 6
ARW & LAM, JPL/Caltech
Large simulation run for FULL OF TIME
'''
from pylab import *
import time
#SCRIPT TO PRODUCE THE HUBLLE ORBIT SIMULATIONS AND ANALYSIS RUNS.

if __name__ == "__main__":
    
	import argparse
	parser=argparse.ArgumentParser(description='fot_hubble_sim routine to simulat hubble observations and calculate delay inference')

	parser.add_argument("-od", "--outputdir", help="output directory for simulation products", type=str)
	args=parser.parse_args()

    # ensure the output director ends with '/' to keep the format consisent.
	if(args.outputdir[-1]!='/'):
		args.outputdir = args.outputdir+'/'

	#simulation parameters
	delay=1.5
	delay_prior=1.5
	delay_prior_min=-3.
	delay_prior_max = 3. 
	delta_mag = 0.3 
	delta_mag_prior =0.3 
	delta_mag_prior_min=-5. 
	delta_mag_prior_max=5.
	sigma = 0.07 
	sigma_prior = 0.07
	sigma_prior_min = 0.0007 
	sigma_prior_max = 7. 
	tau = 121 
	tau_prior = 121 
	tau_prior_min = 5. 
	tau_prior_max = 1000.
	avg_mag = 19. 
	avg_mag_prior = 19. 
	avg_mag_prior_min=10.
	avg_mag_prior_max=100. 
	redshift=0.658 
	photometric_uncertainty=0.02 

	fout = open('run_list.txt', 'w')
	num_orbits_list = [40,50,60,70,80,90,100]
	photom_unc_list = [0.02, 0.05, 0.10]
	num_runs =100
	for num_orbits in num_orbits_list:
  		for photometric_uncertainty in photom_unc_list:
			for k in range(0,num_runs):
				output_tag = 'hubble_sim_%d_orbits_%1.2f_pu_r%d'%(num_orbits,photometric_uncertainty,k)
				output_tag = output_tag.replace('.','p')
				log_fname = '%s.log'%(output_tag)
				com = './fot_hubble_sim.py -dt %1.5f -dtp %1.5f -dtpmin %1.5f -dtpmax %1.5f -dm %1.5f -dmp %1.5f -dmpmin %1.5f -dmpmax %1.5f -s %1.5f -sp %1.5f -spmin %1.5f -spmax %1.5f -t %1.5f -tp %1.5f -tpmin %1.5f -tpmax %1.5f -m %1.5f -mp %1.5f -mpmin %1.5f -mpmax %1.5f -z %1.5f -p %1.5f -n %d -o %s -od %s > %s%s'%(delay, delay_prior, delay_prior_min, delay_prior_max, delta_mag, delta_mag_prior, delta_mag_prior_min, delta_mag_prior_max, sigma, sigma_prior, sigma_prior_min, sigma_prior_max, tau, tau_prior, tau_prior_min, tau_prior_max, avg_mag, avg_mag_prior, avg_mag_prior_min, avg_mag_prior_max, redshift, photometric_uncertainty, num_orbits, output_tag, args.outputdir, args.outputdir, log_fname)
				fout.write(com+'\n')	
	fout.close()

