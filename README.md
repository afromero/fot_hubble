# fot_hubble
A Bayesian inference based quasar light curve delay reconstruction toolkit.

Quick Start

0. Set up the FOTDIR environment variable. This should be the directory of your fot installation.

1. Under FOTDIR, create a directory called outputs

2. Run the gen_hubble_runs.py script to generate a list of runs. This will generate a text file of commands called run_list.txt

3. Execute the runs in the run_list.txt

4. To view the results of any individual run use view_run.py, you need to specify the simulated data file name and emcee chain output file name. 
   The code allows you to specify non-default true values and emcee chain parameters. See the python script for details
   > ./view_run.py -fs simulated_data_file_name.npz -fc emcee_chain_output_file_name.nps


