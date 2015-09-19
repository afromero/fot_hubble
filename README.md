# fot_hubble
A Bayesian inference based quasar light curve delay reconstruction toolkit.

Quick Start

1. Run the gen\_hubble\_runs.py script with a directory specified for your outputs to generate a list of runs. This will generate a text file of commands called run_list.txt
   > ./gen\_hubble\_runs.py -od [your output directory]

3. Execute the lines run_list.txt in the command line. A script that uses pythons multiprocessing library is available. To run it
   > ./parallel\_runs.py -np [number of active processes] -ni [set niceness (optional)]

4. To view the results of any individual run use view_run.py, you need to specify the simulated data file name and emcee chain output file name. 
   The code allows you to specify non-default true values and emcee chain parameters. See the python script for details
   > ./view\_run.py -fs simulated\_data\_file_name.npz -fc emcee\_chain\_output\_file\_name.npz


