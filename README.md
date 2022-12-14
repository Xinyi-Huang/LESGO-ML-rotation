By Xinyi Huang on Jun. 30th 2022

Modified version of LESGO for rotating channel using hard-coded neural network.
The original paper by Xinyi Huang can be found in https://doi.org/10.1063/1.5129178.
See https://lesgo.me.jhu.edu/ for some instructions. 
See README-JHU file in code folder for more instructions.

To compile:
 - Make sure the modules are correctly loaded. On ACI system, it is 'module purge; module load intel impi cmake'.
 - Make sure the previous build folder is removed. 'rm -r bld'.
 - Compile by using './build-lesgo' and generate the binary file 'lesgo-mpi'.

To run:
 - Use the binary 'lesgo-mpi', the config file 'lesgo.conf'; load correct modules and submit the pbs file, which is 'job.pbs' on ACI system.
 - Change the configuration according to needs. Make sure the number of cores in 'mpirun' command is the same as the number of cores in 'lesgo.conf' file.

To postprocess:
 - A test postprocessing code 'test_read.m' is attached. The function 'read_snap.m' is the only file you need for getting instantaneous flow fields.
