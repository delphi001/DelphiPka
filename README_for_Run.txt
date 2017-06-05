
This is a brief intro to run DelphiPKA on Palmetto

You need to load the following modules on Palmetto first:
gcc/4.8.1
openmpi/1.6.4
gsl/1.16


To run the program, you need:
1. PDB file
2. run.prm file (A sample can be found in current directory)
3. PBS script (A sample can be found in current directory)


Steps:
1. Copy the run.prm and PBS script with your PDB file in any working directory.

2. Open the run.prm and edit pdb name entry, charge and radius parameter entries with your desire. Currently supports amber, charmm22 and parse parameters. 

3. For other entries, ref the manual and edit with your own desire.

4. Change the sample.pbs to your desired job name. And qsub the pbs script.






