# MultiSE
Code to design 3D k-space trajectories and optimized excitation pulses. The following MRM paper used this to design the presented waveforms: A multi spin echo pulse sequence with optimized excitation pulses and a 3D cone readout for hyperpolarized 13C imaging. 

Cone trajectory desing:

cone_trajectory_example.m shows the design example used in the paper with 3.2cm FOV and 2mm resolution and writes the resulted wavefroms into files using the VnmrJ format. For generating the trajectory itself coneTrajDesign.m function is called which uses the functions from the paper: 

Gurney PT, Hargreaves BA, Nishimura DG. Design and analysis of a practical 3D cones trajectory. Magnetic Resonance in Medicine 2006;55(3):575-582.

******************************************************************************************************************************************
These functions are not provided, please download them from http://mrsrl.stanford.edu/~ptgurney/cones.php and add to path. Please keep in mind that the wcc.c mex file is set to the gyromagnetic ratio of proton!!! If used for other nuclei replace all the occurences of the number 4258 accordingly then compile the file.
******************************************************************************************************************************************

Once the cone segments are designed they have to be rewinded, i.e returned to the center of the k-space. For this there are two options available. First, gradRewinder3D.m which tries to find a least time solution "analytically". The resulted rewinder is very close to time optimal and fast to compute although in some cases (i.e. for some FOV and resolution pairs) might give unnecessarily long waveform. The waveforms presented in the paper used this function and worked fine for the FOV-resolution pairs but check on the plots all the time!
The second option is based on convex optimizaton and the cvx solver which is easy to set up and free to download from: http://cvxr.com/cvx/download/
Once added to path the last (10th !!!) argument or the last defaul value in the coneTrajDesign.m function should be changed to 'true'. This believed to give time-optimal result although takes more time to compute.


RF pulse design:

fminsearch_optim_pulseDesign.m uses a conventional SLR pulse designed by means of the auxiliary functions gen_alpha.m and iSLR.m and fine-tunes it with numerical optimization in the fmin_pulse_design_cost.m function. For evaluating the cost function the Bloch-simulator from http://www-mrsrl.stanford.edu/~brian/blochsim/ is used. Once freely downloaded add it to path and change the gyromagnetic ratio according to the nucleus of interest then compile. 

******************************************************************************************************************************************
Due to computational cost this method is not recommended for routine pulse design, rather for specific problems only!!!
******************************************************************************************************************************************

No responsibility taken.

Vencel Somai 2020. -> vs460@cam.ac.uk
