An Auxilliary Particle Filter (APF) implementation for state and parameter estimation of non-linear and multi-modal dynamic systems with added Non-Gaussian noise.

To get more context as to what APFs are, please refer to Dr. McNames' YouTube channel: https://www.youtube.com/@JamesMcNames/videos

I have not made some simplifications as shown in Dr. McNames' videos because of how the optimality of the algorithm would have suffered. For example, integrations are needed at times in the algorithm which are usually simplified in the literature on APFs but I felt the need to carry them out numerically instead of analytical approximations. 
Parallel Pooling is needed to run this code on MATLAB for efficiency, therefore, the Parallel Computing License is needed for that. If you do not have it, you can revert the "parfor"s to standard for-loops.

If you have any issues with the codes or questions, please contact me at mohammadshedayati@gmail.com or hedayat2@uwindsor.ca
