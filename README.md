# Schatten_norm_estimation
It gives an estimate of k-Schatten norms of a partially observed symmetric matrix for any integer k upto 7.
It is a MATLAB implementation of an unbiased estimator of k-Schatten norm of a symmetric matrix with missing entries. 
The mathematical details of the project can be found in the paper titled 'Matrix norm estimation from a few entries' co-authored by Ashish Khetan and Professor Sewoong Oh, University of Illinois Urbana Champaign, USA.
This project primarily contains two MATLAB function files 'sch_est_erdos_renyi_sampling.m' and 'sch_est_graph_sampling.m'.
The function 'sch_est_erdos_renyi_sampling.m' takes three inputs: a real symmetric matrix 'M' and a positive integer 'k' less than equal to 7 and a real number 'p' between 0 and 1. Suppose, each entry of 'M' is observed i.i.d. with probability 'prob' from a complete matrix cM, and the missing entries of cM are replaced by zeroes in M. This function returns an estimate of the 'k'-Schatten norm raised to the power 'k' of the complete matrix cM. The estimator is unbiased under Erdos Renyi sampling assumption.
The function 'sch_est_graph_sampling' takes two inputs: a real symmetric matrix 'M' and a positive integer 'k' less than equal to 7. 
Suppose, 'M' is observed from a complete matrix cM, where the missing entries of cM are replaced by zeroes in 'M'. This function returns an estimate of the 'k'-Schatten norm raised to the power 'k' of the complete matrix cM. The estimator is unbiased under the assumption that 'M' is sampled using graph based sampling.
There is a MATLAB script file 'test.m' two test these two functions. It shows how to call the two functions. It generates a random symmetric matrix and randomly sample its entries, and calls the two estimator functions two get an estimate of k-Schatten norms of the matrix.  
For any issues regarding the use of these functions please reach me at khetan2@illinois.edu. I'd be happy to answer any queries and concerns.
