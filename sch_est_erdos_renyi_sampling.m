%% Author: Ashish Khetan, email: khetan2@illinois.edu
%This function takes three inputs: a real symmetric matrix 'M' and a positive
%integer 'k' less than equal to 7 and a real number 'prob' between 0 and 1.
%Suppose, each entry of 'M' is observed i.i.d. with probability 'prob' from the complete matrix cM, 
%and the missing entries of cM are replaced by zeroes in M.
%This function returns an estimate of the 'k'-Schatten norm raised to the power 'k' of the complete matrix cM that is
%(\sum_{i=1}^d {\sigma_i^k})$, where $\sigma_i$'s are eigenvalues of the complete matrix cM and d is its dimension. 

function [sch_est] = sch_est_erdos_renyi_sampling(M,prob,k)
counts = unique_edge_counts(k); %computing counts of unique edges in all pseudographs of length k
sch_est = nansum(weighted_pseudograph_counts(M,k)./((prob*(ones(length(counts),1))).^counts)); %computing the unbiased estimate of k-schatten norm
end