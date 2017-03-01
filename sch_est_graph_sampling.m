%% Author: Ashish Khetan, email: khetan2@illinois.edu
%This function takes two inputs: a real symmetric matrix 'M' and a positive integer 'k' less than equal to 7. 
%Suppose, 'M' is observed from a complete matrix cM, where the missing entries of cM are replaced by zeroes in 'M'. 
%This function returns an estimate of the 'k'-Schatten norm raised to the power 'k' of the complete matrix cM that is 
%(\sum_{i=1}^d {\sigma_i^k})$, where $\sigma_i$'s are eigenvalues of the complete matrix cM and d is its dimension.
%The estimator is unbiased under the assumption that 'M' is sampled using graph based sampling.

function [sch_est] = sch_est_graph_sampling(M,k)
M_index = M; M_index(M~=0) = 1; % matrix represing indices of cM that have been observed 

t1 = weighted_pseudograph_counts(M,k); %computing weighted count of all pseudographs of length k in the partially observed matrix M
t2 = weighted_pseudograph_counts(M_index,k); %computing unweighted count of all pseudographs of length k in the partially observed matrix M
t3 = weighted_pseudograph_counts(ones(size(M,1)),k); %computing unweighted count of all pseudographs of length k in a complete graph over d vertices

sch_est = nansum(t1.*(t3./t2)); %computing the unbiased estimate of k-schatten norm assuming graph sampling
end