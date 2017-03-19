%% Author: Ashish Khetan, email: khetan2@illinois.edu
%This is a test file to run two functions: (a) sch_est_erdos_renyi_sampling, (b) sch_est_graph_sampling.
%These functions return an estimate of k-Schatten norm of a partially observed matrix. 
%For test purpose a random symmetric matrix is generated and each entry of
%it is sampled i.i.d. with probability 'p'. 

clc; clear;
d=100; %size of matrix M
r = 20; %rank of matrix M
k = 7; %estimates first k Schatten norms.

%% Erdos Renyi sampling of entries of matrix M
p = 1; % sampling probability

eigvals = unifrnd(1,2,r,1); %generating eigenvalues of M randomly
sigma = diag(eigvals);
U = unifrnd(-1,1,d,r); U = orth(U); % generating eigenvectors of matrix M
M = U*sigma*U'; %constructing matrix M

oM = M.*(rand(d) <= p); %randomly sampling entries of M
oM = (triu(oM) + transpose(triu(oM,1))); %symmeterizing the sampling
Tr = zeros(1,k); esTr_ER = Tr; esTr_GR = Tr;
for sch = 1:k
    esTr_ER(sch) = (sch_est_erdos_renyi_sampling(oM,p,sch)); %computing an estimate of Schatten norm assuming Erdos Renyi sampling
    esTr_GR(sch) = (sch_est_graph_sampling(oM,sch));  % computing an estimate of Schatten norm assuming graph pattern based sampling
    Tr(sch) = (sum(eigvals.^sch)); %computing true Schatten norms
    display(['true ',num2str(sch),' Schatten norm: ',num2str(Tr(sch))]);
    display(['estimated ',num2str(sch),' Schatten norm assuming Erdos Renyi sampling: ',num2str(esTr_ER(sch))])
    display(['estimated ',num2str(sch),' Schatten norm assuming graph based sampling: ',num2str(esTr_GR(sch))])
    fprintf('\n')
end
