%% Author: Ashish Khetan, email: khetan2@illinois.edu
%This function returns weighted counts of all pseudographs of length k, for
%any k upto 7, in a given symmetric matrix M. 
function [psuedograph] = weighted_pseudograph_counts(M,k)
n = size(M,2);
O = M - diag(diag(M)); % M without diagonal entries
D = diag(diag(M)); % diagonal entries of M
R = ones(n) - diag(ones(n,1)); % matrix of all-ones except on diagonals
switch k
    case 1
        psuedograph = sum(sum(D,2));
    case 2
        psuedograph = zeros(2,1);
        psuedograph(1) = sum(sum(O.^2,2));
        psuedograph(2) = sum(sum(D.^2,2));
    case 3
        psuedograph = zeros(3,1);
        psuedograph(1) = trace(D^3);
        psuedograph(2) = trace(3*D*O^2);
        psuedograph(3) = trace(O^3);
    case 4
        psuedograph = zeros(7,1);
        psuedograph(1) = sum(sum(D.^4,2));
        psuedograph(2) = sum(sum(O.^4,2));
        psuedograph(3) = 4*trace(O*O*D*D);
        psuedograph(4) = 2*sum(sum((O.^2*O.^2).*R,2));
        psuedograph(5) = 2*trace(O*D*O*D);
        sch_all = trace(M*M*M*M) - sum(psuedograph(1:5));
        psuedograph(6) = trace(O*O*O*O) - sum(psuedograph([2,4])); 
        psuedograph(7) = sch_all - psuedograph(6); 
    case 5
        psuedograph = zeros(12,1);
        psuedograph(1) = trace(D.^5);
        psuedograph(2) = 5*sum(sum(D*O.^4,2));
        psuedograph(3) = 5*sum(sum(D.^3*O.^2,2));
        psuedograph(4) = 3*(10/6)*trace(O.^3*O*O);
        psuedograph(5) = 5*sum(sum(D*O.^2*D.^2,2));
        psuedograph(6) = 5*sum(sum((O.^2*D*O.^2).*R,2));
        psuedograph(7) = 5*sum(sum((D*O.^2*O.^2).*R,2));
        psuedograph(8) = 5*trace(O*O*O*D.^2);
        psuedograph(9) = 5*(sum(diag(O^3).*sum(O.^2,2))-2*trace(O.^3*O*O));
        sch_all = trace(M^5) - sum(psuedograph(1:9));
        psuedograph(10) = trace(O*O*O*O*O) - sum(psuedograph([4,9])); %without diagonals
        sch_all_diag = sch_all - psuedograph(10); %with diagonals
        psuedograph(11) = 5*trace(O*D*O*D*O);
        psuedograph(12) = sch_all_diag - psuedograph(11);
    case 6
        psuedograph = zeros(32,1);
        psuedograph(1) = sum(sum(D.^6,2));
        psuedograph(2) = sum(sum(O.^6,2));
        psuedograph(3) = 6*sum(sum((O.^2*O.^4).*R,2));
        psuedograph(4) = 6*sum(sum((O.^2*D.^4).*R,2));
        psuedograph(5) = 9*sum(sum((D.^2*O.^4).*R,2));
        psuedograph(6) = 3*sum(sum((D.^2*O.^2*D.^2).*R,2));
        psuedograph(7) = 6*sum(sum((D.^2*O.^2*O.^2).*R,2));
        psuedograph(8) = 9*sum(sum((O.^2*D.^2*O.^2).*R,2));
        psuedograph(9) = 6*sum(sum((D.^3*O.^2*D).*R,2));
        psuedograph(10) = 6*sum(sum((D*O.^4*D).*R,2));
        psuedograph(11) = 3*sum(sum((O.^2*O.^2).*R,2).*sum(O.^2,2) - sum((O.^4*O.^2).*R,2) - diag((O.^2*O.^2*O.^2)));
        psuedograph(12) = 4*trace((O.^2*O.^2*O.^2));
        psuedograph(13) = 2*sum(sum(O.^2,2).^3 - sum(O.^6,2) - 3*(sum(O.^4,2).*sum(O.^2,2)-sum(O.^6,2)));
        psuedograph(14) = 3*sum(sum((D*O.^2*O.^2*D).*R,2));
        psuedograph(15) = 12*sum(sum((D*O.^2*D*O.^2).*R,2));
        psuedograph(16) = 6*sum(sum((O.^3*O).*R.*(O*O),2) - sum((O.^4*O.^2).*R,2));
        psuedograph(17) = 6*trace(D.^3*O*O*O);
        psuedograph(18) = 24*trace(D*O.^3*O*O);
        psuedograph(19) = 6*trace(D*O*O.^3*O);
        psuedograph(20) = 6*(sum(sum(((O^2.*R).*((O*D.^2*O).*R)),2)) - sum(sum((O.^2*D.^2*O.^2).*R),2));
        psuedograph(21) = 12*trace(O*D.^2*O*D*O);
        psuedograph(22) = 6*(sum(sum(((O*O).*R.*(O*O) - (O.^2*O.^2).*R),2).*sum(O.^2,2))...
            - 2*sum(sum(((O.^3*O).*R.*(O*O) - (O.^4*O.^2).*R),2))...
            - sum(sum(((O*O).*R.*(O*O) - (O.^2*O.^2).*R).*(O.^2),2)));
        psuedograph(23) = 9*sum(sum(((O*O).*R.*(O*O) - (O.^2*O.^2).*R).*(O.^2),2)); 
        psuedograph(24) = 12*sum(diag(O*D*O*O).*sum(O.^2,2) - diag(O.^3*D*O*O) - diag(O.^3*O*D*O));
        psuedograph(25) = 6*sum(diag(O*O*O).*sum(O.^2*D,2) - 2*diag(O.^3*D*O*O));
        psuedograph(26) = 12*sum(diag(O*O*O).*diag(D).*sum(O.^2,2) - 2*diag(O.^3*O*O).*diag(D));
        sch_all = trace(M^6) - sum(psuedograph(1:26));
        sch_all_wdiag = trace(O*O*O*O*O*O) - sum(psuedograph([2,3,11,12,13,16,22,23])); %without diagonals
        psuedograph(27) = 3*sum(diag(O*O*O).^2 - 2*diag(O.^2*O.^2*O.^2)) - (48/36)*psuedograph(23); 
        psuedograph(28) = sch_all_wdiag - psuedograph(27);
        sch_all_diag = sch_all - sch_all_wdiag; %with diagonals
        psuedograph(29) = 2*trace(D*O*D*O*D*O);
        psuedograph(30) = 3*sum(sum((O*D*O).*R.*(O*D*O),2) - sum((O.^2*D.^2*O.^2).*R,2));
        psuedograph(31) = 6*sum(sum((O*D*O*D).*R.*(O*O),2) - sum((O.^2*D*O.^2*D).*R,2));
        psuedograph(32) = sch_all_diag - psuedograph(29) - psuedograph(30) - psuedograph(31);
    case 7
        psuedograph = zeros(69,1);
        psuedograph(1) = sum(diag(D.^7));
        psuedograph(2) = 7*sum(sum(O.^2*D.^5,2)) ;
        psuedograph(3) = 7*sum(sum((D.^2*O.^2*D.^3).*R,2)) ;        
        psuedograph(4) = 14*sum(sum(O.^4*D.^3,2));
        psuedograph(5) = 7*sum(sum(O.^6*D.^1,2));
        psuedograph(6) = 7*sum(sum((D*O.^2*D.^4).*R,2));
        psuedograph(7) = 21*sum(sum((D*O.^4*D.^2).*R,2));
        psuedograph(8) = 7*sum(sum((O.^2*O.^2*D.^3).*R,2));
        psuedograph(9) = 14*sum(sum((O.^2*D.^3*O.^2).*R,2));
        psuedograph(10) = 7*sum(sum((O.^4*O.^2*D).*R,2));
        psuedograph(11) = 21*sum(sum((O.^4*D*O.^2).*R,2));
        psuedograph(12) = 14*sum(sum((D*O.^4*O.^2).*R,2));
        psuedograph(13) = 7*trace(O.^5*O*O);  
        psuedograph(14) = 14*trace(O.^3*O*O.^3); 
        psuedograph(15) = 7*sum(sum((O.^2*O.^2).*R,2).*sum(O.^2*D,2) - sum((O.^4*D*O.^2).*R,2) - diag((O.^2*D*O.^2*O.^2)));
        psuedograph(16) = 14*sum((sum((O.^2*O.^2).*R,2).*sum(O.^2,2) - sum((O.^4*O.^2).*R,2) - diag((O.^2*O.^2*O.^2))).*diag(D));
        psuedograph(17) = 7*sum((sum(O.^2,2).^3 - sum(O.^6,2) - 3*(sum(O.^4,2).*sum(O.^2,2)-sum(O.^6,2))).*diag(D));
        Z1 = 0.5*(sum(O.^2,2).^2 - sum(O.^4,2)); % sum of pairwise product of square of edges
        psuedograph(18) = 14*sum(sum(O.^2*D,2).*Z1 - sum(O.^4*D,2).*sum(O.^2,2) + sum(O.^6*D,2));
        psuedograph(19) = 28*sum(diag(O.^2*O.^2*O.^2).*diag(D));
        psuedograph(20) = 21*sum(sum((D*O.^2*D.^2*O.^2).*R,2));
        psuedograph(21) = 14*sum(sum((D.^2*O.^2*D*O.^2).*R,2));
        psuedograph(22) = 7*sum(sum((D*O.^2*O.^2*D.^2).*R,2));
        psuedograph(23) = 7*sum(diag(O*O*O).*diag(D.^4));
        psuedograph(24) = 28*sum(diag(O.^3*O*O).*sum(O.^2,2) - diag(O.^5*O*O) -diag(O.^3*O*O.^3));
        psuedograph(25) = 7*sum(diag(O*O.^3*O).*sum(O.^2,2) - 2*diag(O.^3*O.^3*O));
        psuedograph(26) = 7*sum(diag(O*O.^3*O).*diag(D.^2));
        psuedograph(27) = 42*sum(diag(O.^3*O*O).*diag(D.^2));
        psuedograph(28) = 7*sum(diag(O*O*O).*sum(O.^4,2) - 2*diag(O.^5*O*O));
        psuedograph(29) = 7*sum(sum((D*O.^2*D*O.^2*D).*R,2));
        psuedograph(30) = 28*sum(diag(O*D*O.^3*O).*diag(D));
        psuedograph(31) = 28*trace(O*D*O.^3*D*O);
        psuedograph(32) = 14*sum(diag(O*D.^2*O*O).*sum(O.^2,2) - diag(O.^3*O*D.^2*O) - diag(O.^3*D.^2*O*O));
        psuedograph(33) = 14*sum(diag(O*D*O*O).*diag(D.^3));
        psuedograph(34) = 7*trace(O*D.^2*O*D.^2*O);
        psuedograph(35) = 7*(sum(sum(((O^2.*R).*((O*D.^3*O).*R)),2)) - sum(sum((O.^2*D.^3*O.^2).*R),2)); 
        psuedograph(36) = 14*sum(sum((O.^3*O).*R.*(O*D*O),2) - sum((O.^4*D*O.^2).*R,2));
        psuedograph(37) = 28*sum(sum((O.^3*D*O).*R.*(O*O),2) - sum((O.^4*D*O.^2).*R,2));
        Z2 = (((O*O).*R)*O - O.*(ones(n,1)*sum(O.^2,1)- O.^2)).*R; % paths of 3 unique edges
        Z3 = (O.*((O*O).*R)).*R; % triangles between two nodes
        Z4 = (O.*((O.^5*O).*R)).*R;
        Z6 = (O.^3.*((O*O).*R)).*R; Z7 = (O.*((O.^3*O.^3).*R)).*R;
        psuedograph(38) = 7*sum(sum(((O.^3*O).*R.*Z2 - ((O.^4*Z3).*R - Z4) - ((Z6*O.^2).*R - Z7)),2)); 
        Z7 = 0.5*sum(sum(O.*((O.^2*O.^2).*R).*((O*O).*R) - O.*((O.^3*O.^3).*R)),2);
        psuedograph(39) = 7*(sum(sum((O.*((O*O).*R).*(sum(O.^2,2)*ones(1,n) - O.^2).*(ones(n,1)*sum(O.^2,1) - O.^2)),2))...
            - sum(sum((O.*((O.^3*O).*R).*(ones(n,1)*sum(O.^2,1) - O.^2)),2))...
            - sum(sum((O.*((O*O.^3).*R).*(sum(O.^2,2)*ones(1,n) - O.^2)),2))...
            + sum(sum((O.*((O.^3*O.^3).*R)),2))) -14*Z7;
        psuedograph(40) = 21*sum(diag(D.^2*O*O*O).*sum(O.^2,2) - 2*diag(D.^2*O.^3*O*O));
        psuedograph(41) = 7*sum(diag(O*O*O).*sum(O.^2*D.^2,2) - 2*diag(O.^3*D.^2*O*O));
        psuedograph(42) = 7*(sum(diag(O*O*O).*sum((O.^2*O.^2).*R,2) - 2*diag(O.^3*O.^3*O))...
            - 2*sum(diag(O.^3*O*O).*sum(O.^2,2) - diag(O.^5*O*O) - diag(O.^3*O*O.^3))) - 28*Z7;
        psuedograph(43) = 14*sum(diag(O*O*O).*Z1 - 2*(diag(O.^3*O*O).*sum(O.^2,2) - diag(O.^5*O*O) - 0.5*diag(O.^3*O*O.^3))); 
        psuedograph(44) = 56*Z7;
        Z8 = (O.*((O.^3*O).*R)).*R; Z9 = (O.^1.*((O*O).*R)).*R; Z10 = (O.*((O.^1*O.^3).*R)).*R;
        Z11 = ((O*O).*R.*Z2 - ((O.^2*Z3).*R - Z8) - ((Z9*O.^2).*R - Z10)); % pentagons between two edges
        psuedograph(45) = 14*(sum(0.5*sum(Z11,2).*sum(O.^2,2)) - (2/14)*psuedograph(38) - sum(sum((O.^2).*Z11,2))); 
        psuedograph(46) = 21*sum(sum((O.^2).*Z11,2));
        psuedograph(47) = 7*sum(sum(Z11,2).*diag(D.^2));
        psuedograph(48) = 7*trace(D.^2*O*D*O*D*O); 
        psuedograph(49) = 14*sum(diag(D*O*O*O).*sum(O.^2*D,2) - 2*diag(D*O.^3*D*O*O));
        psuedograph(50) = 14*sum(diag(O*O*D*O).*sum(O.^2*D,2) - diag(O.^3*D*O*D*O) - diag(O.^3*D.^2*O*O));
        psuedograph(51) = 28*sum(diag(D*O*D*O*O).*sum(O.^2,2) - diag(D*O.^3*D*O*O) - diag(D*O.^3*O*D*O));
        psuedograph(52) = 7*sum(diag(O*D*O*D*O).*sum(O.^2,2) - 2*diag(O.^3*D*O*D*O));
        psuedograph(53) = 14*sum((sum((((O^2.*R).*((O*D*O).*R)) - (O.^2*D*O.^2).*R),2)).*diag(D.^2)); 
        psuedograph(54) = 7*sum(sum(((((O*D*O).*R).*((O*D.^2*O).*R)) - (O.^2*D.^3*O.^2).*R),2));
        Z12 = sum(0.5*sum((((O^2.*R).*(O^2.*R))- ((O.^2*O.^2).*R)).*(O.^2*D),2));
        Z13 = sum(sum(((((O.^3*D*O).*R).*((O*O).*R)) - (O.^4*D*O.^2).*R),2)); 
        Z14 = 0.5*sum(sum(((((O*D*O).*R).*(O^2.*R))- ((O.^2*D*O.^2).*R)).*(O.^2),2));
        Z15 = sum(sum(((((O.^3*O).*R).*((O*D*O).*R)) -(O.^4*D*O.^2).*R),2)); 
        psuedograph(55) = 14*(sum(0.5*(sum((((O^2.*R).*((O^2).*R)) - (O.^2*O.^2).*R),2)).*sum(O.^2*D,2)) -Z13 - Z12);
        psuedograph(56) = 28*(sum(0.5*(sum((((O^2.*R).*((O^2).*R)) - (O.^2*O.^2).*R),2)).*sum(D*O.^2,2)) -Z13 - Z12);
        psuedograph(57) = 14*(sum((sum(((((O*D*O).*R).*((O^2).*R)) - (O.^2*D*O.^2).*R),2)).*sum(O.^2,2))...
            -Z13 - Z15 -2*Z14);
        psuedograph(58) = 14*(sum(0.5*sum(((((O^2.*R).*(O^2.*R))- ((O.^2*O.^2).*R))*D),2).*sum(O.^2,2)) - Z15 - Z12);
        psuedograph(59) = 84*Z12;
        psuedograph(60) = 42*Z14;
        sch_all = trace(M^7) - sum(psuedograph(1:60));
        sch_all_wdiag = trace(O^7) - sum(psuedograph([13,14,24,25,28,38,39,42:46])); %without diagonals
        Z16 = (1/6)*((O*O.*R).^3 - (O.^3*O.^3.*R) -3*((O.^2*O.^2.*R).*(O*O.*R) - (O.^3*O.^3.*R)));
        psuedograph(61) = 42*sum(sum(Z16.*O,2));
        Z17 = sum(sum(0.5*((O*O.*R).^2 - (O.^2*O.^2.*R)),2).*(0.5*diag(O*O*O)));
        psuedograph(62) = 28*(Z17 - (6/84)*psuedograph(61) - (2/42)*psuedograph(46) - (3/56)*psuedograph(44));
        psuedograph(63) = sch_all_wdiag - psuedograph(61) - psuedograph(62);
        sch_all_diag = sch_all - sch_all_wdiag; %with diagonals
        psuedograph(64) = 7*sum(sum((D*O*D*O*D.*R).*(O*O.*R),2)) - 7*sum(sum(D*O.^2*D*O.^2*D.*R,2));
        psuedograph(65) = 7*sum(sum(D*Z11*D,2));
        Z18 = sum(((O*O).*R.*(O*O) - (O.^2*O.^2).*R).*(O.^2),2);
        psuedograph(66) = 7*sum((diag(O*O*O).^2 - 2*diag(O.^2*O.^2*O.^2) - 4*Z18).*diag(D));
        Z20 = 0.5*sum(sum(((O*O.*R).*(O*D*O.*R) - (O.^2*D*O.^2)).*O.^2,2));
        psuedograph(67) = 14*(sum(diag(O*O*O).*diag(O*O*D*O) - 2*diag(O.^2*O.^2*D*O.^2)) - 2*sum(Z18.*diag(D)) - 4*Z20);
        Z21 = (((O*D*O*D).*R)*O - O.*(ones(n,1)*sum(D*O.^2*D,1)- D*O.^2*D)).*R; % paths of 3 unique edges with two diags
        Z22 = (O.*((O*D*O).*R)).*R; % triangles between two nodes with a diag
        Z23 = (O.*((D*O.^3*D*O).*R)).*R;  Z24 = (O.*((O*D*O.^3*D).*R)).*R;
        psuedograph(68) = 7*sum(sum(((O*O).*R.*Z21 - ((O.^2*D*Z22).*R - Z23) - ((Z22*D*O.^2).*R - Z24)),2)); % pentagons between two vertices
        psuedograph(69) = sch_all_diag - sum(psuedograph(64:68));
end


