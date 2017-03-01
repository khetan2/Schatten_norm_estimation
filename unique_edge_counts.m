%% Author: Ashish Khetan, email: khetan2@illinois.edu
%This function returns counts of unique edges in all pseudographs of length k, for any k upto 7. 
function [psuedograph] = unique_edge_counts(k)
switch k
    case 1
        psuedograph = 1;
    case 2
        psuedograph = zeros(2,1);
        psuedograph(1) = 1;
        psuedograph(2) = 1;
    case 3
        psuedograph = zeros(3,1);
        psuedograph(1) = 1;
        psuedograph(2) = 2;
        psuedograph(3) = 3;
    case 4
        psuedograph = zeros(7,1);
        psuedograph(1) = 1;
        psuedograph(2) = 1;
        psuedograph(3) = 2;
        psuedograph(4) = 2;
        psuedograph(5) = 3;
        psuedograph(6) = 4;
        psuedograph(7) = 4;
    case 5
        psuedograph = zeros(12,1);
        psuedograph(1) = 1;
        psuedograph(2) = 2;
        psuedograph(3) = 2;
        psuedograph(4) = 3;
        psuedograph(5) = 3;
        psuedograph(6) = 3;
        psuedograph(7) = 3;
        psuedograph(8) = 4;
        psuedograph(9) = 4;
        psuedograph(10) = 5; 
        psuedograph(11) = 5;
        psuedograph(12) = 5;
    case 6
        psuedograph = zeros(32,1);
        psuedograph(1) = 1;
        psuedograph(2) = 1;
        psuedograph(3) = 2;
        psuedograph(4) = 2;
        psuedograph(5) = 2;
        psuedograph(6) = 3;
        psuedograph(7) = 3;
        psuedograph(8) = 3;
        psuedograph(9) = 3;
        psuedograph(10) = 3;
        psuedograph(11) = 3;
        psuedograph(12) = 3;
        psuedograph(13) = 3;
        psuedograph(14) = 4;
        psuedograph(15) = 4;
        psuedograph(16) = 4;
        psuedograph(17) = 4;
        psuedograph(18) = 4;
        psuedograph(19) = 4;
        psuedograph(20) = 5;
        psuedograph(21) = 5;
        psuedograph(22) = 5;
        psuedograph(23) = 5;
        psuedograph(24) = 5;
        psuedograph(25) = 5;
        psuedograph(26) = 5;
        psuedograph(27) = 6;
        psuedograph(28) = 6;
        psuedograph(29) = 6;
        psuedograph(30) = 6;
        psuedograph(31) = 6;
        psuedograph(32) = 6;
    case 7
        psuedograph = zeros(69,1);
        psuedograph(1) = 1;
        psuedograph(2) = 2;
        psuedograph(3) = 3;        
        psuedograph(4) = 2;
        psuedograph(5) = 2;
        psuedograph(6) = 3;
        psuedograph(7) = 3;
        psuedograph(8) = 3;
        psuedograph(9) = 3;
        psuedograph(10) = 3;
        psuedograph(11) = 3;
        psuedograph(12) = 3;
        psuedograph(13) = 3; 
        psuedograph(14) = 3;
        psuedograph(15) = 4;
        psuedograph(16) = 4;
        psuedograph(17) = 4;
        psuedograph(18) = 4;
        psuedograph(19) = 4;
        psuedograph(20) = 4;
        psuedograph(21) = 4;
        psuedograph(22) = 4;
        psuedograph(23) = 4;
        psuedograph(24) = 4;
        psuedograph(25) = 4;
        psuedograph(26) = 4;
        psuedograph(27) = 4;
        psuedograph(28) = 4;
        psuedograph(29) = 5;
        psuedograph(30) = 5; 
        psuedograph(31) = 5;
        psuedograph(32) = 5;
        psuedograph(33) = 5;
        psuedograph(34) = 5;
        psuedograph(35) = 5;
        psuedograph(36) = 5;
        psuedograph(37) = 5;
        psuedograph(38) = 5;
        psuedograph(39) = 5;
        psuedograph(40) = 5;
        psuedograph(41) = 5;
        psuedograph(42) = 5;
        psuedograph(43) = 5;
        psuedograph(44) = 5;
        psuedograph(45) = 6;
        psuedograph(46) = 6;
        psuedograph(47) = 6;
        psuedograph(48) = 6;
        psuedograph(49) = 6;
        psuedograph(50) = 6;
        psuedograph(51) = 6;
        psuedograph(52) = 6;
        psuedograph(53) = 6;
        psuedograph(54) = 6;
        psuedograph(55) = 6;
        psuedograph(56) = 6;
        psuedograph(57) = 6;
        psuedograph(58) = 6;
        psuedograph(59) = 6;
        psuedograph(60) = 6;
        psuedograph(61) = 7;
        psuedograph(62) = 7;
        psuedograph(63) = 7;
        psuedograph(64) = 7;
        psuedograph(65) = 7;
        psuedograph(66) = 7;
        psuedograph(67) = 7;
        psuedograph(68) = 7;
        psuedograph(69) = 7;
end


