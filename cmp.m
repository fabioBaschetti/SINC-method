% This script compares two diffent strategies to FT pricing within the SINC 
% approach. They correspond to functions 'SINC_discFT' and 'SINC_discFT_sbs'
% respectively. 
% The enfasis is in the CPU time for both functions to complete their tasks
% and the relative accordance on the output prices. You can change the
% number of Fourier modes and the multiplier L (that determines the length
% of the truncation range) depending on yor needs, and evaluate the
% relative performance of the two strategies.
% Experiments are under the rough Heston model, again, and the forward variance 
% curve assumed to be flat (for simplicity).

clc
close all
clear 

%% FT pricing for a given strike (CoN, AoN and put options)

% rough Heston parameters 
H   =  0.05;
rho = -0.65;
nu  =  0.40;
params_rHeston = [H nu rho];

% (flat) forward variance curve \xi_0(t)  
c = 0.0256;
j = 100;
xi = repmat(c,1,j);

% option's specifics
fwd = 1;
t = 1;
K = 1;

% truncation range for the PDF of s_T
bl = -80;
bh = -bl;
Nt = 2^16;
L = 100;
[Xc,Xm] = truncMeasures(params_rHeston,j,xi,t,bl,bh,Nt,L);

% Fourier modes for put options
N = 2^16;

% digital and plain vanilla put price
tic
[con1,aon1,put1] = SINC_discFT(fwd,t,K,params_rHeston,j,xi,Xc,Xm,N);
toc
tic
[con2,aon2,put2] = SINC_discFT_sbs(fwd,t,K,params_rHeston,j,xi,Xc,Xm,N);
toc

abs(con1-con2)
abs(aon1-aon2)
abs(put1-put2)