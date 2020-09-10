% This script provides a minimal working example for two exercises with the 
% SINC approach:
% 1) pricing digital and European put options for a given strike based on 
%    the FT-version of the method
% 2) pricing all options in a smile concurrently with the FFT-version of the
%    method.
% The underlying model is rough Heston, the forward variance curve assumed
% to be flat (for simplicity).

% RUN SECTIONS 1) AND 2) SEPARATELY

clc
close all
clear 

%% 1) FT pricing for a given strike (CoN, AoN and put options)

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
M = 2^16;
L = 100;
[Xc,Xm] = truncMeasures(params_rHeston,j,xi,t,bl,bh,M,L);

% Fourier modes for put options
N = 2^16;

% digital and plain vanilla put price
[con,aon,put] = SINC_discFT(fwd,t,K,params_rHeston,j,xi,Xc,Xm,N);

clc
close all
clear

%% 2) FFT pricing of one entire smile (CoN, AoN and put options)

smile = readtable('smile.csv');

% rough Heston parameters  
H   =  0.1216;
rho = -0.6714;
nu  =  0.2910;
params_rHeston = [H nu rho];

% (flat) forward variance curve \xi_0(t)  
c = 0.0256;
j = 100;
xi = repmat(c,1,j);

% option's specifics
t = unique(smile.t);
fwd = unique(smile.fwd);
K = smile.K;

clear smile

% truncation range for the PDF of s_T
bl = -80;
bh = -bl;
M = 2^16;
L = 100;
[Xc,Xm] = truncMeasures(params_rHeston,j,xi,t,bl,bh,M,L);

% Fourier modes for put options (N=2^pow)
pow = 16;

% digital and plain vanilla put price
[con,aon,put] = SINC_fastFT(fwd,t,K,params_rHeston,j,xi,Xc,Xm,pow);
