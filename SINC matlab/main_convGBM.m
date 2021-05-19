%% this produces tables in subsection 4.1 [GBM]
%  you can play around with the maturity to move from one table to another

clear
clc
close all

model = 'GBM';

CP = -1;        % +1: call  -1:put
type = 1;       % 1: plain vanilla   2: cash or nothing   3: asset or nothing

% option's specifics
S = 1;
% t = 0.01;
t = 0.10;
% t = 1.00;
% t = 10;
K = S*(0.60:0.10:1.40);

IR = 0.1;
DY = 0;

% model parameters
sigma = 0.25;
params = sigma;

% calculated Xc
if t == 0.01
    Xc = 1.3863;
elseif t == 0.10
    Xc = 2.0105;
elseif t == 1.00
    Xc = 6.3577;
elseif t == 10
    Xc = 20.1036;
end
Xm = 0;

% grid for the number of evaluations of the CF
N_vec = [20 40 60 80 100 120];
m = length(N_vec);

Q = length(K);

% populate a table of SINC and COS prices where N = N_vec
convtab = zeros(2*Q,m);
for q = 1:Q
    for i = 1:m
        NC = N_vec(i);
        if type == 1
            NS = 2*NC;
        elseif type == 2 || type == 3
            NS = 4*NC;
        end
        convtab((q-1)*2+1,i) = SINC_discFT(S,t,K(q),IR,DY,model,params,[],[],Xc,Xm,NS,CP,type);
        convtab((q-1)*2+2,i) = COSm_discFT(S,t,K(q),IR,DY,model,params,[],[],Xc,Xm,NC,CP,type);
    end
end

% produce a benchmark
benchmark = BS_closedForm(S,K,t,IR,DY,sigma,CP,type);

benchtab  = [];
for q = 1:Q
    benchtab = [benchtab; benchmark(q)*ones(2,m)];
end

% compute absolute errors
abstab = abs(benchtab-convtab);

% compute relative errors
reltab = abstab ./ benchtab;

% the final output
for a = 1:2*length(K)
    for b = 1:length(N_vec)
        if reltab(a,b) >= 1
            reltab(a,b) = 1;
        end
    end 
end
