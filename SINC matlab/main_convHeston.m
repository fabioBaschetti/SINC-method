%% this produces tables in subsection 4.2 [Heston]
%  you can play around with the maturity to move from one table to another
%  additional comments can be found in main_convGBM (all scripts
%  main_convMETHOD follow the same logic)

clear
clc
close all

model = 'Heston';

CP = -1;
type = 1;

S = 1;
% t = 0.01;             
% t = 0.10;
t = 1.00;
% t = 10;
K = S*(0.60:0.10:1.40);

IR = 0;
DY = 0;

lambda = 1.5768;
v_bar  = 0.0398;
eta    = 0.5751;
v0     = 0.0175;
rho    = -.5711;
params = [lambda v_bar eta v0 rho];

if t == 0.01
    Xc = 1.3863;
elseif t == 0.10
    Xc = 2.0499;
elseif t == 1.00
    Xc = 12.1802;
elseif t == 10
    Xc = 32.8824;
end
Xm = 0;

if t == 0.1;
    N_vec = [64 128 192 256 384 512];
elseif t == 1
    N_vec = [128 192 256 384 512 768];
end
m = length(N_vec);

Q = length(K);

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

NC = 2^18;
if type == 1
    NS = 2*NC;
elseif type == 2 || type == 3
    NS = 4*NC;
end
benchmarkS = SINC_discFT(S,t,K,IR,DY,model,params,[],[],Xc,Xm,NS,CP,type);
benchmarkC = COSm_discFT(S,t,K,IR,DY,model,params,[],[],Xc,Xm,NC,CP,type);

diff_SINCvsCOS = abs(benchmarkS - benchmarkC);

d = max(diff_SINCvsCOS);
o = max(order(d)+1,-10);

benchmark = fix(benchmarkS*10^(-o))*10^o;

benchtab  = [];
for q = 1:Q
    benchtab = [benchtab; benchmark(q)*ones(2,m)];
end

abstab = abs(benchtab-convtab);

reltab = abstab ./ benchtab;

convtab = fix(convtab*10^(-o))*10^o;

for a = 1:2*length(K)
    for b = 1:length(N_vec)
        if convtab(a,b) == benchtab(a,b)
            reltab(a,b) = 0;
        elseif reltab(a,b) >= 1
            reltab(a,b) = 1;
        else
            % do nothing
        end
    end 
end
