%% this produces tables in subsection 4.4 [rough Heston]
%  you can play around with the maturity to move from one table to another
%  additional comments can be found in main_convGBM (all scripts
%  main_convMETHOD follow the same logic)

clear
clc
close all

model = 'rough Heston';

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
    
H   =  0.05;
rho = -0.65;
nu  =  0.40;
params= [H nu rho];

c = 0.0256;
n = 100;
xi = repmat(c,1,n);

if t == 0.01
    Xc = 2.7074;
elseif t == 0.10
    Xc = 7.6532;
elseif t == 1.00
    Xc = 18.9469;
elseif t == 10
    Xc = 61.4549;
end
Xm = 0;

NS = 2^17;

N_vec = [256 512 768 1024 1536 2048];
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
        convtab((q-1)*2+1,i) = SINC_discFT(S,t,K(q),IR,DY,model,params,xi,n,Xc,Xm,NS,CP,type);
        convtab((q-1)*2+2,i) = COSm_discFT(S,t,K(q),IR,DY,model,params,xi,n,Xc,Xm,NC,CP,type);
    end
end

NC = 2^16;
if type == 1
    NS = 2*NC;
elseif type == 2 || type == 3
    NS = 4*NC;
end
benchmarkS = SINC_discFT(S,t,K,IR,DY,model,params,xi,n,Xc,Xm,NS,CP,type);
benchmarkC = COSm_discFT(S,t,K,IR,DY,model,params,xi,n,Xc,Xm,NC,CP,type);

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
