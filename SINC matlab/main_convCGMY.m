%% this produces tables in subsection 4.3 [CGMY]
%  you can play around with the maturity and parameter Y to move from one
%  table to another
%  additional comments can be found in main_convGBM (all scripts
%  main_convMETHOD follow the same logic)

clear
clc
close all

model = 'CGMY';

CP = -1;    
type = 1;

S = 1;
% t = 0.01;
% t = 0.10;
t = 1.00;
% t = 10;
K = S*(0.60:0.10:1.40);

IR = 0.1;
DY = 0;

C = 1;
G = 5;
M = 5;
% Y = 0.5;
Y = 1.5;
% Y = 1.98;
params = [C G M Y];

if t == 0.01
    if Y == 0.5
        Xc = 12.0723;
    elseif Y == 0.9
        Xc = 11.7596;
    elseif Y == 1.1
        Xc = 11.6256;
    elseif Y == 1.5
        Xc = 11.4582;
    elseif Y == 1.98
        Xc = 24.9357;
    end
elseif t == 0.10
    if Y == 0.5
        Xc = 14.0518;
    elseif Y == 0.9
        Xc = 13.9021;
    elseif Y == 1.1
        Xc = 13.9529;
    elseif Y == 1.5
        Xc = 14.9722;
    elseif Y == 1.98
        Xc = 78.7193;
    end
elseif t == 1
    if Y == 0.5
        Xc = 18.3512;
    elseif Y == 0.9
        Xc = 20.0115;
    elseif Y == 1.1
        Xc = 21.9519;
    elseif Y == 1.5
        Xc = 33.0891;
    elseif Y == 1.98
        Xc = 248.9047;
    end
elseif t == 10
    if Y == 0.5
        Xc = 36.1317;
    elseif Y == 0.9
        Xc = 47.7702;
    elseif Y == 1.1
        Xc = 58.2017;
    elseif Y == 1.5
        Xc = 101.5438;
    elseif Y == 1.98
        Xc = 0;
    end
end
Xm = 0;

if t == 0.01
    if Y == 0.5
        N_vec = [256 512 1024 2048 4096 8192];
    elseif Y == 1.5
        N_vec = [16 32 64 128 256 512];
    elseif Y == 1.98;
        N_vec = [16 32 48 64 96 128];
    end
elseif t == 1
    if Y == 0.5
        N_vec = [16 32 64 128 256 512];
    elseif Y == 1.5
        N_vec = [16 32 48 64 96 128];
    elseif Y == 1.98;
        N_vec = [16 32 48 64 96 128];
    end
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
if t == 0.01
    NC = 2^20;
end
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
