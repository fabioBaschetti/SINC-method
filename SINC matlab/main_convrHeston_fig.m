%% study step-by-step convergence of FT-SINC and COS (as explained in subsection 4.4)
%  [experiments in the rough Heston model]
%  when t = 1.00 this produces figure 2 in the paper
%  when t = 0.01 this produces figure 5 in the paper

clc
close all
clear

model = 'rough Heston';

S = 1.00;
t = 0.01;
% t = 1.00;
K = 0.80;

IR = 0;
DY = 0;

H   =  0.05;
rho = -0.65;
nu  =  0.40;
params= [H nu rho];

j  = 100;
xi = 0.0256*ones(1,j);

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

N = 2^16;

[conBS,aonBS,putBS,~] = SINC_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,4*N);
[conBC,aonBC,putBC,~] = COSm_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,  N);

con_diff = abs(conBS - conBC);
aon_diff = abs(aonBS - aonBC);
put_diff = abs(putBS - putBC); 

d = max([con_diff,aon_diff,put_diff]);
o = max(order(d)+1,-10);

conB  = fix(conBS*10^(-o))*10^o;
aonB  = fix(aonBS*10^(-o))*10^o;
putB  = fix(putBS*10^(-o))*10^o;

clear conBS conBC aonBS aonBC putBS putBC con_diff aon_diff put_diff d  

tol = 1e-08; 

[conS(1),aonS(1),putS(1),~] = SINC_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,4);
nS(1) = 1;

on = 1;
N  = 2;
while on == 1 
    [conS(N),aonS(N),putS(N),~] = SINC_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,4*N);
    nS(N) = N;
    if (abs(conS(N)-conB)<tol) && (abs(aonS(N)-aonB)<tol)
        on = 0;
    end
    N = N +1;
end

[conC(1),aonC(1),putC(1),~] = COSm_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,1);
nC(1) = 1;

on = 1;
N  = 2;
while on == 1 
    [conC(N),aonC(N),putC(N),~] = COSm_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,N);
    nC(N) = N;
    if (abs(conC(N)-conB)<tol) && (abs(aonC(N)-aonB)<tol)
        on = 0;
    end
    N = N +1;
end

% plot convergence results for SINC
figure
grid on
plot(log2(nC),zeros(1,size(nC,2)),'w','LineWidth',1)
hold on
plot(log2(nS),repmat(conB,1,size(nS,2)),'c','LineWidth',0.5)
hold on
plot(log2(nS),repmat(aonB,1,size(nS,2)),'c','LineWidth',0.5)
hold on
plot(log2(nS),repmat(putB,1,size(nS,2)),'c','LineWidth',0.5)
hold on
CoN = plot(log2(nS),conS,'r--','LineWidth',1);
hold on
AoN = plot(log2(nS),aonS,'b-.','LineWidth',1);
hold on
PUT = plot(log2(nS(2:2:end)),putS(1:floor(nS(end)/2)),'k','LineWidth',1);
xlabel('${\it} \log_2 N_{F}$','Interpreter','Latex','FontSize', 15)
ylabel('SINC price')
ylim([-0.1 0.5])
legend([CoN,AoN,PUT],'CoN','AoN','PUT')

% plot convergence results for COS
figure
grid on
plot(log2(nC),zeros(1,size(nC,2)),'w','LineWidth',1)
hold on
plot(log2(nC),repmat(conB,1,size(nC,2)),'c','LineWidth',0.5)
hold on
plot(log2(nC),repmat(aonB,1,size(nC,2)),'c','LineWidth',0.5)
hold on
plot(log2(nC),repmat(putB,1,size(nC,2)),'c','LineWidth',0.5)
hold on
CoN = plot(log2(nC),conC,'r--','LineWidth',1);
hold on
AoN = plot(log2(nC),aonC,'b-.','LineWidth',1);
hold on
PUT = plot(log2(nC),putC,'k','LineWidth',1);
xlabel('${\it} \log_2 N_{F}$','Interpreter','Latex','FontSize', 15)
ylabel('COS price')
ylim([-0.1 0.5])
legend([CoN,AoN,PUT],'CoN','AoN','PUT')
