%% draw PDFs for CGMY and rough Heston model [Figure 1 and 7 in the paper]

clear
clc
close all

%% -------------------------------- CGMY -------------------------------- %
model = 'CGMY';

% t = 0.01;
t = 1.00;

C = 1;
G = 5;
M = 5;
Y = 0.5;
params = [C G M Y];

if t == 0.01
    Xc = 12.0723;
elseif t == 1
    Xc = 18.3512;
end
Xm = 0;

N = 2^20;

[x,pdf] = chf2pdf(model,params,t,[],[],Xc,Xm,N);

figure
plot(x,real(pdf),'k','Linewidth',1)
xlim([-3 3])
xlabel('${\it} s_T$','Interpreter','Latex','FontSize', 15)
ylabel('PDF')


%% ---------------------------- rough Heston ---------------------------- %
model = 'rough Heston';

% t = 0.01;
t = 1.00;

H   =  0.05;
rho = -0.65;
nu  =  0.40;
params= [H nu rho];

c = 0.0256;
j = 100;
xi = repmat(c,1,j);

if t == 0.01
    Xc = 2.7074;
elseif t == 1.00
    Xc = 18.9469;
end
Xm = 0;

N = 2^18;

[x,pdf] = chf2pdf(model,params,t,xi,j,Xc,Xm,N);

figure
plot(x,real(pdf),'k','Linewidth',1)
xlim([-2 2])
if t == 0.01
    ylim([0 60])
end
xlabel('${\it} s_T$','Interpreter','Latex','FontSize', 15)
ylabel('PDF')
