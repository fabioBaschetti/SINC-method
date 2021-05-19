%% this relates to the second set of experiments in section 5 [SPX 20170719]  
%  table 8 is based on this
%  you can read vectors errVol_METHOD to check average absolute errors over
%  each smile
%  suffix d,f,r associate to FT,FFT,frFFT methods, respectively

clear
clc
close all

model = 'rough Heston';

% parameter set
H   = 0.05;
rho = -.65;
nu  = 0.40;
params = [H nu rho];

% forward variance curve (assume: flat)
j = 100;
xi = 0.0256*ones(1,j);

IR = 0;
DY = 0;

surface17 = readtable('spx_surf20170519.csv');

Xc = [ 2.455837865542;  3.156996686346;  3.709870562801;  4.578616822386;  4.763278907139;  5.104646394964;
       5.560795007528;  5.835434846020;  6.090941161903;  6.444406529772;  6.663247943194;  6.870617852062;
       7.162990754165;  7.347021253525;  7.523383143325;  7.775128620353;  7.935313060025;  8.090036165588;
       8.594476977920;  9.051906585519;  9.472397938240; 10.022999158715; 10.904012452716; 11.475423791431;
      12.084309867048; 12.613916316535; 13.717700918689; 15.090868558866; 15.486298567033; 16.056942721061;
      17.470790598461; 17.780505063733; 19.525948065635; 23.209012598363; 29.146196447375];
Xm = 0;

expiries = unique(surface17.t,'stable');
L = size(expiries,1);

counter = 1:L;

M = size(surface17.fwd,1);

CP = -1;

% % -- compute benchmark put and impVol to pupulate cols 3 and 4 in surface13.csv -- %
% % -- we check the agreement of FT-SINC vs COS and use the latter as a benchmark -- % 
% Nb = 2^18;
% 
% n   = 0;
% put_SINC = zeros(M,1);
% put_COSm = zeros(M,1);
% surfVol_SINC = zeros(M,1);
% surfVol_COSm = zeros(M,1);
% for i = 1:length(counter)
%     slice  = counter(i);
%     expiry = expiries(slice);
%     pick = surface17.t == expiry;
%     temp = surface17(pick,:);
%     fwd  = unique(temp.fwd);
%     t    = unique(temp.t);
%     K    = temp.K;
%     m    = length(K);
%     smile_SINC = SINC_discFT(fwd,t,K,IR,DY,model,params,xi,j,Xc(slice),Xm,2*Nb,CP,1);
%     smile_COSm = COSm_discFT(fwd,t,K,IR,DY,model,params,xi,j,Xc(slice),Xm,  Nb,CP,1);
%     put_SINC(n+1:n+m) = smile_SINC;
%     put_COSm(n+1:n+m) = smile_COSm;
%     smileVol_SINC = BSImpliedVolPut(fwd,K,t,0,smile_SINC);
%     surfVol_SINC(n+1:n+m) = smileVol_SINC;
%     smileVol_COSm = BSImpliedVolPut(fwd,K,t,0,smile_COSm);
%     surfVol_COSm(n+1:n+m) = smileVol_COSm;
%     n    = n + m;
% end
% % ---------------------------------------------------------------------- %


%% ---------------------------- FT  methods ---------------------------- %
powSd = 11;
NSd = 2^powSd;

powCd = 11;
NCd = 2^powCd;

powLd = 13;
NLd = 2^powLd;
aa_Ld = 1.5;
ub_Ld = (NLd-1)./(2*Xc*aa_Ld);

powMd = 14;
NMd = 2^powMd;
aa_Md = 2.975;
ub_Md = (NMd-1)./(2*Xc*aa_Md);
dump = 0.4;

n   = 0;
put_SINCd = zeros(M,1);
put_COSmd = zeros(M,1);
put_Lewid = zeros(M,1);
put_Carrd = zeros(M,1);
surfVol_SINCd = zeros(M,1);
surfVol_COSmd = zeros(M,1);
surfVol_Lewid = zeros(M,1);
surfVol_Carrd = zeros(M,1);
errVol_SINCd = zeros(L,1);
errVol_COSmd = zeros(L,1);
errVol_Lewid = zeros(L,1);
errVol_Carrd = zeros(L,1);
for i = 1:length(counter)
    slice  = counter(i);
    expiry = expiries(slice);
    pick = surface17.t == expiry;
    temp = surface17(pick,:);
    fwd  = unique(temp.fwd);
    t    = unique(temp.t);
    K    = temp.K;
    m    = length(K);
    impVol_bench = temp.impVol;
    smile_SINCd = SINC_discFT(fwd,t,K,IR,DY,model,params,xi,j,Xc(slice),Xm,2*NSd,CP,1);
    put_SINCd(n+1:n+m) = smile_SINCd;
    smileVol_SINCd = BSImpliedVolPut(fwd,K,t,IR,smile_SINCd);
    surfVol_SINCd(n+1:n+m) = smileVol_SINCd;
    errVol_SINCd(i) = mean(abs(impVol_bench - smileVol_SINCd));
    smile_COSmd = COSm_discFT(fwd,t,K,IR,DY,model,params,xi,j,Xc(slice),Xm,  NCd,CP,1);
    put_COSmd(n+1:n+m) = smile_COSmd;
    smileVol_COSmd = BSImpliedVolPut(fwd,K,t,IR,smile_COSmd);
    surfVol_COSmd(n+1:n+m) = smileVol_COSmd;
    errVol_COSmd(i) = mean(abs(impVol_bench - smileVol_COSmd));
    smile_Lewid = Lewis_discFT(fwd,t,K,IR,DY,model,params,xi,j,ub_Ld(slice),NLd,CP);
    put_Lewid(n+1:n+m) = smile_Lewid;
    smileVol_Lewid = BSImpliedVolPut(fwd,K,t,IR,smile_Lewid);
    surfVol_Lewid(n+1:n+m) = smileVol_Lewid;
    errVol_Lewid(i) = mean(abs(impVol_bench - smileVol_Lewid));
    smile_Carrd = Carr_discFT(fwd,t,K,IR,DY,model,params,xi,j,dump,ub_Md(slice),NMd,CP);
    put_Carrd(n+1:n+m) = smile_Carrd;
    smileVol_Carrd = BSImpliedVolPut(fwd,K,t,IR,smile_Carrd);
    surfVol_Carrd(n+1:n+m) = smileVol_Carrd;
    errVol_Carrd(i) = mean(abs(impVol_bench - smileVol_Carrd));
    n    = n + m;
end
% ---------------------------------------------------------------------- %


%% ---------------------------- FFT methods ---------------------------- %
powSf = 14;

powLf = 16;
NLf = 2^powLf;
aa_Lf = 1;
ub_Lf = (NLf-1)./(2*Xc*aa_Lf);

powMf = 16;
NMf = 2^powMf;
aa_Mf = 2.5;
ub_Mf = (NMf-1)./(2*Xc*aa_Mf);
dump = 0.4;

n   = 0;
put_SINCf = zeros(M,1);
put_Lewif = zeros(M,1);
put_Carrf = zeros(M,1);
surfVol_SINCf = zeros(M,1);
surfVol_Lewif = zeros(M,1);
surfVol_Carrf = zeros(M,1);
errVol_SINCf = zeros(L,1);
errVol_Lewif = zeros(L,1);
errVol_Carrf = zeros(L,1);
for i = 1:length(counter)
    slice  = counter(i);
    expiry = expiries(slice);
    pick = surface17.t == expiry;
    temp = surface17(pick,:);
    fwd  = unique(temp.fwd);
    t    = unique(temp.t);
    K    = temp.K;
    m    = length(K);
    impVol_bench = temp.impVol;
    smile_SINCf = SINC_fastFT(fwd,t,K,IR,DY,model,params,xi,j,2*Xc(slice),Xm,powSf+1,CP);
    put_SINCf(n+1:n+m) = smile_SINCf;
    smileVol_SINCf = BSImpliedVolPut(fwd,K,t,IR,smile_SINCf);
    surfVol_SINCf(n+1:n+m) = smileVol_SINCf;
    errVol_SINCf(i) = mean(abs(impVol_bench - smileVol_SINCf));
    smile_Lewif = Lewis_fastFT(fwd,t,K,IR,DY,model,params,xi,j,ub_Lf(slice),powLf,CP);
    put_Lewif(n+1:n+m) = smile_Lewif;
    smileVol_Lewif = BSImpliedVolPut(fwd,K,t,IR,smile_Lewif);
    surfVol_Lewif(n+1:n+m) = smileVol_Lewif;
    errVol_Lewif(i) = mean(abs(impVol_bench - smileVol_Lewif));
    smile_Carrf = Carr_fastFT(fwd,t,K,IR,DY,model,params,xi,j,dump,ub_Mf(slice),powMf,CP);
    put_Carrf(n+1:n+m) = smile_Carrf;
    smileVol_Carrf = BSImpliedVolPut(fwd,K,t,IR,smile_Carrf);
    surfVol_Carrf(n+1:n+m) = smileVol_Carrf;
    errVol_Carrf(i) = mean(abs(impVol_bench - smileVol_Carrf));
    n    = n + m;
end
% ---------------------------------------------------------------------- %


%% --------------------------- frFT  methods --------------------------- %
powSr = 11;
e_Sr = 0.15;

powLr = 13;
NLr = 2^powLr;
aa_Lr = 1.17;
ub_Lr = (NLr-1)./(2*Xc*aa_Lr);
e_Lr = 0.02;

powMr = 14;
NMr = 2^powMr;
aa_Mr = 3.168;
ub_Mr = (NMr-1)./(2*Xc*aa_Mr);
e_Mr = 0.02;
dump = 0.4;

n   = 0;
put_SINCr = zeros(M,1);
put_Lewir = zeros(M,1);
put_Carrr = zeros(M,1);
surfVol_SINCr = zeros(M,1);
surfVol_Lewir = zeros(M,1);
surfVol_Carrr = zeros(M,1);
errVol_SINCr = zeros(L,1);
errVol_Lewir = zeros(L,1);
errVol_Carrr = zeros(L,1);
for i = 1:length(counter)
    slice  = counter(i);
    expiry = expiries(slice);
    pick = surface17.t == expiry;
    temp = surface17(pick,:);
    fwd  = unique(temp.fwd);
    t    = unique(temp.t);
    K    = temp.K;
    m    = length(K);
    impVol_bench = temp.impVol;
    smile_SINCr = SINC_fracFT(fwd,t,K,IR,DY,model,params,xi,j,Xc(slice),Xm,powSr+1,e_Sr,CP);
    put_SINCr(n+1:n+m) = smile_SINCr;
    smileVol_SINCr = BSImpliedVolPut(fwd,K,t,IR,smile_SINCr);
    surfVol_SINCr(n+1:n+m) = smileVol_SINCr;
    errVol_SINCr(i) = mean(abs(impVol_bench - smileVol_SINCr));
    smile_Lewir = Lewis_fracFT(fwd,t,K,IR,DY,model,params,xi,j,ub_Lr(slice),powLr,e_Lr,CP);
    put_Lewir(n+1:n+m) = smile_Lewir;
    smileVol_Lewir = BSImpliedVolPut(fwd,K,t,IR,smile_Lewir);
    surfVol_Lewir(n+1:n+m) = smileVol_Lewir;
    errVol_Lewir(i) = mean(abs(impVol_bench - smileVol_Lewir));
    smile_Carrr = Carr_fracFT(fwd,t,K,IR,DY,model,params,xi,j,dump,ub_Mr(slice),powMr,e_Mr,CP);
    put_Carrr(n+1:n+m) = smile_Carrr;
    smileVol_Carrr = BSImpliedVolPut(fwd,K,t,IR,smile_Carrr);
    surfVol_Carrr(n+1:n+m) = smileVol_Carrr;
    errVol_Carrr(i) = mean(abs(impVol_bench - smileVol_Carrr));
    n    = n + m;
end
% ---------------------------------------------------------------------- %
