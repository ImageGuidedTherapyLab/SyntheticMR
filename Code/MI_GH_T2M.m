%% Measurement model
% clear all
% close all
% format shortg

function [Mmodel,MI,Hz]=MI_GH_T2M()

%% signal model parameters
%%          GM       WM        CSF     Tumor
T1mean = [1400   , 1000    ,  4000    , 250  ]./1000; % s
T1stdd = [100    ,  100    ,   200    ,  50  ]./1000; % s
T2mean = [95    ,   75    ,   600    , 250  ]./1000; % s
T2stdd = [  5    ,    5    ,    30    ,  50  ]./1000; % s
M0mean = [  0.9  ,    0.8  ,     1.0  ,   1.2];       % relative intensity
M0stdd = [   .05 ,     .05 ,      .05 ,    .4];       % relative intensity

TE = [0:1:1000]./1000; % s

%% Quadrature Points
NumGP=10;
[x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumGP,T2mean(1),T2stdd(1));

%% loop over delay time and echo time
lte=length(TE);

for gp = 1:NumGP
    for iii = 1:lte
        disp(sprintf('Model eval: %d of %d, %d of %d',gp,NumGP,iii,lte))
        Mmodel(iii,gp) = M0mean(1).*exp(-TE(iii)./xm_t2(gp));
    end
end

msize=size(Mmodel);
msize=[msize(1),1];
Ezr=zeros(msize);
Sigrr=zeros(msize);
for iii=1:NumGP
    disp(sprintf('Quadrature: %d of %d',iii,NumGP))
    Ezr=Ezr+wn_t2(iii).*Mmodel(:,iii);
    Sigrr=Sigrr+wn_t2(iii).*Mmodel(:,iii).^2;
end

N=1;
signu=0.1;
varz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2);
Hz=0.5.*log(2*pi*2.7183.*varz);
Hzmu=0.5.*log(2*pi*2.7183.*signu.^2);
MI=Hz-Hzmu;

% save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas_mi_results_nq%d.mat',NumGP),'kspace','MI','Hz','-v7.3');
