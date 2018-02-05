
global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

clearvars;

% Fixed parameters
flipangle = 120/180.*pi; % radian
excitationangle =  90/180.*pi; % radian
repetitiontime  = 4000; % ms

% Reference points (Calculate MI at each reference point in parameter
% space)
delaytime       = [0:200:4000 ] ; % ms
echotime        = [0:7:140 ] ; % ms

tic;
% experimental data to be fitted
necho = [1; 2; 3; 4; 5];
Mmeas = [889; 2134; 3021; 3684; 3654];     
sigmeas = [0.1; 0.1; 0.1; 0.1; 0.1];

data{1} = necho;
data{2} = Mmeas;
data{3} = diag(sigmeas.^2); %covariance matrix

Nlive = 10;
Nmcmc = 0; %input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_gaussian;
model = @mnqalasmodel;
% parameters and range
prior = {'M0', 'uniform', 3500, 3550, 'fixed'; ...
         'T1', 'uniform', 0, .1, 'fixed'; ...
         'T2', 'uniform', 0, 1, 'fixed'};
extraparams = {}; % no extra parameters beyond what's in the prior

% use nested_sampler
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

toc;

%% Measurement model
% clear all
% close all
% format shortg

function [Kspace,KspaceCenter,MI,H_z]=MagicSignalModel1p(materialID,NumGP)

%% signal model parameters
%%         GM      WM      CSF     Tumor
T1mean = [1400  , 1000   ,  4000   , 250  ]; % ms
T1stdd = [100   ,  100   ,  200    ,  50  ]; % ms
% T1stdd = [0,0,0,0];
T2mean = [100   , 75   , 600   , 250  ]; % ms
T2stdd = [5    ,  5   ,  30   ,  50  ]; % ms
% T2stdd=[0,0,0,0];
M0mean = [0.9   , 0.8   , 1.0   , 1.2  ]; % relative intensity
M0stdd = [.05    ,  .05   ,  .05   ,  .4  ]; % relative intensity
noisestdd = 0.01;

%% Quadrature Points
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumGP,[mean(T1mean(1:2)),T1mean(3),mean(T2mean(1:2)),T2mean(3)],[mean(T1stdd(1:2)),T1stdd(3),mean(T2stdd(1:2)),T2stdd(3)]);
[xz,xnz,xmz,wz,wnz]=GaussHermiteNDGauss(NumGPz,[0,0],[noisestdd,noisestdd]);

%% loop over delay time and echo time

%% Generate parameters maps from mean
T1map = eps * (materialID==0) + xm(gp1,1) * or(materialID == 1 , materialID == 2) + xm(gp2,2) * (materialID == 3) ;
T2map = eps * (materialID==0) + xm(gp3,3) * or(materialID == 1 , materialID == 2) + xm(gp4,4) * (materialID == 3) ;
M0map = M0mean(1) * (materialID == 1) + M0mean(2) * (materialID == 2) + M0mean(3) * (materialID == 3) ;

disp( sprintf( '%d of %d, %d of %d',iii ,length(delaytime ) ,jjj,length(echotime        )))
eTD = exp(-delaytime(iii)./T1map );
eTR = exp(-repetitiontime./T1map );
eTE = exp(-echotime(jjj)./T2map );
% compute magnetization image in image domain
Mmodel = M0map.*( 1 - (1 - cos(flipangle)) * eTD  - cos(flipangle) * eTR ).*((1- cos(flipangle)*eTR*cos(excitationangle) ).^(-1) ).*eTE ;
% tranform to fourier space. note this is a function of TD and TE
Kspace{iii,jjj,gp1,gp2,gp3,gp4} = fftshift(fftn(Mmodel));
ksc=floor(size(Kspace)/2+1);
KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4)=Kspace{iii,jjj,gp1,gp2,gp3,gp4}(1);
