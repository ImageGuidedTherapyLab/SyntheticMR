
%% Make In Vivo Predictions
T1mean = [1200,  900, 900, 1200]./1000; % s
T1stdd = [ 100,  100,  100,  150]./1000; % s
T2mean = [ 100,   80, 80,  110]./1000; % s
T2stdd = [   5,    4,   4,   10]./1000; % s
M0mean = [ 0.9,  0.9,  0.9,  0.9];       % relative intensity
M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];
tisinputnovar=[M0mean;0.001*ones(size(M0stdd));T1mean;0.001*ones(size(T1stdd));T2mean;0.001*ones(size(T2stdd))];

flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.871;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
TD = [0.5,0.5,0.5,0.5,0.5];          % s
signu=3.4762E-4;
% signu=0.40;
acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD,signu];

lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
tmptissue = load_untouch_nii(lfname);
materialID=int32(tmptissue.img(15:165,20:200,92));

% Theoretical Optimal Case
mi0=MI_GHQuad_3DQALAS(TD,tisinput,acqparam,materialID);
%popt=MI_GHQuad_3DQALAS_driver(acqparam,tisinput,materialID);
%mipred_topt=MI_GHQuad_3DQALAS(popt,tisinput,acqparam,materialID);
