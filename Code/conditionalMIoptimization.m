% This is a function to calculate conditional mutual information from low
% resolution images to optimize acquisition parameters

patientnum=1;

%% Default Acquisition Parameters

flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
TD = [0.5,0.5,0.5,0.5];          % s

acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];

%% Population Tissue Statistics
% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
T1mean = [1200,  900, 4000, 1200]./1000; % s
T1stdd = [ 100,  100,  200,  150]./1000; % s

T2mean = [ 100,   80, 1000,  110]./1000; % s
T2stdd = [   5,    4,   50,   10]./1000; % s

M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity

tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];

%% Load Population Model
lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
tmptissue=load_untouch_nii(lfname);
% materialID=tmptissue.img;
materialID=int32(tmptissue.img(15:165,20:200,92));

%% Load Patient
load('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/syntheticPatientPopulation.mat');
patientmatID_hires=syntheticPatientPop{patientnum,1}(15:165,20:200,92);
patientmatID_lores=syntheticPatientPop{patientnum,2}(4:42,5:50,23);
patient_tisinput=syntheticPatientPop{patientnum,3};

%% Optimization Parameters
pspacelabels={'TD(1)','TD(2)','TD(3)','TD(4)'};
subsmpllabels={};   %{'variance'};
% pinit=tconstrain/4*ones([4,1]);
% pAeq=zeros(length(pinit));
% pAeq(2,:)=[1,1,1,1];
% pbeq=zeros(size(pinit));
% pbeq(2)=tconstrain;
% pmin=[0,0,0,0]';
% pmax=tconstrain*[1,1,1,1]';
% findiffrelstep=1.e-6;
% tolx=1.e-2;%1.e-5;
% tolfun=1.e-2;%1.e-5;
% maxiter=500;

if exist('opt_history.txt','file')==2
    delete('opt_history.txt');
end

%% Calculate MI on Population Model
% tic;
% [popt,fval,exitflag,output,lambda,grad,hessian]=fmincon(@(x) conditionalMI_objfun(x,pspacelabels,subsmpllabels,tisinput,acqparam,materialID),...
%     pinit,pAeq,pbeq,[],[],pmin,pmax,[],...
%     optimset('FinDiffRelStep',findiffrelstep,'TolX',tolx,'TolFun',tolfun,'MaxIter',maxiter,'Display','iter-detailed','OutputFcn',@outfun));
% toc;
% [popt,fval]=fmincon(@(x) conditionalMI_objfun(x,pspacelabels,subsmpllabels,tisinput,acqparam,materialID),...
%     pinit,[],[],[],[],[],[],[],...
%     optimset('Display','iter-detailed'));
%% Calculate Conditional MI on Low Res Patient Image
delayit=0:.02:1;
for iii=1:10
    patientmatID_lores=syntheticPatientPop{iii,2}(4:42,5:50,23);
    patient_tisinput=syntheticPatientPop{iii,3};
    for jjj=1:51
    iii
    jjj
        pinit=delayit(jjj)*ones([1,4]);
        mi(iii,jjj)=conditionalMI_objfun(pinit,pspacelabels,subsmpllabels,patient_tisinput,acqparam,patientmatID_lores);
    end
end