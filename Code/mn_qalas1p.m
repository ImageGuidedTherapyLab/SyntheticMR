% example: (from hogg et al., 1008.4686)
% fit a line to data omitting outlier points

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

pdv=-1;
tconoverride=1;
ttotal=2;

% read sample data for fitting a line
% data = readdata_line;

%% Optimization Space Acquisition Parameters
geometrycase=2;
lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
switch geometrycase
    case 1
        materialID = int32(1);
    case 1.5
        materialID = int32([1,2,3]);
    case 1.75
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID=int32(tmptissue.img(:,:,ceil(size(tmptissue.img,3)/2)));
    case 2
        tmptissue = load_untouch_nii(lfname);
        materialID=int32(tmptissue.img(15:165,20:200,92));
    case 2.5
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img);
    case 0
        %                 system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 1x1x3vox -interpolation NearestNeighbor -resample 256x192x64 -o resampleimg.nii.gz',lfname));
        %                 tmptissue = load_untouch_nii('resampleimg.nii.gz');
        %                 system('rm resampleimg.nii.gz');
        %                 materialID = int32(tmptissue.img);
        materialID=1;
    otherwise
        tmptissue = load_untouch_nii(lfname);
        materialID = int32(tmptissue.img);
end

%% Tissue Properties
% M0/T1/T2 Variance Flags
M0varflag = 1;
T1varflag = 1;
T2varflag = 1;

% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
T1mean = [1200,  900, 4000, 1200]./1000; % s
if T1varflag~=0
    T1stdd = [ 100,  100,  200,  150]./1000; % s
else
    T1stdd = [   0,    0,    0,    0];
end
T2mean = [ 100,   80, 1000,  110]./1000; % s
if T2varflag~=0
    T2stdd = [   5,    4,   50,   10]./1000; % s
else
    T2stdd = [   0,    0,    0,    0];
end
M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
if M0varflag~=0
    M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
else
    M0stdd = [   0,    0,    0,    0];
end

overwritecsf=1;
if overwritecsf==1
    T1mean = [1200,  900, 900, 1200]./1000; % s
    T1stdd = [ 100,  100,  100,  150]./1000; % s
    
    T2mean = [ 100,   80, 80,  110]./1000; % s
    T2stdd = [   5,    4,   4,   10]./1000; % s
    
    M0mean = [ 0.9,  0.9,  0.9,  0.9];       % relative intensity
    M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
end

tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];

%% Create synthetic data
stdmapT1=normrnd(0,1,size(materialID));
stdmapT2=normrnd(0,1,size(materialID));
stdmapM0=normrnd(0,1,size(materialID));

synthdataT1=(materialID==1).*(tisinput(3,1)+tisinput(4,1)*stdmapT1)+(materialID==2).*(tisinput(3,2)+tisinput(4,2)*stdmapT1)+(materialID==3).*(tisinput(3,3)+tisinput(4,3)*stdmapT1);
synthdataT2=(materialID==1).*(tisinput(5,1)+tisinput(6,1)*stdmapT2)+(materialID==2).*(tisinput(5,2)+tisinput(6,2)*stdmapT2)+(materialID==3).*(tisinput(5,3)+tisinput(6,3)*stdmapT2);
synthdataM0=(materialID==1).*(tisinput(1,1)+tisinput(2,1)*stdmapM0)+(materialID==2).*(tisinput(1,2)+tisinput(2,2)*stdmapM0)+(materialID==3).*(tisinput(1,3)+tisinput(2,3)*stdmapM0);
synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;

goldstandardT1=(materialID==1).*tisinput(3,1)+(materialID==2).*tisinput(3,2)+(materialID==3).*tisinput(3,3);
goldstandardT2=(materialID==1).*tisinput(5,1)+(materialID==2).*tisinput(5,2)+(materialID==3).*tisinput(5,3);
goldstandardM0=(materialID==1).*tisinput(1,1)+(materialID==2).*tisinput(1,2)+(materialID==3).*tisinput(1,3);
goldstandardT1(goldstandardT1==0)=nan; goldstandardT2(goldstandardT2==0)=nan; goldstandardM0(goldstandardM0==0)=nan;

%% Default Acquisition Parameters
flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
TD = [10,10,10,10];%[0.5,0.5,0.5,0.5];          % s

acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];

% Create consistent noise
% noise=normrnd(0,1,[size(materialID,1),size(materialID,2),size(materialID,3),nacq]);

pdxaccel=1;
pdyaccel=1;
pdarg=[pdxaccel,pdyaccel,pdv];
if tconoverride==1
    tconstrain=ttotal;
else
    if pdarg(3)==-1
        subsmplconstrain=ones([size(materialID,1),size(materialID,2)]);
    else
        subsmplconstrain=bart(sprintf('poisson -Y %i -Z %i -y %f -z %f -V %f',size(materialID,1),size(materialID,2),pdarg(1),pdarg(2),pdarg(3)));
    end
    tconstrain=ttotal*60/ceil(sum(subsmplconstrain(:))/100)-TE_T2prep-TDpT2-nacq*Tacq-TDinv; % seconds
end
% 80% elliptical sampling?

% optFA=0;
% if optFA==1
%     B1inhomflag=1;
%     pspacelabels={'flipAngle','TD(1)','TD(2)','TD(3)','TD(4)'};
%     subsmpllabels={};   %{'variance'};
%     %         pinit=[4,1,1,1,1]';
%     pinit=[4;tconstrain/4*ones([4,1])];
%     pAeq=zeros(length(pinit));
%     pAeq(2,2:end)=[1,1,1,1];
%     pbeq=zeros(size(pinit));
%     pbeq(2)=tconstrain;
%     pmin=[0,0,0,0,0]';
%     pmax=[180,2,2,2,2]';
%     findiffrelstep=1.e-6;
%     tolx=1.e-2;%1.e-5;
%     tolfun=1.e-2;%1.e-5;
%     maxiter=500;
% else
%     B1inhomflag=0;
%     pspacelabels={'TD(1)','TD(2)','TD(3)','TD(4)'};
%     subsmpllabels={};   %{'variance'};
%     %         pinit=[4,1,1,1,1]';
%     pinit=tconstrain/4*ones([4,1]);
%     pAeq=zeros(length(pinit));
%     pAeq(2,:)=[1,1,1,1];
%     pbeq=zeros(size(pinit));
%     pbeq(2)=tconstrain;
%     pmin=[0,0,0,0]';
%     pmax=tconstrain*[1,1,1,1]';
%     findiffrelstep=1.e-6;
%     tolx=1.e-2;%1.e-5;
%     tolfun=1.e-2;%1.e-5;
%     maxiter=500;
% end

        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,Mmodel_GM]=qalas1p(tisinput(1,1),tisinput(1,1),tisinput(3,1),tisinput(5,1),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,Mmodel_WM]=qalas1p(tisinput(1,2),tisinput(1,2),tisinput(3,2),tisinput(5,2),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,Mmodel_CSF]=qalas1p(tisinput(1,3),tisinput(1,3),tisinput(3,3),tisinput(5,3),TR,TE_T2prep,flipAngle,nacq,dt);
% std of patient csf = 9.8360; max signal in patient brain = 500; 
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
signu=3.4762E-4;

% Gray matter data
data{1}=1:5; %1:15;
% White matter data
data{2}=Mmodel_GM; %[Mmodel_GM,Mmodel_WM,Mmodel_CSF];
% CSF data
data{3}=signu*ones(size(data{2}));

% % omit data corresponding to outliers (first 4 points)
% data{1} = data{1}(5:end); % x_i
% data{2} = data{2}(5:end); % y_i
% data{3} = data{3}(5:end); % sigma_yi

% convert sigmas to a covariance matrix and reassign to data{3}
C = diag(data{3}.^2);
data{3} = C;

% define nested sampling parameters
Nlive = 100;
Nmcmc = 0; %input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_gaussian;
model = @qalas_model;
prior = {'M0', 'uniform', 0, 2, 'fixed'; ...
    'T1', 'uniform', 0, 5, 'fixed'; ...
    'T2', 'uniform', 0, 1, 'fixed'};
extraparams = {'TR',TR; ...
    'TE_T2prep',TE_T2prep; ...
    'flipAngle',flipAngle; ...
    'nacq',nacq; ...
    'dt',dt}; 

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
    tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

% plot posterior distributions
wp = [1];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [2];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [1 2];
posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});

