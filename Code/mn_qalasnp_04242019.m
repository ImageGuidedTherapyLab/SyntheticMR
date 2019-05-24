
function [logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(TDin,tisinput,acqparam,materialID,Nlive,Nmcmc)

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

%% Optimization Space Acquisition Parameters
% materialID=1;

%% Tissue Properties
% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
% T1mean = [1200,  900, 4000, 1200]./1000; % s
%     T1stdd = [ 100,  100,  200,  150]./1000; % s
% 
% T2mean = [ 100,   80, 1000,  110]./1000; % s
%     T2stdd = [   5,    4,   50,   10]./1000; % s
% 
% M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
%     M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
% 
% overwritecsf=1;
% if overwritecsf==1
%     T1mean = [1200,  900, 900, 1200]./1000; % s
%     T1stdd = [ 100,  100,  100,  150]./1000; % s
%     
%     T2mean = [ 100,   80, 80,  110]./1000; % s
%     T2stdd = [   5,    4,   4,   10]./1000; % s
%     
%     M0mean = [ 0.9,  0.9,  0.9,  0.9];       % relative intensity
%     M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
% end
% 
% tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];

%% Default Acquisition Parameters
% flipAngle = 4;           % deg
% TR = 0.005;              % s
% TE_T2prep = 0.100;       % s
% Tacq = 0.700;            % s
% TDpT2 = TDin(1);             % s
% TDinv = 0.03;            % s
% nacq = 5;
% TD = TDin(2:end); %[10,10,10,10]; %[0.5,0.5,0.5,0.5];          % s
% 
% acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];
% Default Parameters
flipAngle=acqparam(1);
TR=acqparam(2);
TE_T2prep=acqparam(3);
Tacq=acqparam(4);
TDpT2=acqparam(5);
TDinv=acqparam(6);
nacq=acqparam(7);
TD=acqparam(8:6+nacq);

TDpT2=TDin(1);
TD=TDin(2:end);

dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
[~,Mmodel_GM]=qalas1p(tisinput(1,1),tisinput(1,1),tisinput(3,1),tisinput(5,1),TR,TE_T2prep,flipAngle,nacq,dt);
[~,Mmodel_WM]=qalas1p(tisinput(1,2),tisinput(1,2),tisinput(3,2),tisinput(5,2),TR,TE_T2prep,flipAngle,nacq,dt);
[~,Mmodel_CSF]=qalas1p(tisinput(1,3),tisinput(1,3),tisinput(3,3),tisinput(5,3),TR,TE_T2prep,flipAngle,nacq,dt);

% std of patient csf = 9.8360; max signal in patient brain = 500;
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;

% signu=3.4762E-4;
signu=acqparam(end);

% Dummy variables
data{1}=1:15;
% Mmeas data
switch materialID
    case 1
        data{2}=Mmodel_GM;
    case 2
        data{2}=Mmodel_WM;
    case 3
        data{2}=Mmodel_CSF;
end
% data{2}=[Mmodel_GM,Mmodel_WM,Mmodel_CSF];
% Mmeas variance
data{3}=signu*ones(size(data{2}));

% convert sigmas to a covariance matrix and reassign to data{3}
data{3} = data{3}.^2;

% define nested sampling parameters
% Nlive = 500;
% Nmcmc = 0; %input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_tis2img_gaussian_04242019;
model = @qalasnp_model_04242019;
% for ivox=1:numel(materialID)
%     prior((3*ivox-2):(3*ivox),:) = ...
%         {sprintf('M0_%i',ivox), 'uniform', 0, 2, 'fixed'; ...
%         sprintf('T1_%i',ivox), 'uniform', 0, 5, 'fixed'; ...
%         sprintf('T2_%i',ivox), 'uniform', 0, 1, 'fixed'};
% end
switch materialID
    case 1
        prior = {'M0_GM', 'gaussian', tisinput(1,1), tisinput(2,1), ''; ...
        'T1_GM', 'gaussian', tisinput(3,1), tisinput(4,1), ''; ...
        'T2_GM', 'gaussian', tisinput(5,1), tisinput(6,1), ''};
    case 2
        prior = {'M0_WM', 'gaussian', tisinput(1,2), tisinput(2,2), ''; ...
        'T1_WM', 'gaussian', tisinput(3,2), tisinput(4,2), ''; ...
        'T2_WM', 'gaussian', tisinput(5,2), tisinput(6,2), ''};
    case 3
        prior = {'M0_CSF', 'gaussian', tisinput(1,3), tisinput(2,3), ''; ...
        'T1_CSF', 'gaussian', tisinput(3,3), tisinput(4,3), ''; ...
        'T2_CSF', 'gaussian', tisinput(5,3), tisinput(6,3), ''};
end
% prior = {'M0_GM', 'gaussian', tisinput(1,1), tisinput(2,1), ''; ...
%         'T1_GM', 'gaussian', tisinput(3,1), tisinput(4,1), ''; ...
%         'T2_GM', 'gaussian', tisinput(5,1), tisinput(6,1), ''; ...
%         'M0_WM', 'gaussian', tisinput(1,2), tisinput(2,2), ''; ...
%         'T1_WM', 'gaussian', tisinput(3,2), tisinput(4,2), ''; ...
%         'T2_WM', 'gaussian', tisinput(5,2), tisinput(6,2), ''; ...
%         'M0_CSF', 'gaussian', tisinput(1,3), tisinput(2,3), ''; ...
%         'T1_CSF', 'gaussian', tisinput(3,3), tisinput(4,3), ''; ...
%         'T2_CSF', 'gaussian', tisinput(5,3), tisinput(6,3), ''};
extraparams = {'materialID',materialID; ...
    'TR',TR; ...
    'TE_T2prep',TE_T2prep; ...
    'flipAngle',flipAngle; ...
    'nacq',nacq; ...
    'dt',dt};

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
    tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

% plot posterior distributions
plotflag=0;
if plotflag==1
    wp = [1];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [2];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [3];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [1 2];
    posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});
    wp = [4];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [5];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [6];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [4 5];
    posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});
    wp = [7];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [8];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [9];
    posteriors(post_samples, wp, {prior{wp,1}});
    wp = [7 8];
    posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});
end

end