
function [logZ,nest_samples,post_samples,prior]=cond_mn_qalasnp(materialID,tisinput,TDin,noisein,Nlive,Nmcmc)

% in: materialID, tisinput

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

% pdv=-1;
% tconoverride=1;
% ttotal=2;

%% Default Acquisition Parameters
flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
TD = TDin*ones([1,4]); %[10,10,10,10]; %[0.5,0.5,0.5,0.5];          % s

acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];

% Create consistent noise
% noise=normrnd(0,1,[size(materialID,1),size(materialID,2),size(materialID,3),nacq]);

% pdxaccel=1;
% pdyaccel=1;
% pdarg=[pdxaccel,pdyaccel,pdv];
subsamplemask=ones([nacq,size(materialID,1),size(materialID,2)]);

% if tconoverride==1
%     tconstrain=ttotal;
% else
%     tconstrain=ttotal*60/ceil(sum(subsamplemask(:))/100)-TE_T2prep-TDpT2-nacq*Tacq-TDinv; % seconds
% end
% 80% elliptical sampling?

dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
[~,Mmodel_GM]=qalas1p(tisinput(1,1),tisinput(1,1),tisinput(3,1),tisinput(5,1),TR,TE_T2prep,flipAngle,nacq,dt);
[~,Mmodel_WM]=qalas1p(tisinput(1,2),tisinput(1,2),tisinput(3,2),tisinput(5,2),TR,TE_T2prep,flipAngle,nacq,dt);
[~,Mmodel_CSF]=qalas1p(tisinput(1,3),tisinput(1,3),tisinput(3,3),tisinput(5,3),TR,TE_T2prep,flipAngle,nacq,dt);
% Mmodel=repmat(permute(materialID,[ndims(materialID)+1,1:ndims(materialID)]),[5,1]);
% Mmodel=(Mmodel==1)*
for ii=1:length(Mmodel_GM)
    Mmodel{ii}=(materialID==1)*Mmodel_GM(ii)+(materialID==2)*Mmodel_WM(ii)+(materialID==3)*Mmodel_CSF(ii);
    kmodel{ii}=bart(sprintf('fft %i',sum(0:(ndims(materialID)-1))),Mmodel{ii});
end
% std of patient csf = 9.8360; max signal in patient brain = 500;
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
% signu=3.4762E-4;
signu=noisein;

% Dummy variables
data{1}=1:15;
% Mmeas data
data{2}=[Mmodel_GM,Mmodel_WM,Mmodel_CSF];
% Mmeas variance
data{3}=signu*ones(size(data{2}));

% convert sigmas to a covariance matrix and reassign to data{3}
data{3} = data{3}.^2;

% define nested sampling parameters
% Nlive = 500;
% Nmcmc = 0; %input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_tis2img_gaussian;
model = @qalasnp_model;
% for ivox=1:numel(materialID)
%     prior((3*ivox-2):(3*ivox),:) = ...
%         {sprintf('M0_%i',ivox), 'uniform', 0, 2, 'fixed'; ...
%         sprintf('T1_%i',ivox), 'uniform', 0, 5, 'fixed'; ...
%         sprintf('T2_%i',ivox), 'uniform', 0, 1, 'fixed'};
% end
prior = {'M0_GM', 'gaussian', tisinput(1,1), tisinput(2,1), ''; ...
        'T1_GM', 'gaussian', tisinput(3,1), tisinput(4,1), ''; ...
        'T2_GM', 'gaussian', tisinput(5,1), tisinput(6,1), ''; ...
        'M0_WM', 'gaussian', tisinput(1,2), tisinput(2,2), ''; ...
        'T1_WM', 'gaussian', tisinput(3,2), tisinput(4,2), ''; ...
        'T2_WM', 'gaussian', tisinput(5,2), tisinput(6,2), ''; ...
        'M0_CSF', 'gaussian', tisinput(1,3), tisinput(2,3), ''; ...
        'T1_CSF', 'gaussian', tisinput(3,3), tisinput(4,3), ''; ...
        'T2_CSF', 'gaussian', tisinput(5,3), tisinput(6,3), ''};
extraparams = {'materialID',materialID; ...
    'subsamplemask',subsamplemask; ...
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