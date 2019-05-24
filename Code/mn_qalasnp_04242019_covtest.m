
function [logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,materialID,Nlive,Nmcmc,priorcase)

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

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

signu=3.4762E-4;

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

switch priorcase
    case 'gaussian'
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
        
    case 'uniform'
        switch materialID
            case 1
                prior = {'M0_GM', 'uniform', tisinput(1,1)-2*tisinput(2,1), tisinput(1,1)+2*tisinput(2,1), ''; ...
                    'T1_GM', 'uniform', tisinput(3,1)-2*tisinput(4,1), tisinput(3,1)+2*tisinput(4,1), ''; ...
                    'T2_GM', 'uniform', tisinput(5,1)-2*tisinput(6,1), tisinput(5,1)+2*tisinput(6,1), ''};
            case 2
                prior = {'M0_WM', 'uniform', tisinput(1,2)-2*tisinput(2,2), tisinput(1,2)+2*tisinput(2,2), ''; ...
                    'T1_WM', 'uniform', tisinput(3,2)-2*tisinput(4,2), tisinput(3,2)+2*tisinput(4,2), ''; ...
                    'T2_WM', 'uniform', tisinput(5,2)-2*tisinput(6,2), tisinput(5,2)+2*tisinput(6,2), ''};
            case 3
                prior = {'M0_CSF', 'uniform', tisinput(1,3)-2*tisinput(2,3), tisinput(1,3)+2*tisinput(2,3), ''; ...
                    'T1_CSF', 'uniform', tisinput(3,3)-2*tisinput(4,3), tisinput(3,3)+2*tisinput(4,3), ''; ...
                    'T2_CSF', 'uniform', tisinput(5,3)-2*tisinput(6,3), tisinput(5,3)+2*tisinput(6,3), ''};
        end
        
    case 'cov'
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
end

extraparams = {'materialID',materialID; ...
    'TR',TR; ...
    'TE_T2prep',TE_T2prep; ...
    'flipAngle',flipAngle; ...
    'nacq',nacq; ...
    'dt',dt};

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
    tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

end