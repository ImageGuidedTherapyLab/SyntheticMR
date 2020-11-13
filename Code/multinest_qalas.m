
function [logZ, nest_samples, post_samples] = multinest_qalas(TR, flipAngle)

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

% clearvars;

tic;
% experimental data to be fitted
necho = [1; 2; 3; 4; 5];
Mmeas = [889; 2134; 3021; 3684; 3654];     
sigmeas = [0.1; 0.1; 0.1; 0.1; 0.1];

data{1} = necho;
data{2} = Mmeas;
data{3} = diag(sigmeas.^2); %covariance matrix

Nlive = 100;
Nmcmc = 0; %input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_gaussian;
model = @mnqalasmodel;
% parameters and range
prior = {'M0', 'uniform', 3500, 3550, 'fixed'; ...
         'T1', 'uniform', 0, .1, 'fixed'; ...
         'T2', 'uniform', 0, 1, 'fixed'};
extraparams = {'TR', TR; ...
               'flipAngle', flipAngle}; % no extra parameters beyond what's in the prior

% use nested_sampler
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

toc;

end