function logL = logL_tis2img_gaussian(data, model, parnames, parvals)

% logL = logL_gaussian(data, model, parnames, parvals)
%
% This function will compute the log likelihood of a multivariate
% gaussian:
%
%     L = 1/sqrt((2 pi)^N det C)
%         exp[-0.5*(y - model(x,params))^T * inv(C) * (y - model(x,params))]
%
% The input parameters are:
%     data - a cell array with three columns
%            { x values, y values, C: covariance matrix }
%     NOTE: if C is a single number, convert to a diag covariance matrix
%     model - the function handle for the signal model.
%     parnames - a cell array listing the names of the model parameters
%     parvals - a cell array containing the values of the parameters given
%         in parnames. These must be in the same order as in parnames.
%         If parvals is an empty vector the noise-only likelihood will be
%         calculated.
%
% -------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
% -------------------------------------------------------------------------

% check whether model is a string or function handle
if ischar(model)
    fmodel = str2func(model);
elseif isa(model, 'function_handle')
    fmodel = model;
else
    error('Error... Expecting a model function!');
end

% get image structure
lpn = length(parnames);
for ii=1:lpn
    switch parnames{ii}
        case 'materialID'
            materialID = parvals{ii};
        case 'subsamplemask'
            subsamplemask = parvals{ii};
    end
end

% get data values from cell array
x = data{1};
y = data{2};
C = data{3};
N = length(x);

% evaluate the model
if isempty(parvals)
    % if parvals is not defined get the null likelihood (noise model
    % likelihood)
    md = 0;
else
    md = feval(fmodel, x, parnames, parvals);

    % if the model returns a NaN then set the likelihood to be zero (e.g.
    % loglikelihood to be -inf
    if isnan(md)
        logL = -inf;
        return;
    end
end

% get inverse of covariance matrix and the log of the determinant
invC = C.^(-1);
% calculate the log likelihood
C_GM=C(1:5);
C_WM=C(6:10);
C_CSF=C(11:15);
invC_GM=invC(1:5);
invC_WM=invC(6:10);
invC_CSF=invC(11:15);

% ymd_GM=y(1:5)-md(1:5);
% ymd_WM=y(6:10)-md(6:10);
% ymd_CSF=y(11:15)-md(11:15);
% tic;
y=kron(materialID(:)==1,y(1:5))+kron(materialID(:)==2,y(6:10))+kron(materialID(:)==3,y(11:15));
y=reshape(y,[size(materialID),5]);
y=permute(y,[ndims(y),1:ndims(y)-1]);
y=bart(sprintf('fft %i',sum(2.^(1:ndims(y)-1))),y);
y=y.*subsamplemask;
% toc;

% tic;
md=kron(materialID(:)==1,md(1:5))+kron(materialID(:)==2,md(6:10))+kron(materialID(:)==3,md(11:15));
md=reshape(md,[size(materialID),5]);
md=permute(md,[ndims(md),1:ndims(md)-1]);
md=bart(sprintf('fft %i',sum(2.^(1:ndims(md)-1))),md);
md=md.*subsamplemask;
% toc;

C_tis=kron(ones(size(materialID(:))),C_GM);%+kron(materialID(:)==2,C_GM)+kron(materialID(:)==3,C_CSF);
invC_tis=kron(ones(size(materialID(:))),invC_GM);%+kron(materialID(:)==2,invC_GM)+kron(materialID(:)==3,invC_CSF);
% ymd_tis=kron(materialID(:)==1,ymd_GM)+kron(materialID(:)==2,ymd_WM)+kron(materialID(:)==3,ymd_CSF);
C_tis=[C_tis(:);C_tis(:)];
invC_tis=[invC_tis(:);invC_tis(:)];
yc=zeros([2*length(y(:)),1]);
yc(1:2:end)=real(y(:));
yc(2:2:end)=imag(y(:));
mdc=zeros([2*length(md(:)),1]);
mdc(1:2:end)=real(md(:));
mdc(2:2:end)=imag(md(:));
ymd_tis=yc(:)-mdc(:);
% ymd_tis=ymd_tis(:);

logL = -0.5*ymd_tis' .* invC_tis' * ymd_tis;

lDetC = log(prod(C_tis));
% calculate the log likelihood
% logL = logL - 0.5*length(invC_tis)*log(2*pi) - 0.5*lDetC;
logL = logL - 0.5*length(invC_tis)*log(2*pi);% - 0.5*lDetC;

if isnan(logL)
    error('Error: log likelihood is NaN!');
end

return
