
% MI objective function for fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun]=MI_objfun_fminsearch_1d(xin)

% FAin=xin(1);
% FAin=4;
TD=xin;

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

%% Acquisition Parameters
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
nacq = 5;
flipAngle = 4;           % deg
TDinv=0.03;              % s
TDpT2=1;                 % s
% TD=[0.4,TDin,TDin,TDin,TDin];          % s

%% Loading tissue types
% disp('loading tissue types');
materialID = int32([1,2,3]);

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[T1mean(1:3),T2mean(1:3)],[T1stdd(1:3),T2stdd(1:3)]);
lqp=length(xn{1}(:));
parfor qp=1:lqp
%     disp(sprintf('Model eval: %d of %d',qp,lqp))
    dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
    [~,Mmodel_GM(:,qp)]=qalas1p(M0mean(1),M0mean(1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_WM(:,qp)]=qalas1p(M0mean(2),M0mean(2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_CSF(:,qp)]=qalas1p(M0mean(3),M0mean(3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
end

% disp('Performing quadrature...')
kspace=mean(cat(3,Mmodel_GM,Mmodel_WM,Mmodel_CSF),3)';
ksr=real(kspace);
ksi=imag(kspace);
wnmult=repmat(wn(:),[1,size(kspace,2)]);

Ezr=sum(ksr.*wnmult,1);
Ezi=sum(ksi.*wnmult,1);
Sigrr=sum(ksr.^2.*wnmult,1);
Sigii=sum(ksi.^2.*wnmult,1);
Sigri=sum(ksr.*ksi.*wnmult,1);

N=2;
signu=1E-4;
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

MIobjfun=-sum(MI);
