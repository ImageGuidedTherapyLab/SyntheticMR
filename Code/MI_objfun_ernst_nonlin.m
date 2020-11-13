
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun]=MI_objfun_ernst_nonlin(flipAngle,tisinput,acqparam,signu)

% signu=3.4762E-4;

%% Assign Acquisition Parameters
% Default Parameters
K=acqparam(1);
H=acqparam(2);
TR=acqparam(3);
TE=acqparam(4);

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[~,xn,~,~,wn]=GaussHermiteNDGauss(NumQP,tisinput(1),tisinput(3));
lqp=length(xn{1}(:));
parfor qp=1:lqp
    %     disp(sprintf('Model eval: %d of %d',qp,lqp))
    
    S(qp)=spoiledgre(K,H,flipAngle,TR,TE,xn{1}(qp),tisinput(2));
    
%     [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
%     [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
%     [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
end

[~,xn2,~,~,wn2]=GaussHermiteNDGauss(NumQP,0,signu);
S2=repmat(S,[size(xn2{1},1),1])+repmat(xn2{1},[1,size(S,2)]);

lnterm=log(sum(repmat(wn(:)',[size(xn2{1},1),1]).*S2,2));
pterm=sum(repmat(wn(:)',[size(xn2{1},1),1]).*S2,2);
hz=-sum(wn2.*pterm.*lnterm);

%% Gauss-Hermite Quadrature MI Approximation
% N-D k-space, no averaging
% imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
% for iii=1:length(wn(:))
%     kspace=S(iii);
%     ksr=real(kspace);
%     ksi=imag(kspace);
%     
%     Ezr=Ezr+ksr*wn(iii);
%     Ezi=Ezi+ksi*wn(iii);
%     Sigrr=Sigrr+ksr.^2*wn(iii);
%     Sigii=Sigii+ksi.^2*wn(iii);
%     Sigri=Sigri+ksr.*ksi*wn(iii);
% end
% end

% N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500;
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;

% detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
% detSigz=Sigrr;
% Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);

% Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=hz;%-Hzmu;

% szmi=size(MI);
MIobjfun=-real(MI);
end

function S = spoiledgre(K,H,flipAngle,TR,TE,T1,T2star)
    S = (K.*H.*sind(flipAngle).*(1-exp(-TR./T1)).*exp(-TE./T2star))./(1-exp(-TR./T1).*cosd(flipAngle));
end