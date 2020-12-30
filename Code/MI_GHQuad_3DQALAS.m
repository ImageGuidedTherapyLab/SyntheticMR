% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

% N = 3 (3 tissue properties, M0, T1, and T2)
% P = 5 (5 measurements from 3D QALAS acquisition)
% Q_1 = ... = Q_N = S_1 = ... = S_N = 5 (number of quadrature points for eta)
% K_1 = ... = K_P = 5 (number of quadrature points for nu)

function [MIobjfun]=MI_GHQuad_3DQALAS(xopt,tisinput,acqparam,materialID)

%% Tissue Parameters
% Tissue parameters define prior distribution of tissue properties (M0, T1,
% T2)

% M0 mean
% tisinput(1)=0.9;

% M0 var
% tisinput(2)=0.05;

% T1 mean
% tisinput(3)=1.200; %s

% T1 var
% tisinput(4)=0.100; %s

% T2 mean
% tisinput(5)=0.100; %s

% T2 var
% tisinput(6)=0.005; %s

%% Acquisition Parameters
% Acquisition parameter variables are read in and renamed from the acqparam
% array.
flipAngle=acqparam(1);
TR=acqparam(2);
TE_T2prep=acqparam(3);
Tacq=acqparam(4);
% TDpT2=acqparam(5);
TDinv=acqparam(6);
nacq=acqparam(7);
% TD=acqparam(8:6+nacq);
signu=acqparam(end);

% Array of delay times is input as a separate variable for the optimization
% function.
TD=xopt;

%% Calculate and Store Quadrature Points and 3D QALAS Evaluations
% Compute quadrature points for eta (omega_q1, ..., omega_qN, omega_s1,
% ..., omega_SN)
% Note that these can be recycled by using the same quadrature points for both
% instances where quadrature is used to integrate over eta. This also
% significantly reduces the number of signal model evaluations required.
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(1:2:5)],[tisinput(2:2:6)]);
lqp=length(xn{1}(:));
% xn{1}(iii) contains the value for eta(1), which is M0, at quadrature
% point iii; xn{2}(iii) contains the value for eta(2), which is T1, at
% quadrature point iii; xn{3}(iii) contains the value for eta(3), which is
% T2, at quadrature point iii.
% wn(iii) contains the product of quadrature weights for quadrature point
% iii; e.g. if wn(iii) corresponds to
% wn(1,2,3)=omega_(q1=1)*omega_(q2=2)*omega_(q3=3).

% Evaluate 3D QALAS model at each quadrature point (G_mu(eta_q),
% G_mu(eta_s))
for qp=1:lqp
    dt=[0,TE_T2prep,Tacq,TD(1),0,TDinv,Tacq,TD(2),Tacq,TD(3),Tacq,TD(4),Tacq,TD(5)];
    [~,Mmodel(:,qp)]=qalas1p(xn{1}(qp),xn{1}(qp),xn{2}(qp),xn{3}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
end

% Compute quadrature points for nu (omega_k1, ..., omega_kP)
NumQP2=5;
[x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NumQP2,zeros([1,5]),signu*ones([1,5]));
lqp2=length(xn2{1}(:));
% xn{1}(iii) contains the value for sqrt(2)*Sigma_nu*x_k(1) at quadrature
% point iii, etc.
% wn(iii) contains the product of quadrature weights for quadrature point
% iii; e.g. if wn(iii) corresponds to
% wn(1,2,3,4,5)=omega_(k1=1)*omega_(k2=2)*omega_(k3=3)*omega_(k4=4)*omega_(k5=5).

%% Compute MI
% Define covariance matrix Sigma_nu
covnu=signu*ones([1,5]);

% Compute the ln term as a function of q1,...,qN QPs and k1,...,kP QPs by
% summing over s1,...,sN QPs
for iii=1:lqp
    for jjj=1:lqp2
        znu=[xn2{1}(jjj),xn2{2}(jjj),xn2{3}(jjj),xn2{4}(jjj),xn2{5}(jjj)];
        lntermtmp=0;
        for kkk=1:lqp
            lntermtmp=lntermtmp+wn(kkk)*mvnpdf(znu+Mmodel(:,iii)',Mmodel(:,kkk)',covnu);
        end
        lnterm(iii,jjj)=log(lntermtmp)+log(pi^(-1.5));
    end
end

% Store all permutations of weighted ln term
Hztmp=wn(:)*wn2(:)'.*lnterm(iii,jjj);
% Compute H(z) by summing over q1,...,qN QPs and k1,...,kP QPs
Hz=-pi^(-1.5-2.5)*sum(Hztmp(:));

% Compute H(z|eta), which is a function of measurement noise
Hzeta=.5*log((2*pi*2.7183)^5.*signu.^5);

% Compute mutual information as the difference in entropy of the
% evidence, H(z), and entropy of the likelihood, H(z|eta).
MI=Hz-Hzeta;

% Make objective function negative, because optimization is minimizing
% objective function and maximizing mutual information.
MIobjfun=-MI;

end
