%% Measurement model
% clear all
% close all
% format shortg

function [signal_lib,wn_t1_lib,wn_t2_lib,Msize]=MI_GH_QALAS_5acq_3tissue(nparspace)

M0varflag = 1;
T1varflag = 1;
T2varflag = 1;
%% signal model parameters
%%          GM       WM        CSF     Tumor
T1mean = [1200   , 900     ,  4000    , 250  ]./1000; % s
if T1varflag~=0
    T1stdd = [100    ,  100    ,   200    ,  50  ]./1000; % s
else
    T1stdd = [0,0,0,0];
end
T2mean = [100    ,   80    ,  1000    , 250  ]./1000; % s
if T2varflag~=0
    T2stdd = [  5    ,    4    ,    50    ,  50  ]./1000; % s
else
    T2stdd = [0,0,0,0];
end
M0mean = [  0.9  ,    0.9  ,     1.0  ,   1.2];       % relative intensity
if M0varflag~=0
    M0stdd = [   .05 ,     .05 ,      .05 ,    .4];       % relative intensity
else
    M0stdd = [0,0,0,0];
end

TR = 0.005;              % s
TE_T2prep = 0.100;        % s
Tacq = 0.500;             % s
nacq = 5;
flipAngle = 4;          % deg
TDinv=0.100;            % s

switch nparspace
    case 1
        TDpT2 = [0:20:3000]./1000;
        TDpT1 = [0:20:3000]./1000;
    case 2
        TD = repmat(10:300:3000,[nacq-1,1])./1000;
        TD1 = [800:30:1100]./1000;
        TD2 = [500:20:700]./1000;
        TD3 = [500:20:700]./1000;
        TD4 = [600:20:800]./1000;
        TD5 = [800:40:1200]./1000;
        TD = [900:40:1300;
            0:10:100;
            0:10:100;
            0:10:100;
            0:10:100]./1000;
        % TD = repmat(10:20:200,[nacq-1,1])./1000; % s
    case 3
        % flipAngle = 5:1:15;
        flipAngle = 1:.5:10;      % degrees
        TD1 = [0:10:200]./1000;   % s
end

NumQP=5;
for labelindex=1:3
    %% Quadrature Points
    switch labelindex
        case 1
            [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,T1mean(1),T1stdd(1));
            [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,T2mean(1),T2stdd(1));
        case 2
            [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,T1mean(2),T1stdd(2));
            [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,T2mean(2),T2stdd(2));
        case 3
            [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,T1mean(3),T1stdd(3));
            [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,T2mean(3),T2stdd(3));
    end
    % [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,[mean(T1mean(1:2)),T1mean(3)],[mean(T1stdd(1:2)),T1stdd(3)]);
    % [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,[mean(T2mean(1:2)),T2mean(3)],[mean(T2stdd(1:2)),T2stdd(3)]);
    
    [p1,p2,x1,x2]=ndgrid(TDpT2,TDpT1,xm_t1,xm_t2);
    p1=p1(:);
    p2=p2(:);
    x1=x1(:);
    x2=x2(:);
    lqp=length(p1);
    
    parfor qp=1:lqp
        disp(sprintf('Model eval: %d of %d',qp,lqp))
        dt = [0,0,Tacq,p1(qp),0,TDinv,Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp)];
        [~,Mmodel(:,qp)]=qalas1p(M0mean(1),M0mean(1),x1(qp),x2(qp),TR,TE_T2prep,p1(qp),nacq,dt);
    end
    
    signal_lib(:,labelindex)=Mmodel(:);
    wn_t1_lib(:,labelindex)=wn_t1;
    wn_t2_lib(:,labelindex)=wn_t2;
    
end

Msize=[nacq,length(TDpT2),length(TDpT1),size(wn_t1_lib,1)^2];%size(wn_t1_lib,1)^size(wn_t1_lib,2)*size(wn_t2_lib,1)^size(wn_t2_lib,2)];
disp('Saving...')
save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas5acq3tissue_nq%d_siglib.mat',NumQP),'signal_lib','wn_t1_lib','wn_t2_lib','Msize','-v7.3');

end