
% MI in just delay time 1 space

% 1000 samples of gray matter

% MI as function of delay time 1

% compare to ->

% measure 95% confidence interval as function of delay time 1

function [varstats,meanstats,medianstats,M0pred,T1pred,T2pred] = syntheticevalvar_np (tislabel,parspace,pslabels,vardecay)

%% Load tissue labels
materialID = tislabel*ones([1000,1]);

%% Signal model parameters
%          GM       WM        CSF     Tumor
T1mean = [1200   ,  900    ,  4000    , 250  ]./1000;  % s
T1stdd = [100    ,  100    ,   200    ,  50  ]./1000; % s
T2mean = [100    ,   80    ,  1000    , 250  ]./1000; % s
T2stdd = [  5    ,    4    ,    50    ,  50  ]./1000; % s
M0mean = [  0.9  ,    0.9  ,     1.0  ,   1.2];       % relative intensity
M0stdd = [   .05 ,     .05 ,      .05 ,    .4];       % relative intensity

TR = 0.005;             % s
TE_T2prep = 0.100;      % s
Tacq = 0.500;           % s
TDpT2 = 1.0;            % s
TD1=TDpT2;
TDpT1 = 0.05;            % s
TD2=TDpT1; TD3=TDpT1; TD4=TDpT1; TD5=TDpT1;
flipAngle = 4;          % deg
TDinv=0.100;            % s
nacq = 5;
nrepeat = 1;

signu=1E-3;

% Load and replace parameter space points
for iii=1:length(pslabels)
    eval(sprintf('%s=parspace(iii)',pslabels{iii}));
end

dt = [0,0,Tacq,TD1,0,TDinv,Tacq,TD2,Tacq,TD3,Tacq,TD4,Tacq,TD5];

%% Create synthetic data
if vardecay~=0
    stdmapT1=normrnd(0,1,size(materialID));
    stdmapT2=normrnd(0,1,size(materialID));
    stdmapM0=normrnd(0,1,size(materialID));
    
    synthdataT1=(materialID==1).*(T1mean(1)+T1stdd(1)*stdmapT1)+(materialID==2).*(T1mean(2)+T1stdd(2)*T1stdd(2)*stdmapT1)+(materialID==3).*(T1mean(3)+T1stdd(3)*stdmapT1);
    synthdataT2=(materialID==1).*(T2mean(1)+T2stdd(1)*stdmapT2)+(materialID==2).*(T2mean(2)+T2stdd(2)*T2stdd(2)*stdmapT2)+(materialID==3).*(T2mean(3)+T2stdd(3)*stdmapT2);
    synthdataM0=(materialID==1).*(M0mean(1)+M0stdd(1)*stdmapM0)+(materialID==2).*(M0mean(2)+M0stdd(2)*M0stdd(2)*stdmapM0)+(materialID==3).*(M0mean(3)+M0stdd(3)*stdmapM0);
    synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
else
    synthdataT1=(materialID==1).*T1mean(1)+(materialID==2).*T1mean(2)+(materialID==3).*T1mean(3);
    synthdataT2=(materialID==1).*T2mean(1)+(materialID==2).*T2mean(2)+(materialID==3).*T2mean(3);
    synthdataM0=(materialID==1).*M0mean(1)+(materialID==2).*M0mean(2)+(materialID==3).*M0mean(3);
    synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
end

%% Create synthetic QALAS measurements
[~,Mmeas]=qalas(synthdataM0,synthdataM0,synthdataT1,synthdataT2,TR,TE_T2prep,flipAngle,nacq,dt);
stdmapmeas=normrnd(0,signu,size(materialID));
Mmeas=Mmeas+stdmapmeas;

%% Reconstruct synthetic QALAS measurements
% Optimization solution for M0 and T1 prediction
xinit=[mean(M0mean(1:3)),mean(T1mean(1:3)),mean(M0mean(1:3))];%xinit=[mean(M0mean(1:3)),mean(T1mean(1:3))];
smeas=size(Mmeas);
Mmeasvec=reshape(Mmeas,[prod(smeas(1:3)),smeas(4:end)]);
mmvsize=size(Mmeasvec,1);
parfor iii=1:size(Mmeasvec,1)
    if sum(isnan(squeeze(Mmeasvec(iii,:))))>0
        M0predvec(iii)=nan;
        T1predvec(iii)=nan;
        T2predvec(iii)=nan;
    else
        xm=fminsearch(@(x) qalasobjfun(x,squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt),xinit);
        M0predvec(iii)=xm(1);
        T1predvec(iii)=xm(2);
        T2predvec(iii)=xm(3);
    end
    fprintf('Element: %d of %d\n',iii,mmvsize)
end
M0pred(:,:,:)=reshape(M0predvec,smeas(1:3));
T1pred(:,:,:)=reshape(T1predvec,smeas(1:3));
T2pred(:,:,:)=reshape(T2predvec,smeas(1:3));

M0varmeas=quantile(M0pred,0.975)-quantile(M0pred,0.025);
T1varmeas=quantile(T1pred,0.975)-quantile(T1pred,0.025);
T2varmeas=quantile(T2pred,0.975)-quantile(T2pred,0.025);
varstats=[M0varmeas,T1varmeas,T2varmeas];

meanstats=[mean(M0pred(:)),mean(T1pred(:)),mean(T2pred(:))];
medianstats=[median(M0pred(:)),median(T1pred(:)),median(T2pred(:))];

end
