
function [predstats,M0pred,T1pred,T2pred] = createsyntheticdata (TDpT2,TDpT1,vardecay)

%% Load tissue labels
labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
tissuelabel = load_untouch_nii(labelfilename);
materialID = int32(tissuelabel.img);
materialID = materialID(:,:,89:91);

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
flipAngle = 4;          % deg
TDinv=0.100;            % s
nacq = 5;
nrepeat = 1;

signu=1E-3;
% signu=mean(abs(Mmeas(~isnan(Mmeas))))/50;

% TDpT2 = 0;%1;%.5;   % s
% TDpT1 = 0;%3.5;%.2;   % s
dt = [0,0,Tacq,TDpT2,0,TDinv,Tacq,TDpT1,Tacq,TDpT1,Tacq,TDpT1,Tacq,0];

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
% for lll=1:nrepeat
    % T2 prediction    
%     T2pred(:,:,:,lll)=TE_T2prep./log(Mmeas(:,:,:,end,lll)./Mmeas(:,:,:,1,lll));
    
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
% end

% Optimization averages
% M0pred=mean(M0pred,4);
% T1pred=mean(T1pred,4);
% T2pred=mean(T2pred,4);

T1predmean=[mean(T1pred(materialID==1)),mean(T1pred(materialID==2)),mean(T1pred(materialID==3))];
T1predstd=[std(T1pred(materialID==1)),std(T1pred(materialID==2)),std(T1pred(materialID==3))];
T2predmean=[mean(T2pred(materialID==1)),mean(T2pred(materialID==2)),mean(T2pred(materialID==3))];
T2predstd=[std(T2pred(materialID==1)),std(T2pred(materialID==2)),std(T2pred(materialID==3))];
M0predmean=[mean(M0pred(materialID==1)),mean(M0pred(materialID==2)),mean(M0pred(materialID==3))];
M0predstd=[std(M0pred(materialID==1)),std(M0pred(materialID==2)),std(M0pred(materialID==3))];
predstats=[T1predmean;T1predstd;T2predmean;T2predstd;M0predmean;M0predstd];

end
