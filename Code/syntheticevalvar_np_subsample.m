
% MI in just delay time 1 space

% 1000 samples of gray matter

% MI as function of delay time 1

% compare to ->

% measure 95% confidence interval as function of delay time 1

function [varstats,meanstats,medianstats,M0pred,T1pred,T2pred] = syntheticevalvar_np_subsample(pspace,pspacelabels,subsmpllabels,tisinput,acqparam,materialID,subsmplmask)

%% Assign Acquisition Parameters
% Default Parameters
flipAngle=acqparam(1);
TR=acqparam(2);
TE_T2prep=acqparam(3);
Tacq=acqparam(4);
TDpT2=acqparam(5);
TDinv=acqparam(6);
nacq=acqparam(7);
TD=acqparam(8:6+nacq);

% Optimization Parameters
for iii=1:length(pspacelabels)
    eval(sprintf('%s=pspace(iii);',pspacelabels{iii}));
end

% Subsampling Parameters
for iii=length(pspacelabels)+1:length(subsmpllabels)
    eval(sprintf('%s=pspace(iii);',subsmpllabels{iii}));
end

signu=1E-3;

dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];

%% Create synthetic data
vardecay=1;
if vardecay~=0
    stdmapT1=normrnd(0,1,size(materialID));
    stdmapT2=normrnd(0,1,size(materialID));
    stdmapM0=normrnd(0,1,size(materialID));
    
    synthdataT1=(materialID==1).*(tisinput(3,1)+tisinput(4,1)*stdmapT1)+(materialID==2).*(tisinput(3,2)+tisinput(4,2)*tisinput(4,2)*stdmapT1)+(materialID==3).*(tisinput(3,3)+tisinput(4,3)*stdmapT1);
    synthdataT2=(materialID==1).*(tisinput(5,1)+tisinput(6,1)*stdmapT2)+(materialID==2).*(tisinput(5,2)+tisinput(6,2)*tisinput(6,2)*stdmapT2)+(materialID==3).*(tisinput(5,3)+tisinput(6,3)*stdmapT2);
    synthdataM0=(materialID==1).*(tisinput(1,1)+tisinput(2,1)*stdmapM0)+(materialID==2).*(tisinput(1,2)+tisinput(2,2)*tisinput(2,2)*stdmapM0)+(materialID==3).*(tisinput(1,3)+tisinput(2,3)*stdmapM0);
    synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
else
    synthdataT1=(materialID==1).*tisinput(3,1)+(materialID==2).*tisinput(3,2)+(materialID==3).*tisinput(3,3);
    synthdataT2=(materialID==1).*tisinput(5,1)+(materialID==2).*tisinput(5,2)+(materialID==3).*tisinput(5,3);
    synthdataM0=(materialID==1).*tisinput(1,1)+(materialID==2).*tisinput(1,2)+(materialID==3).*tisinput(1,3);
    synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
end

%% Create synthetic QALAS measurements
[~,Mmeas]=qalas(synthdataM0,synthdataM0,synthdataT1,synthdataT2,TR,TE_T2prep,flipAngle,nacq,dt);
stdmapmeas=normrnd(0,signu,size(materialID));
Mmeas=Mmeas+stdmapmeas;

%% Subsample synthetic measurements
Mmeassub=Mmeas;
if subsmplmask~=0
    Mmeassub(isnan(Mmeassub))=0;
    % kmeas=zeros(size(Mmeassub));
    % for iii=1:size(Mmeassub,4)
    %     kmeas(:,:,:,iii)=fftshift(fftn(Mmeassub(:,:,:,iii)));
    % end
    kmeas=bart('fft 3',Mmeassub);
    subsample=squeeze(subsmplmask(1,:,:))*10000;%permute(subsmplmask,[2,3,4,1])*100000;
    subsample=subsample>=rand(size(subsample));
    subsample=repmat(subsample,[1,1,size(kmeas,3),size(kmeas,4)]);
    % pics recon - replicate data in coil channel direction and input ones
    % for sens maps
    Mmeassub=bart('fft -i 3',kmeas.*subsample)*size(kmeas,4)/numel(kmeas);
    Mmeassub=real(Mmeassub);
    Mmeassub(materialID==0)=nan;
end
% for iii=1:size(kmeas,4)
%     Mmeassub(:,:,:,iii)=ifftn(kmeas(:,:,:,iii).*subsample(:,:,:,1));
% end

%% Reconstruct synthetic QALAS measurements
% Optimization solution for M0 and T1 prediction
xinit=[mean(tisinput(1,1:3)),mean(tisinput(3,1:3)),mean(tisinput(1,1:3))];%xinit=[mean(tisinput(1,1:3)),mean(tisinput(3,1:3))];
smeas=size(Mmeassub);
Mmeasvec=reshape(Mmeassub,[prod(smeas(1:3)),smeas(4:end)]);
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
