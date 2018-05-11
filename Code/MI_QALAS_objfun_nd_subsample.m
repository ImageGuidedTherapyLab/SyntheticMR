
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun]=MI_QALAS_objfun_nd_subsample(pspace,pspacelabels,subsmpllabels,tisinput,synthdataT1,synthdataT2,synthdataM0,noise,acqparam,materialID,pdarg,B1inhomflag,filename)

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
variance=-1;

% Optimization Parameters
for iii=1:length(pspacelabels)
    eval(sprintf('%s=pspace(iii);',pspacelabels{iii}));
end

% Subsampling Parameters
for iii=length(pspacelabels)+1:length(subsmpllabels)
    eval(sprintf('%s=pspace(iii);',subsmpllabels{iii}));
end

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(3,1:3),tisinput(5,1:3)],[tisinput(4,1:3),tisinput(6,1:3)]);
lqp=length(xn{1}(:));
if B1inhomflag==2
    parfor qp=1:lqp
        %     disp(sprintf('Model eval: %d of %d',qp,lqp))
        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle+flipAngle*randn/(3*1.96),nacq,dt);
        [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle+flipAngle*randn/(3*1.96),nacq,dt);
        [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle+flipAngle*randn/(3*1.96),nacq,dt);
    end
else
    parfor qp=1:lqp
        %     disp(sprintf('Model eval: %d of %d',qp,lqp))
        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    end
end

%% Following Section Possibly More Efficient than Preceeding Section for Larger Quadrature Point Space
% NumQP=5;
% [x_GM,xn_GM,xm_GM,w_GM,wn_GM]=GaussHermiteNDGauss(NumQP,[tisinput(3,1),tisinput(5,1)],[tisinput(4,1),tisinput(6,1)]);
% [x_WM,xn_WM,xm_WM,w_WM,wn_WM]=GaussHermiteNDGauss(NumQP,[tisinput(3,2),tisinput(5,2)],[tisinput(4,2),tisinput(6,2)]);
% [x_CSF,xn_CSF,xm_CSF,w_CSF,wn_CSF]=GaussHermiteNDGauss(NumQP,[tisinput(3,3),tisinput(5,3)],[tisinput(4,3),tisinput(6,3)]);
% lqp=length(xn_GM{1}(:));
% parfor qp=1:lqp
%     fsprintf('Model eval: %d of %d',qp,lqp)
%     dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
%     [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn_GM{1}(qp),xn_GM{2}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
%     [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn_WM{1}(qp),xn_WM{2}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
%     [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn_CSF{1}(qp),xn_CSF{2}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
% end


%% Gauss-Hermite Quadrature MI Approximation
% disp('Performing quadrature...')
% if geometryflag==1
%     % Averaging into one voxel
%     kspace=sum(cat(3,sum(materialID(:)==1)*Mmodel_GM,sum(materialID(:)==2)*Mmodel_WM,sum(materialID(:)==3)*Mmodel_CSF),3)'./numel(materialID);
%     ksr=real(kspace);
%     ksi=imag(kspace);
%     wnmult=repmat(wn(:),[1,size(kspace,2)]);
%     
%     Ezr=sum(ksr.*wnmult,1);
%     Ezi=sum(ksi.*wnmult,1);
%     Sigrr=sum(ksr.^2.*wnmult,1);
%     Sigii=sum(ksi.^2.*wnmult,1);
%     Sigri=sum(ksr.*ksi.*wnmult,1);
% else
    % N-D k-space, no averaging
    nd=ndims(materialID);
    evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
    imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
    materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
    for iii=1:length(wn(:))
        imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:,iii)+(materialIDtemp==2).*Mmodel_WM(:,iii)+(materialIDtemp==3).*Mmodel_CSF(:,iii);
%         if optcase~=0
            for jjj=1:size(imspace,1)
                eval(evalstr);
                %                 kspace(jjj,:,:)=fftshift(fftn(squeeze(imspace(jjj,:,:))));
            end
            ksr=real(kspace);
            ksi=imag(kspace);
%         else
%             ksr=real(imspace);
%             ksi=imag(imspace);
%         end
        
        Ezr=Ezr+ksr*wn(iii);
        Ezi=Ezi+ksi*wn(iii);
        Sigrr=Sigrr+ksr.^2*wn(iii);
        Sigii=Sigii+ksi.^2*wn(iii);
        Sigri=Sigri+ksr.*ksi*wn(iii);
    end
% end

N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500; 
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
signu=3.4762E-4;
if B1inhomflag==1
    signu=signu*(1+flipAngle/1.2);
end
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

% switch optcase
%     case 0
%         % Sum over image space
%         MIobjfun=-sum(MI(:));
%     case 1
%         % Sum over k-space
%         MIobjfun=-sum(MI(:));
%     case 2
%         % Center of k-space
%         szmi=size(MI);
%         centercoord=num2cell(ceil(szmi/2));
%         MIobjfun=-sum(MI(:,sub2ind(szmi(2:end),centercoord{2:end})));
%     case 3
%         % Sum over Poisson disc subsample
%         if variance==-1
%             subsmplmask=0;
%             MIobjfun=-sum(MI(:));
%         else
%             subsmpltype='continuousgauss';
%             switch subsmpltype
%                 case 'importmask'
%                     szmi(2:1+ndims(subsmplmask))=1;
%                     subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%                     subsmplmask=repmat(subsmplmask,szmi);
%                     MIobjfun=-sum(MI(:).*subsmplmask(:)/nseed);
%                 case 'poissondisc'
%                     nseed=100;
%                     szmi=size(MI);
%                     if szmi>=3
%                         poissonmask=zeros(szmi(2:max(3,end-1)));
%                         for iii=1:nseed
%                             poissonpts=poissonDisc(szmi(2:max(3,end-1)),pdspace);
%                             poissonpts=round(poissonpts);
%                             poissonind=zeros([size(poissonpts,1),1]);
%                             for iii=1:size(poissonpts,1)
%                                 tmpcell=num2cell(poissonpts(iii,:));
%                                 poissonind(iii)=sub2ind(szmi(2:max(3,end-1)),tmpcell{:});
%                             end
%                             poissonmasktmp=zeros(szmi(2:max(3,end-1)));
%                             poissonmasktmp(poissonind)=1;
%                             poissonmask=poissonmask+poissonmasktmp;
%                         end
%                         szmi(2:1+ndims(poissonmask))=1;
%                         poissonmask=permute(poissonmask,[ndims(poissonmask)+1,1:ndims(poissonmask)]);
%                         poissonmask=repmat(poissonmask,szmi);
%                         MIobjfun=-sum(MI(:).*poissonmask(:)/nseed);
%                     else
%                         MIobjfun=-sum(MI(:));
%                     end
%                     %         szmi=size(MI);
%                     %         if szmi>=3
%                     %             szmi(2:1+ndims(subsamplemask))=1;
%                     %             subsamplemask=permute(subsamplemask,[ndims(subsamplemask)+1,1:ndims(subsamplemask)]);
%                     %             subsamplemask=repmat(subsamplemask,szmi);
%                     %             MIobjfun=-sum(MI(:).*subsamplemask(:));
%                     %         else
%                     %             MIobjfun=-sum(MI(:));
%                     %         end
%                 case 'discretegauss'
%                     szmi=size(MI);
%                     for iii=1:2%max(2,ndims(MI)-2)
%                         gaussln{iii}=exp(-variance).*besseli(-ceil(szmi(iii+1)/2-1):floor(szmi(iii+1)/2),variance);
%                     end
%                     subsmplmask=kron(gaussln{1}',gaussln{2});
%                     szmi(2:1+ndims(subsmplmask))=1;
%                     subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%                     subsmplmask=repmat(subsmplmask,szmi);
%                     MIobjfun=-sum(MI(:).*subsmplmask(:));
%                 case 'continuousgauss'
%                     szmi=size(MI);
%                     [xg,yg]=ndgrid(-ceil(szmi(2)/2-1):floor(szmi(2)/2),-ceil(szmi(3)/2-1):floor(szmi(3)/2));
%                     sig=diag([variance,variance]);
%                     mu=[0,0];
%                     for iii=1:numel(xg)
%                         gauss2dc(iii)=(det(2*pi*sig))^(-0.5)*exp(-0.5*(mu-[xg(iii),yg(iii)])*inv(sig)*(mu-[xg(iii),yg(iii)])');
%                     end
%                     gauss2dc=reshape(gauss2dc,size(xg));
%                     subsmplmask=gauss2dc/sum(gauss2dc(:));
%                     szmi(2:1+ndims(subsmplmask))=1;
%                     subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%                     subsmplmask=repmat(subsmplmask,szmi);
%                     MIobjfun=-sum(MI(:).*subsmplmask(:));
%             end
%         end
%     case 4
%         szmi=size(MI);
%         subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',szmi(2),szmi(3),1,1,1));
%         subsmplmask=double(squeeze(subsmplmask));
%         szmi(2:1+ndims(subsmplmask))=1;
%         subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%         subsmplmask=repmat(subsmplmask,szmi);
%         MIobjfun=-sum(MI(:).*subsmplmask(:));
%     case 5
%         szmi=size(MI);
%         subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',szmi(2),szmi(3),1,1,3));
%         subsmplmask=double(squeeze(subsmplmask));
%         szmi(2:1+ndims(subsmplmask))=1;
%         subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%         subsmplmask=repmat(subsmplmask,szmi);
%         MIobjfun=-sum(MI(:).*subsmplmask(:));
%     case 6
%         szmi=size(MI);
%         subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',szmi(2),szmi(3),1,1,10));
%         subsmplmask=double(squeeze(subsmplmask));
%         szmi(2:1+ndims(subsmplmask))=1;
%         subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%         subsmplmask=repmat(subsmplmask,szmi);
%         MIobjfun=-sum(MI(:).*subsmplmask(:));
%     case 7
%         szmi=size(MI);
%         subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',szmi(2),szmi(3),1,1,20));
%         subsmplmask=double(squeeze(subsmplmask));
%         szmi(2:1+ndims(subsmplmask))=1;
%         subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%         subsmplmask=repmat(subsmplmask,szmi);
%         MIobjfun=-sum(MI(:).*subsmplmask(:));
%     case 8
%         szmi=size(MI);
%         subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',szmi(2),szmi(3),1,1,25));
%         subsmplmask=double(squeeze(subsmplmask));
%         szmi(2:1+ndims(subsmplmask))=1;
%         subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
%         subsmplmask=repmat(subsmplmask,szmi);
%         MIobjfun=-sum(MI(:).*subsmplmask(:));
%     case 9
%         % Morphological operations on segmentation, 20 patients
%         MIobjfun=-sum(MI(:));
%     otherwise
%         MIobjfun=-sum(MI(:));
% end

szmi=size(MI);
if pdarg(3)==-1
    subsmplmask=ones([1,szmi(2),szmi(3)]);
else
    subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',szmi(2),szmi(3),pdarg(1),pdarg(2),pdarg(3)));
    subsmplmask=double(squeeze(subsmplmask));
end
szmi(2:1+ndims(subsmplmask))=1;
subsmplmask=permute(subsmplmask,[ndims(subsmplmask)+1,1:ndims(subsmplmask)]);
subsmplmask=repmat(subsmplmask,szmi);
MIobjfun=-sum(MI(:).*subsmplmask(:));

%% Reconstruct to Compute Variances
reconflag=1;
if reconflag~=0
%     %% Create synthetic data
%     vardecay=1;
%     if vardecay~=0
%         stdmapT1=normrnd(0,1,size(materialID));
%         stdmapT2=normrnd(0,1,size(materialID));
%         stdmapM0=normrnd(0,1,size(materialID));
%         
%         synthdataT1=(materialID==1).*(tisinput(3,1)+tisinput(4,1)*stdmapT1)+(materialID==2).*(tisinput(3,2)+tisinput(4,2)*tisinput(4,2)*stdmapT1)+(materialID==3).*(tisinput(3,3)+tisinput(4,3)*stdmapT1);
%         synthdataT2=(materialID==1).*(tisinput(5,1)+tisinput(6,1)*stdmapT2)+(materialID==2).*(tisinput(5,2)+tisinput(6,2)*tisinput(6,2)*stdmapT2)+(materialID==3).*(tisinput(5,3)+tisinput(6,3)*stdmapT2);
%         synthdataM0=(materialID==1).*(tisinput(1,1)+tisinput(2,1)*stdmapM0)+(materialID==2).*(tisinput(1,2)+tisinput(2,2)*tisinput(2,2)*stdmapM0)+(materialID==3).*(tisinput(1,3)+tisinput(2,3)*stdmapM0);
%         synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
%     else
%         synthdataT1=(materialID==1).*tisinput(3,1)+(materialID==2).*tisinput(3,2)+(materialID==3).*tisinput(3,3);
%         synthdataT2=(materialID==1).*tisinput(5,1)+(materialID==2).*tisinput(5,2)+(materialID==3).*tisinput(5,3);
%         synthdataM0=(materialID==1).*tisinput(1,1)+(materialID==2).*tisinput(1,2)+(materialID==3).*tisinput(1,3);
%         synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
%     end
    
    %% Create synthetic QALAS measurements
    dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
    [~,Mmeas]=qalas(synthdataM0,synthdataM0,synthdataT1,synthdataT2,TR,TE_T2prep,flipAngle,nacq,dt);
%     stdmapmeas=normrnd(0,signu,size(materialIDtemp)); %materialIDtemp?
%     stdmapmeas=permute(stdmapmeas,[2,3,4,1]);
    stdmapmeas=signu*noise;
    Mmeas=Mmeas+stdmapmeas;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % THE FOLLOWING SECTION ONLY WORKS FOR 2D!
    % FIX IT FOR N-D!!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Subsample synthetic measurements
    Mmeassub=Mmeas;
%     if and(optcase==3,subsmplmask~=0)
%         Mmeassub(isnan(Mmeassub))=0;
%         kmeas=bart('fft 3',Mmeassub);
%         subsample=squeeze(subsmplmask(1,:,:))*10000;%permute(subsmplmask,[2,3,4,1])*100000;
%         subsample=subsample>=rand(size(subsample));
%         subsample=repmat(subsample,[1,1,size(kmeas,3),size(kmeas,4)]);
%         % pics recon - replicate data in coil channel direction and input ones
%         % for sens maps
%         Mmeassub=bart('fft -i 3',kmeas.*subsample)*size(kmeas,4)/numel(kmeas);
%         Mmeassub=real(Mmeassub);
%         Mmeassub(materialID==0)=nan;
%     end
%     if and(optcase>3,optcase<9)
if pdarg(3)~=-1
        Mmeassub(isnan(Mmeassub))=0;
        kmeas=bart('fft 3',Mmeassub);
        subsample=squeeze(subsmplmask(1,:,:));
        subsample=repmat(subsample,[1,1,size(kmeas,3),size(kmeas,4)]);
%         Mmeassub=bart('pics',kmeas.*subsample,ones(size(kmeas)));
        Mmeassub=bart('fft -i 3',kmeas.*subsample)*size(kmeas,4)/numel(kmeas);
        Mmeassub=real(Mmeassub);
        Mmeassub(repmat(materialID,[1,1,size(Mmeassub,3),size(Mmeassub,4)])==0)=nan;
end
%     end
    
    %% Reconstruct synthetic QALAS measurements
    % Optimization solution for M0 and T1 prediction
    xinit=[mean(tisinput(1,1:3)),mean(tisinput(3,1:3))];%,mean(tisinput(5,1:3))];
    smeas=size(Mmeassub);
    Mmeasvec=reshape(Mmeassub,[prod(smeas(1:3)),smeas(4:end)]);
    mmvsize=size(Mmeasvec,1);
    parfor iii=1:size(Mmeasvec,1)
        if sum(isnan(squeeze(Mmeasvec(iii,:))))>0
            M0predvec(iii)=0;%nan;
            T1predvec(iii)=0;%nan;
            T2predvec(iii)=0;%nan;
        else
%             xm=fminsearch(@(x) qalasobjfun(x,squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt),xinit);
%             M0predvec(iii)=xm(1);
%             T1predvec(iii)=xm(2);
%             T2predvec(iii)=qalasT2calc(M0predvec(iii),T1predvec(iii),squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt);
%             T2predvec(iii)=xm(3);
            [M0predvec(iii),T1predvec(iii),T2predvec(iii)=qalasrecon(squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt);
        end
        %         fprintf('Element: %d of %d\n',iii,mmvsize)
    end
    M0pred(:,:,:)=reshape(M0predvec,smeas(1:3));
    T1pred(:,:,:)=reshape(T1predvec,smeas(1:3));
    T2pred(:,:,:)=reshape(T2predvec,smeas(1:3));
    
    M0p1set=M0pred(:).*(materialID(:)==1); M0p1set=M0p1set(M0p1set~=0);
    M0p2set=M0pred(:).*(materialID(:)==2); M0p2set=M0p2set(M0p2set~=0);
    M0p3set=M0pred(:).*(materialID(:)==3); M0p3set=M0p3set(M0p3set~=0);
    T1p1set=T1pred(:).*(materialID(:)==1); T1p1set=T1p1set(T1p1set~=0);
    T1p2set=T1pred(:).*(materialID(:)==2); T1p2set=T1p2set(T1p2set~=0);
    T1p3set=T1pred(:).*(materialID(:)==3); T1p3set=T1p3set(T1p3set~=0);
    T2p1set=T2pred(:).*(materialID(:)==1); T2p1set=T2p1set(T2p1set~=0);
    T2p2set=T2pred(:).*(materialID(:)==2); T2p2set=T2p2set(T2p2set~=0);
    T2p3set=T2pred(:).*(materialID(:)==3); T2p3set=T2p3set(T2p3set~=0);
    
    M0varmeas=[quantile(M0p1set,0.975)-quantile(M0p1set,0.025)...
        quantile(M0p2set,0.975)-quantile(M0p2set,0.025)...
        quantile(M0p3set,0.975)-quantile(M0p3set,0.025)];
    T1varmeas=[quantile(T1p1set,0.975)-quantile(T1p1set,0.025)...
        quantile(T1p2set,0.975)-quantile(T1p2set,0.025)...
        quantile(T1p3set,0.975)-quantile(T1p3set,0.025)];
    T2varmeas=[quantile(T2p1set,0.975)-quantile(T2p1set,0.025)...
        quantile(T2p2set,0.975)-quantile(T2p2set,0.025)...
        quantile(T2p3set,0.975)-quantile(T2p3set,0.025)];
    varstats=[M0varmeas;T1varmeas;T2varmeas];
    
    meanstats=[mean(M0p1set),mean(M0p2set),mean(M0p3set);
        mean(T1p1set),mean(T1p2set),mean(T1p3set);
        mean(T2p1set),mean(T2p2set),mean(T2p3set)];
    
    medianstats=[median(M0p1set),median(M0p2set),median(M0p3set);
        median(T1p1set),median(T1p2set),median(T1p3set);
        median(T2p1set),median(T2p2set),median(T2p3set)];
    
    %% Save statistics
%     filename='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/MI_QALAS_subsample_poptrecons.mat';
    if exist(filename,'file')==2
        load(filename);
        pspacesave=cat(ndims(pspace)+1,pspacesave,pspace);
        MIsave=cat(ndims(MIobjfun)+1,MIsave,MIobjfun);
        varsave=cat(ndims(varstats)+1,varsave,varstats);
        meansave=cat(ndims(meanstats)+1,meansave,meanstats);
        mediansave=cat(ndims(medianstats)+1,mediansave,medianstats);
        M0save=cat(ndims(M0pred)+1,M0save,M0pred);
        T1save=cat(ndims(T1pred)+1,T1save,T1pred);
        T2save=cat(ndims(T2pred)+1,T2save,T2pred);
        Mmeassave=cat(ndims(Mmeas)+1,Mmeassave,Mmeas);
    else
        pspacesave=pspace;
        MIsave=MIobjfun;
        varsave=varstats;
        meansave=meanstats;
        mediansave=medianstats;
        M0save=M0pred;
        T1save=T1pred;
        T2save=T2pred;
        Mmeassave=Mmeas;
    end
    save(filename,'pspacesave','MIsave','varsave','meansave','mediansave','M0save','T1save','T2save','Mmeassave','-v7.3');
end

end