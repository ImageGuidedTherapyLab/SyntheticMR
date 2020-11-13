
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun]=MI_QALAS_objfun_nd(pspace,pspacelabels,tisinput,acqparam,materialID,subsamplemask,optcase,geometryflag,B1inhomflag)

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

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(3,1:3),tisinput(5,1:3)],[tisinput(4,1:3),tisinput(6,1:3)]);
lqp=length(xn{1}(:));
if B1inhomflag==2
    parfor qp=1:lqp
        %     disp(sprintf('Model eval: %d of %d',qp,lqp))
        dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle+flipAngle*randn/(3*1.96),nacq,dt);
        [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle+flipAngle*randn/(3*1.96),nacq,dt);
        [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle+flipAngle*randn/(3*1.96),nacq,dt);
    end
else
    parfor qp=1:lqp
        %     disp(sprintf('Model eval: %d of %d',qp,lqp))
        dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
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
if geometryflag==1
    % Averaging into one voxel
    kspace=sum(cat(3,sum(materialID(:)==1)*Mmodel_GM,sum(materialID(:)==2)*Mmodel_WM,sum(materialID(:)==3)*Mmodel_CSF),3)'./numel(materialID);
    ksr=real(kspace);
    ksi=imag(kspace);
    wnmult=repmat(wn(:),[1,size(kspace,2)]);
    
    Ezr=sum(ksr.*wnmult,1);
    Ezi=sum(ksi.*wnmult,1);
    Sigrr=sum(ksr.^2.*wnmult,1);
    Sigii=sum(ksi.^2.*wnmult,1);
    Sigri=sum(ksr.*ksi.*wnmult,1);
else
    % N-D k-space, no averaging
    nd=ndims(materialID);
    evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
    imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
    materialID=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
    for iii=1:length(wn(:))
        imspace=zeros(size(materialID))+(materialID==1).*Mmodel_GM(:,iii)+(materialID==2).*Mmodel_WM(:,iii)+(materialID==3).*Mmodel_CSF(:,iii);
        if optcase~=0
            for jjj=1:size(imspace,1)
                eval(evalstr);
                %                 kspace(jjj,:,:)=fftshift(fftn(squeeze(imspace(jjj,:,:))));
            end
            ksr=real(kspace);
            ksi=imag(kspace);
        else
            ksr=real(imspace);
            ksi=imag(imspace);
        end
        
        Ezr=Ezr+ksr*wn(iii);
        Ezi=Ezi+ksi*wn(iii);
        Sigrr=Sigrr+ksr.^2*wn(iii);
        Sigii=Sigii+ksi.^2*wn(iii);
        Sigri=Sigri+ksr.*ksi*wn(iii);
    end
end

N=length(xn);
signu=1E-4;
if B1inhomflag==1
    signu=signu*(1+flipAngle/1.2);
end
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

switch optcase
    case 0
        % Sum over image space
        MIobjfun=-sum(MI(:));
    case 1
        % Sum over k-space
        MIobjfun=-sum(MI(:));
    case 2
        % Center of k-space
        szmi=size(MI);
        centercoord=num2cell(ceil(szmi/2));
        MIobjfun=-sum(MI(:,sub2ind(szmi(2:end),centercoord{2:end})));
    case 3
        % Sum over Poisson disc subsample
%         nseed=100;
%         szmi=size(MI);
%         if szmi>=3
%             poissonmask=zeros(szmi(2:max(3,end-1)));
%             for iii=1:nseed
%                 poissonpts=poissonDisc(szmi(2:max(3,end-1)),1);
%                 poissonpts=round(poissonpts);
%                 poissonind=zeros([size(poissonpts,1),1]);
%                 for iii=1:size(poissonpts,1)
%                     tmpcell=num2cell(poissonpts(iii,:));
%                     poissonind(iii)=sub2ind(szmi(2:max(3,end-1)),tmpcell{:});
%                 end
%                 poissonmasktmp=zeros(szmi(2:max(3,end-1)));
%                 poissonmasktmp(poissonind)=1;
%                 poissonmask=poissonmask+poissonmasktmp;
%             end
%             szmi(2:1+ndims(poissonmask))=1;
%             poissonmask=permute(poissonmask,[ndims(poissonmask)+1,1:ndims(poissonmask)]);
%             poissonmask=repmat(poissonmask,szmi);
%             MIobjfun=-sum(MI(:).*poissonmask(:)/nseed);
            szmi=size(MI);
        if szmi>=3
            szmi(2:1+ndims(subsamplemask))=1;
            subsamplemask=permute(subsamplemask,[ndims(subsamplemask)+1,1:ndims(subsamplemask)]);
            subsamplemask=repmat(subsamplemask,szmi);
            MIobjfun=-sum(MI(:).*subsamplemask(:));
        else
            MIobjfun=-sum(MI(:));
        end
    case 4
        % Morphological operations on segmentation, 20 patients
        MIobjfun=-sum(MI(:));
    otherwise
        MIobjfun=-sum(MI(:));
end
