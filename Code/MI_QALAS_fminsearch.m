
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [signal_lib,wn_t1_lib,wn_t2_lib,Msize,parspace,pslabels]=MI_GH_QALAS_npspace(nparspace)

%% Tissue Properties

% M0/T1/T2 Variance Flags
M0varflag = 1;
T1varflag = 1;
T2varflag = 1;

% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
    T1mean = [1200,  900, 4000,  250]./1000; % s
if T1varflag~=0
    T1stdd = [ 100,  100,  200,   50]./1000; % s
else
    T1stdd = [   0,    0,    0,    0];
end
    T2mean = [ 100,   80, 1000,  250]./1000; % s
if T2varflag~=0
    T2stdd = [   5,    4,   50,   50]./1000; % s
else
    T2stdd = [   0,    0,    0,    0];
end
    M0mean = [ 0.9,  0.9,  1.0,  1.2];       % relative intensity
if M0varflag~=0
    M0stdd = [ .05,  .05,  .05,   .4];       % relative intensity
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
TD=[1,1,1,1,1];          % s

%% Input files
labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation

%% Loading tissue types
disp('loading tissue types');
resized=1.5;
switch resized
    case 1
        materialID = int32(1);
    case 1.5
        materialID = int32([1,2,3]);
    case 2
        tissuelabel = load_untouch_nii(labelfilename);
        materialID=int32(tissuelabel.img(15:165,20:200,92));
    case 2.5
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
        tissuelabel = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tissuelabel.img);
    case 0
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
        tissuelabel = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tissuelabel.img);
    otherwise
        tissuelabel = load_untouch_nii(labelfilename);
        materialID = int32(tissuelabel.img);
end

%%% MI Calc function?
%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[T1mean(1:3),T2mean(1:3)],[T1stdd(1:3),T2stdd(1:3)]);
lqp=length(xn{1}(:));
parfor qp=1:lqp
    disp(sprintf('Model eval: %d of %d',qp,lqp))
    dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
    [~,Mmodel_GM(:,qp)]=qalas1p(M0mean(1),M0mean(1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_WM(:,qp)]=qalas1p(M0mean(2),M0mean(2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_CSF(:,qp)]=qalas1p(M0mean(3),M0mean(3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
end

%% Following Section Possibly More Efficient than Preceeding Section for Larger Quadrature Point Space
% NumQP=5;
% [x_GM,xn_GM,xm_GM,w_GM,wn_GM]=GaussHermiteNDGauss(NumQP,[T1mean(1),T2mean(1)],[T1stdd(1),T2stdd(1)]);
% [x_WM,xn_WM,xm_WM,w_WM,wn_WM]=GaussHermiteNDGauss(NumQP,[T1mean(2),T2mean(2)],[T1stdd(2),T2stdd(2)]);
% [x_CSF,xn_CSF,xm_CSF,w_CSF,wn_CSF]=GaussHermiteNDGauss(NumQP,[T1mean(3),T2mean(3)],[T1stdd(3),T2stdd(3)]);
% lqp=length(xn_GM{1}(:));
% parfor qp=1:lqp
%     fsprintf('Model eval: %d of %d',qp,lqp)
%     dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
%     [~,Mmodel_GM(:,qp)]=qalas1p(M0mean(1),M0mean(1),xn_GM{1}(qp),xn_GM{2}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
%     [~,Mmodel_WM(:,qp)]=qalas1p(M0mean(2),M0mean(2),xn_WM{1}(qp),xn_WM{2}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
%     [~,Mmodel_CSF(:,qp)]=qalas1p(M0mean(3),M0mean(3),xn_CSF{1}(qp),xn_CSF{2}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
% end

disp('Performing quadrature...')
switch resized
    case 1
%         kspace(1,:,:,:,:,:,:,:,:)=reshape(signal_lib(:,1),Msize);
%         ksr=real(kspace);
%         ksi=imag(kspace);
%         wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
%         wn_k=wn_k(:);
%         eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
%         
%         Ezr=sum(ksr.*wnmult,ndims(kspace));
%         Ezi=sum(ksi.*wnmult,ndims(kspace));
%         Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
%         Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
%         Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
    case 1.5
        kspace=mean(cat(3,Mmodel_GM,Mmodel_WM,Mmodel_CSF),3)';
        ksr=real(kspace);
        ksi=imag(kspace);
        wnmult=repmat(wn(:),[1,size(kspace,2)]);
       
        Ezr=sum(ksr.*wnmult,1);
        Ezi=sum(ksi.*wnmult,1);
        Sigrr=sum(ksr.^2.*wnmult,1);
        Sigii=sum(ksi.^2.*wnmult,1);
        Sigri=sum(ksr.*ksi.*wnmult,1);
    case 2
%         siglib=reshape(signal_lib(:,1),Msize);
%         siglib=squeeze(eval(sprintf('siglib(%s,:)',visdim)));
%         ssiglib=size(siglib);
%         siglib=siglib(:);
%         lsiglib=length(siglib);
%         materialID=(materialID~=0);
%         parfor iii=1:lsiglib
%             disp(sprintf('Fourier transform: %d of %d',iii,lsiglib))
%             kspace(:,:,iii) = fftshift(fftn(siglib(iii)*materialID));
%         end
%         kspace=reshape(kspace,[size(kspace,1),size(kspace,2),ssiglib]);
%         ksr=real(kspace);
%         ksi=imag(kspace);
%         wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
%         wn_k=wn_k(:);
%         eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
%         
%         Ezr=sum(ksr.*wnmult,ndims(kspace));
%         Ezi=sum(ksi.*wnmult,ndims(kspace));
%         Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
%         Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
%         Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
    case 2.5
%         wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
%         wn_k=wn_k(:);
%         
%         siglib=reshape(signal_lib,[Msize,size(signal_lib,2)]);
%         siglib=squeeze(eval(sprintf('siglib(%s,:,:)',visdim)));
%         ssiglib=size(siglib);
%         siglib=reshape(siglib,[prod(ssiglib(1:end-2)),ssiglib(end-1),ssiglib(end)]);
% 
%         Ezr=zeros([size(materialID),prod(ssiglib(1:end-2))]);
%         Ezi=Ezr; Sigrr=Ezr; Sigii=Ezr; Sigri=Ezr;
%         sizeijk=ssiglib(end-1);
%         sizeh=prod(ssiglib(1:end-2));
%         for kkk=1:sizeijk
%             for jjj=1:sizeijk
%                 for iii=1:sizeijk
%                     tic;
%                     disp(sprintf('FFT/Quadrature: %d of %d, %d of %d, %d of %d',kkk,ssiglib(end-1),jjj,ssiglib(end-1),iii,ssiglib(end-1)))
%                     for hhh=1:sizeh
%                         disp(sprintf('FFT/Quadrature: %d of %d, %d of %d, %d of %d, %d/%d',kkk,ssiglib(end-1),jjj,ssiglib(end-1),iii,ssiglib(end-1),hhh,prod(ssiglib(1:end-2))))
%                         kspace(:,:,:,hhh) = fftshift(fftn(siglib(hhh,iii,1).*(materialID==1) + siglib(hhh,jjj,2).*(materialID==2) + siglib(hhh,kkk,3).*(materialID==3)));
%                     end
%                     toc;
%                     w=wn_k(iii)*wn_k(jjj)*wn_k(kkk);
%                     Ezr=Ezr+w.*real(kspace);
%                     Ezi=Ezi+w.*imag(kspace);
%                     Sigrr=Sigrr+w.*real(kspace).^2;
%                     Sigii=Sigii+w.*imag(kspace).^2;
%                     Sigri=Sigri+w.*real(kspace).*imag(kspace);
%                 end
%             end
%         end
%         Ezr=reshape(Ezr,[size(materialID),ssiglib(1:end-2)]);
%         Ezi=reshape(Ezi,[size(materialID),ssiglib(1:end-2)]);
%         Sigrr=reshape(Sigrr,[size(materialID),ssiglib(1:end-2)]);
%         Sigii=reshape(Sigii,[size(materialID),ssiglib(1:end-2)]);
%         Sigri=reshape(Sigri,[size(materialID),ssiglib(1:end-2)]);
    otherwise
%         siglib=reshape(signal_lib(:,1),Msize);
%         siglib=squeeze(eval(sprintf('siglib(%s,:)',visdim)));
%         ssiglib=size(siglib);
%         siglib=siglib(:);
%         lsiglib=length(siglib);
%         materialID=(materialID~=0);
%         parfor iii=1:lsiglib
%             disp(sprintf('Fourier transform: %d of %d',iii,lsiglib))
%             kspace(:,:,:,iii) = fftshift(fftn(siglib(iii)*materialID));
%         end
%         kspace=reshape(kspace,[size(kspace,1),size(kspace,2),size(kspace,3),ssiglib]);
%         ksr=real(kspace);
%         ksi=imag(kspace);
%         wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
%         wn_k=wn_k(:);
%         eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
%         
%         Ezr=sum(ksr.*wnmult,ndims(kspace));
%         Ezi=sum(ksi.*wnmult,ndims(kspace));
%         Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
%         Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
%         Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
end

N=2;
signu=1E-4;
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;
