
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MI]=QALAS_synphan_MIcalc1p_objfun(scanArchivePath,fileflag,paraminput)

if nargin==0
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
    fileflag=1;
end
if nargin==1 || nargin==3
    if scanArchivePath==1
        scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
        fileflag=1;
    elseif scanArchivePath==2
        scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
        fileflag=2;
    else
        scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
        fileflag=1;
    end
end
if isempty(scanArchivePath)
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
end
if isempty(fileflag)
    fileflag=1;
end

%% Assign Acquisition Parameters
archive = GERecon('Archive.Load', [scanArchivePath,'.h5']);

flipAngle = archive.DownloadData.rdb_hdr_image.mr_flip; %4;           % deg
TR = archive.DownloadData.rdb_hdr_image.tr/1E6; %0.005;              % s
% TE_T2prep = archive.DownloadData.rdb_hdr_image.t2PrepTE/1E6; %0.100;       % s
TE_T2prep = 0.100;
Tacq = TR*102; %0.500;            % s
TDpT2 = 230632/1E6;             % s
TDinv = archive.DownloadData.rdb_hdr_image.ti/1E6; %0.03;            % s
nacq = 5;
TD = [325168,325168,325168,221636]./1E6;          % s
if fileflag==2
    TDpT2=500000/1E6;
    TD = [142592,142592,142592,500000]./1E6;
end
dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];

GERecon('Archive.Close', archive);


%% Tissue Properties
% Load phantom
% load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
% load([scanArchivePath,'.mat'],'realimg');
% evalmask=load_untouch_nii([scanArchivePath,'_segimg.nii.gz']);
% evalmask=evalmask.img;

% Load synthetic phantom
% lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
% tmptissue = load_untouch_nii(lfname);
% materialID=int32(tmptissue.img(15:165,20:200,92));
materialID=1;
% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
% T1mean = [1200,  900, 4000, 1200]./1000; % s
% T1stdd = [ 100,  100,  200,  150]./1000; % s
% T2mean = [ 100,   80, 1000,  110]./1000; % s
% T2stdd = [   5,    4,   50,   10]./1000; % s
% M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
% M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
% tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];
% 
% if ~isempty(paraminput)
%     tisinput(1,1)=paraminput(1);
%     tisinput(1,2)=paraminput(1);
%     tisinput(3,1)=paraminput(2);
%     tisinput(3,2)=paraminput(2);
%     tisinput(5,1)=paraminput(3);
%     tisinput(5,2)=paraminput(3);
% end

% overwritecsf=1;
% if overwritecsf==1; tisinput(:,3)=tisinput(:,2); end;

%% Generate Quadrature Points for MI Calculation
NumQP=5;
% [x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(1,1),tisinput(3,1),tisinput(5,1)],[tisinput(2,1),tisinput(4,1),tisinput(6,1)]);
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,paraminput(1:3),paraminput(4:6));
% [M,Mmeas]=qalas(xn{1},xn{1},xn{2},xn{3},TR,TE_T2prep,flipAngle,nacq,dt);
[M,Mmeas]=qalas(xn{1},xn{1},xn{2},xn{3},TR,TE_T2prep,flipAngle,nacq,dt);

Mmeasvec=reshape(Mmeas,[prod(size(wn)),size(Mmeas,ndims(Mmeas))]);
Ezr=sum(repmat(wn(:),[1,size(Mmeas,ndims(Mmeas))]).*Mmeasvec,1);
Ezi=0;
Sigrr=sum(repmat(wn(:),[1,size(Mmeas,ndims(Mmeas))]).*Mmeasvec.^2,1);
Sigii=0;
Sigri=0;

%% Gauss-Hermite Quadrature MI Approximation
% disp('Performing quadrature...')
    % N-D k-space, no averaging
%     nd=ndims(materialID);
%     evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
%     imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
%     materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
%     for iii=1:length(wn(:))
%         imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:,iii)+(materialIDtemp==2).*Mmodel_WM(:,iii)+(materialIDtemp==3).*Mmodel_CSF(:,iii);
%             for jjj=1:size(imspace,1)
%                 eval(evalstr);
%                 %                 kspace(jjj,:,:)=fftshift(fftn(squeeze(imspace(jjj,:,:))));
%             end
%             ksr=real(kspace);
%             ksi=imag(kspace);
%         
%         Ezr=Ezr+ksr*wn(iii);
%         Ezi=Ezi+ksi*wn(iii);
%         Sigrr=Sigrr+ksr.^2*wn(iii);
%         Sigii=Sigii+ksi.^2*wn(iii);
%         Sigri=Sigri+ksr.*ksi*wn(iii);
%     end

N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500; 
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
% signu=3.4762E-4;
% if B1inhomflag==1
%     signu=signu*(1+flipAngle/1.2);
% end
signu=paraminput(7);
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

% MI=sum(MI(:));

end