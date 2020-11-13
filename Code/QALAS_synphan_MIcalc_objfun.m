
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIsum,MI]=QALAS_synphan_MIcalc_objfun(scanArchivePath,fileflag,paraminput,erodedilateflag)

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
TE_T2prep = .100;
Tacq = TR*102; %0.500;            % s
TDpT2 = 230632/1E6;             % s
TDinv = archive.DownloadData.rdb_hdr_image.ti/1E6; %0.03;            % s
nacq = 5;
TD = [325168,325168,325168,221636]./1E6;          % s
if fileflag==2
    TDpT2=500000/1E6;
    TD = [142592,142592,142592,500000]./1E6;
end

% TD=tdin./1E6*ones([1,4]);

dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];

GERecon('Archive.Close', archive);


%% Tissue Properties
% Load phantom
% load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
% load([scanArchivePath,'.mat'],'realimg');
% evalmask=load_untouch_nii([scanArchivePath,'_segimg.nii.gz']);
% evalmask=evalmask.img;

% Load synthetic phantom
% erodedilateflag=0;
lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
switch erodedilateflag
    case 1
        tmptissue = load_untouch_nii(lfname);
        materialID=int32(tmptissue.img(15:165,20:200,92));
    case 2
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 1 1x1x1vox -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img(15:165,20:200,92));
        materialID(materialID<0)=0;
    case 3
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -erode 1 1x1x1vox -dilate 2 1x1x1vox -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img(15:165,20:200,92));
        materialID(materialID<0)=0;
    case 4
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 2 1x1x1vox -dilate 3 1x1x1vox -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img(15:165,20:200,92));
        materialID(materialID<0)=0;
    case 5
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -erode 1 1x1x1vox -dilate 2 1x1x1vox -dilate 3 1x1x1vox -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img(15:165,20:200,92));
        materialID(materialID<0)=0;
    case 6
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 1 1x1x1vox -erode 2 1x1x1vox -erode 3 1x1x1vox -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img(15:165,20:200,92));
        materialID(materialID<0)=0;
    case 7
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -erode 1 1x1x1vox -erode 2 1x1x1vox -erode 3 1x1x1vox -o resampleimg.nii.gz',lfname));
        tmptissue = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tmptissue.img(15:165,20:200,92));
        materialID(materialID<0)=0;
end
                
% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
T1mean = [1200,  900, 4000, 1200]./1000; % s
T1stdd = [ 100,  100,  200,  150]./1000; % s
T2mean = [ 100,   80, 1000,  110]./1000; % s
T2stdd = [   5,    4,   50,   10]./1000; % s
M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];

if ~isempty(paraminput)
    tisinput(1,1)=paraminput(1);
    tisinput(1,2)=paraminput(1);
    tisinput(3,1)=paraminput(2);
    tisinput(3,2)=paraminput(2);
    tisinput(5,1)=paraminput(3);
    tisinput(5,2)=paraminput(3);
end

overwritecsf=1;
if overwritecsf==1; tisinput(:,3)=tisinput(:,2); end;

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(3,1:3),tisinput(5,1:3)],[tisinput(4,1:3),tisinput(6,1:3)]);
lqp=length(xn{1}(:));
parfor qp=1:lqp
    %     disp(sprintf('Model eval: %d of %d',qp,lqp))
    dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
    [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
end

%% Gauss-Hermite Quadrature MI Approximation
% disp('Performing quadrature...')
    % N-D k-space, no averaging
    nd=ndims(materialID);
    evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
    imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
    materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
    for iii=1:length(wn(:))
        imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:,iii)+(materialIDtemp==2).*Mmodel_WM(:,iii)+(materialIDtemp==3).*Mmodel_CSF(:,iii);
            for jjj=1:size(imspace,1)
                eval(evalstr);
                %                 kspace(jjj,:,:)=fftshift(fftn(squeeze(imspace(jjj,:,:))));
            end
            ksr=real(kspace);
            ksi=imag(kspace);
        
        Ezr=Ezr+ksr*wn(iii);
        Ezi=Ezi+ksi*wn(iii);
        Sigrr=Sigrr+ksr.^2*wn(iii);
        Sigii=Sigii+ksi.^2*wn(iii);
        Sigri=Sigri+ksr.*ksi*wn(iii);
    end

N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500; 
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
signu=3.4762E-4;
% if B1inhomflag==1
    signu=signu*(1+flipAngle/1.2);
% end
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;
MIsum=sum(MI(:));

end