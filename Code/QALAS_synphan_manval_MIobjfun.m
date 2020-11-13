
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MI]=QALAS_synphan_manval_MIobjfun(scanArchivePath,fileflag)

if scanArchivePath==1
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
    fileflag=1;
elseif scanArchivePath==2
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
    fileflag=2;
end

if isempty(scanArchivePath)
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
end
if isempty(fileflag)
    fileflag=1;
end

load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
M0phan=ones([1,42]);
T1phan=[ones([1,14]),phantomProp.T1element.T30.T1,phantomProp.T2element.T30.T1];
T2phan=[ones([1,14]),phantomProp.T1element.T30.T2,phantomProp.T2element.T30.T2];

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
materialID=1;

%% Generate Quadrature Points for MI Calculation
NumQP=5;
% [x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(1,1),tisinput(3,1),tisinput(5,1)],[tisinput(2,1),tisinput(4,1),tisinput(6,1)]);
% [x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,paraminput(1:3),paraminput(4:6));

for iii=1:42
    
    [x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[M0phan(iii),T1phan(iii),T2phan(iii)],[0.1,0.1,0.005]);
    
    [M,Mmeas]=qalas(xn{1},xn{1},xn{2},xn{3},TR,TE_T2prep,flipAngle,nacq,dt);
    
%     switch tisid
%         case 1
%             [M,Mmeas]=qalas(xn{1},xn{1},repmat(paraminput(2),size(xn{1})),repmat(paraminput(3),size(xn{1})),TR,TE_T2prep,flipAngle,nacq,dt);
%         case 2
%             [M,Mmeas]=qalas(repmat(paraminput(1),size(xn{1})),repmat(paraminput(1),size(xn{1})),xn{1},repmat(paraminput(3),size(xn{1})),TR,TE_T2prep,flipAngle,nacq,dt);
%         case 3
%             [M,Mmeas]=qalas(repmat(paraminput(1),size(xn{1})),repmat(paraminput(1),size(xn{1})),repmat(paraminput(2),size(xn{1})),xn{1},TR,TE_T2prep,flipAngle,nacq,dt);
%     end
    
    Mmeasvec=reshape(Mmeas,[prod(size(wn)),size(Mmeas,ndims(Mmeas))]);
    Ezr=sum(repmat(wn(:),[1,size(Mmeas,ndims(Mmeas))]).*Mmeasvec,1);
    Ezi=0;
    Sigrr=sum(repmat(wn(:),[1,size(Mmeas,ndims(Mmeas))]).*Mmeasvec.^2,1);
    Sigii=0;
    Sigri=0;
    
    %% Gauss-Hermite Quadrature MI Approximation
    N=length(xn);
    % std of patient csf = 9.8360; max signal in patient brain = 500;
    % max approx signal in synthdata = 0.0584
    % std of noise in patient raw data = 17.8574; max signal approx 3000;
    signu=3.4762E-4;
    % if B1inhomflag==1
    %     signu=signu*(1+flipAngle/1.2);
    % end
%     signu=paraminput(7);
    detSigz=(pi^(-N/2)*signu.^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu.^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
    % detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2);
    Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
    Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
    MI(iii,:)=Hz-Hzmu;
    
    % MI=0.5.*log((2*pi*2.7183)^2.*(pi^(-N/2)*signu.^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu.^2))-0.5.*log((2*pi*2.7183)^2.*signu.^4);
    % MI=0.5.*log((pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2)/signu.^4);
    
    % MI=sum(MI(:));
    
end

save([scanArchivePath,'_manvalMI.mat'],'MI','-v7.3');

end