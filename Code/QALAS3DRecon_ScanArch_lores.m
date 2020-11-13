
function[M0pred,T1pred,T2pred] = QALAS3DRecon_ScanArch_lores(scanArchivePath,fileflag)

if nargin==0
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
    fileflag=1;
end
if nargin==1
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

archive = GERecon('Archive.Load', [scanArchivePath,'.h5']);
% realimg=load_untouch_nii([scanArchivePath,'_realimg.nii']);
% realimg=realimg.img;
load([scanArchivePath,'.mat'],'realimg');
% evalmask=load_untouch_nii([scanArchivePath,'_segimg.nii.gz']);
% evalmask=evalmask.img;
load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;

for iii=1:5
    tmpnii=load_untouch_nii([scanArchivePath,'_realimg.nii']);
    tmpnii.img=realimg(:,:,:,iii);
    save_untouch_nii(tmpnii,'tmpnii.nii.gz');
    system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -interpolation NearestNeighbor -resample 50x50x50 -o resampleimg.nii.gz','tmpnii.nii.gz'));
    tmprealimg = load_untouch_nii('resampleimg.nii.gz');
    system('rm resampleimg.nii.gz');
    system('rm tmpnii.nii.gz');
    realimg_lores(:,:,:,iii) = tmprealimg.img;
end
realimg=realimg_lores;

system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -interpolation NearestNeighbor -resample 50x50x50 -o resampleimg.nii.gz',[scanArchivePath,'_segimg.nii.gz']));
tmpevalmask = load_untouch_nii('resampleimg.nii.gz');
system('rm resampleimg.nii.gz');
evalmask = tmpevalmask.img;


% QALAS Fit
% TI = archive.DownloadData.rdb_hdr_image.ti/1E6;                     % s
% TE = archive.DownloadData.rdb_hdr_image.te/1E6;                     % s

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
dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];

GERecon('Archive.Close', archive);

qalasimg=realimg;
qalasimg(evalmask==0)=nan;
[M0pred,T1pred,T2pred]=qalasrecon(qalasimg,3,TR,TE_T2prep,flipAngle,nacq,dt);

% Make figures
M0v=cell([1,42]); T1v=cell([1,42]); T2v=cell([1,42]); M1v=cell([1,42]); M2v=cell([1,42]);
M3v=cell([1,42]); M4v=cell([1,42]); M5v=cell([1,42]);
for iii=1:42
    M0v{iii}(:)=M0pred(evalmask==iii);
    T1v{iii}(:)=T1pred(evalmask==iii);
    T2v{iii}(:)=T2pred(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,1)); M1v{iii}(:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,2)); M2v{iii}(:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,3)); M3v{iii}(:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,4)); M4v{iii}(:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,5)); M5v{iii}(:)=ritmp(evalmask==iii);
end
% M6v=M0v-(M0v-M5v).*exp(-dt(end)./T1v);

% Exclude Gibbs Ringing Voxels
% nstd=1;
% Mv=cat(3,M1v,M2v,M3v,M4v,M5v);
% Mmed=median(Mv,2);
% Mstd=std(Mv,[],2);
% M_ub=repmat(Mmed+nstd*Mstd,[1,size(Mv,2),1]);
% M_lb=repmat(Mmed-nstd*Mstd,[1,size(Mv,2),1]);
% gibbs_mask=(Mv<M_ub)&(Mv>M_lb);
% gibbs_mask=prod(gibbs_mask,3);

% Save results
M0predlr=M0pred;
T1predlr=T1pred;
T2predlr=T2pred;
M0vlr=M0v;
T1vlr=T1v;
T2vlr=T2v;
M1vlr=M1v;
M2vlr=M2v;
M3vlr=M3v;
M4vlr=M4v;
M5vlr=M5v;
save([scanArchivePath,'_parampred_lores.mat'],'M0predlr','T1predlr','T2predlr','M0vlr','T1vlr','T2vlr','M1vlr','M2vlr','M3vlr','M4vlr','M5vlr','-v7.3');%,'M6v','gibbs_mask','-v7.3');
