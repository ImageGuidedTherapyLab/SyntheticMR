function[M0pred,T1pred,T2pred] = QALAS3DRecon_ScanArch(scanArchivePath,fileflag)

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
evalmask=load_untouch_nii([scanArchivePath,'_segimg.nii.gz']);
evalmask=evalmask.img;
load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;

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
for iii=1:42
    M0v(iii,:)=M0pred(evalmask==iii);
    T1v(iii,:)=T1pred(evalmask==iii);
    T2v(iii,:)=T2pred(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,1)); M1v(iii,:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,2)); M2v(iii,:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,3)); M3v(iii,:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,4)); M4v(iii,:)=ritmp(evalmask==iii);
    ritmp=squeeze(realimg(:,:,:,5)); M5v(iii,:)=ritmp(evalmask==iii);
end
M6v=M0v-(M0v-M5v).*exp(-dt(end)./T1v);

% Exclude Gibbs Ringing Voxels
nstd=1;
Mv=cat(3,M1v,M2v,M3v,M4v,M5v);
Mmed=median(Mv,2);
Mstd=std(Mv,[],2);
M_ub=repmat(Mmed+nstd*Mstd,[1,size(Mv,2),1]);
M_lb=repmat(Mmed-nstd*Mstd,[1,size(Mv,2),1]);
gibbs_mask=(Mv<M_ub)&(Mv>M_lb);
gibbs_mask=prod(gibbs_mask,3);

% Save results
save([scanArchivePath,'_parampred.mat'],'M0pred','T1pred','T2pred','M0v','T1v','T2v','M1v','M2v','M3v','M4v','M5v','M6v','gibbs_mask','-v7.3');

gibbs_exclude=1;
if gibbs_exclude~=0
    Mv(~gibbs_mask)=nan;
    M0v(~gibbs_mask)=nan;
    T1v(~gibbs_mask)=nan;
    T2v(~gibbs_mask)=nan;
end

% PD Element Property Distributions
figure;
subplot(2,3,1);
boxplot(M0v(1:14,:)','Positions',1:14);
xlabel('PD Element Number'); ylabel('M0'); title('PD Element M0'); pause(1);
subplot(2,3,2);
boxplot(T1v(1:14,:)','Positions',1:14);
xlabel('PD Element Number'); ylabel('T1 (s)'); title('PD Element T1'); pause(1);
subplot(2,3,3);
boxplot(real(T2v(1:14,:))','Positions',1:14);
xlabel('PD Element Number'); ylabel('T2 (s)'); title('PD Element T2'); pause(1);
axis([0,15,0,2]);

subplot(2,3,4);
boxplot(M0v(1:14,:)'./nanmean(M0v(1:14,:)',1),'Positions',1:14);
xlabel('PD Element Number'); ylabel('M0 Relative Variance'); title('PD Element M0'); pause(1);
subplot(2,3,5);
boxplot(T1v(1:14,:)'./nanmean(T1v(1:14,:)',1),'Positions',1:14);
xlabel('PD Element Number'); ylabel('T1 Relative Variance'); title('PD Element T1'); pause(1);
subplot(2,3,6);
boxplot(real(T2v(1:14,:))'./nanmean(real(T2v(1:14,:))',1),'Positions',1:14);
xlabel('PD Element Number'); ylabel('T2 Relative Variance'); title('PD Element T2'); pause(1);
axis([0,15,0,2]);

% T1 Element Property Distributions
figure;
subplot(2,3,1);
boxplot(M0v(15:28,:)','Positions',1:14);
xlabel('T1 Element Number'); ylabel('M0'); title('T1 Element M0'); pause(1);
subplot(2,3,2); hold on;
boxplot(T1v(15:28,:)','Positions',1:14);
plot(phantomProp.T1element.T30.T1./1000,'gx');
xlabel('T1 Element Number'); ylabel('T1 (s)'); title('T1 Element T1'); pause(1);
subplot(2,3,3); hold on;
boxplot(real(T2v(15:28,:))','Positions',1:14);
plot(phantomProp.T1element.T30.T2./1000,'gx');
xlabel('T1 Element Number'); ylabel('T2 (s)'); title('T1 Element T2'); pause(1);
axis([0,15,0,2]);

subplot(2,3,4);
boxplot(M0v(15:28,:)'./nanmean(M0v(15:28,:)',1),'Positions',1:14);
xlabel('T1 Element Number'); ylabel('M0 Relative Variance'); title('T1 Element M0'); pause(1);
subplot(2,3,5);
boxplot(T1v(15:28,:)'./nanmean(T1v(15:28,:)',1),'Positions',1:14);
xlabel('T1 Element Number'); ylabel('T1 Relative Variance'); title('T1 Element T1'); pause(1);
subplot(2,3,6);
boxplot(real(T2v(15:28,:))'./nanmean(real(T2v(15:28,:))',1),'Positions',1:14);
xlabel('T1 Element Number'); ylabel('T2 Relative Variance'); title('T1 Element T2'); pause(1);
axis([0,15,0,2]);

% T2 Element Property Distributions
figure;
subplot(2,3,1);
boxplot(M0v(29:42,:)','Positions',1:14);
xlabel('T2 Element Number'); ylabel('M0'); title('T2 Element M0'); pause(1);
subplot(2,3,2); hold on;
boxplot(T1v(29:42,:)','Positions',1:14);
plot(phantomProp.T2element.T30.T1./1000,'gx');
xlabel('T2 Element Number'); ylabel('T1 (s)'); title('T2 Element T1'); pause(1);
subplot(2,3,3); hold on;
boxplot(real(T2v(29:42,:))','Positions',1:14);
plot(phantomProp.T2element.T30.T2./1000,'gx');
xlabel('T2 Element Number'); ylabel('T2 (s)'); title('T2 Element T2'); pause(1);
axis([0,15,0,2]);

subplot(2,3,4);
boxplot(M0v(29:42,:)'./nanmean(M0v(29:42,:)',1),'Positions',1:14);
xlabel('T2 Element Number'); ylabel('M0 Relative Variance'); title('T2 Element M0'); pause(1);
subplot(2,3,5);
boxplot(T1v(29:42,:)'./nanmean(T1v(29:42,:)',1),'Positions',1:14);
xlabel('T2 Element Number'); ylabel('T1 Relative Variance'); title('T2 Element T1'); pause(1);
subplot(2,3,6);
boxplot(real(T2v(29:42,:))'./nanmean(real(T2v(29:42,:))',1),'Positions',1:14);
xlabel('T2 Element Number'); ylabel('T2 Relative Variance'); title('T2 Element T2'); pause(1);
axis([0,15,0,2]);

% All Element Property Distributions
figure('pos',[10,10,1410,910]);
subplot(3,3,1);
boxplot(M0v(1:14,:)','Positions',1:14);
xlabel('PD Element Number'); ylabel('M0'); title('PD Element M0'); pause(1);
subplot(3,3,2);
boxplot(T1v(1:14,:)','Positions',1:14);
xlabel('PD Element Number'); ylabel('T1 (s)'); title('PD Element T1'); pause(1);
subplot(3,3,3);
boxplot(real(T2v(1:14,:))','Positions',1:14);
xlabel('PD Element Number'); ylabel('T2 (s)'); title('PD Element T2'); pause(1);
axis([0,15,0,2]);

subplot(3,3,4);
boxplot(M0v(15:28,:)','Positions',1:14);
xlabel('T1 Element Number'); ylabel('M0'); title('T1 Element M0'); pause(1);
subplot(3,3,5); hold on;
boxplot(T1v(15:28,:)','Positions',1:14);
plot(phantomProp.T1element.T30.T1./1000,'gx');
legend('True Value');
xlabel('T1 Element Number'); ylabel('T1 (s)'); title('T1 Element T1'); pause(1);
subplot(3,3,6); hold on;
boxplot(real(T2v(15:28,:))','Positions',1:14);
plot(phantomProp.T1element.T30.T2./1000,'gx');
legend('True Value');
xlabel('T1 Element Number'); ylabel('T2 (s)'); title('T1 Element T2'); pause(1);
axis([0,15,0,2]);

subplot(3,3,7);
boxplot(M0v(29:42,:)','Positions',1:14);
xlabel('T2 Element Number'); ylabel('M0'); title('T2 Element M0'); pause(1);
subplot(3,3,8); hold on;
boxplot(T1v(29:42,:)','Positions',1:14);
plot(phantomProp.T2element.T30.T1./1000,'gx');
legend('True Value');
xlabel('T2 Element Number'); ylabel('T1 (s)'); title('T2 Element T1'); pause(1);
subplot(3,3,9); hold on;
boxplot(real(T2v(29:42,:))','Positions',1:14);
plot(phantomProp.T2element.T30.T2./1000,'gx');
legend('True Value');
xlabel('T2 Element Number'); ylabel('T2 (s)'); title('T2 Element T2'); pause(1);
axis([0,15,0,2]);


% Plot Measurement Distributions
csdt=cumsum(dt);
figure;
for iii=1:14
    subplot(2,7,iii); boxplot(cat(1,M1v(iii,:),M2v(iii,:),M3v(iii,:),M4v(iii,:),M5v(iii,:),M6v(iii,:))','Positions',csdt([2,6:2:end]));%({'M1','M2','M3','M4','M5'});
    xticklabels(strtrim(cellstr(num2str(csdt([2,6:2:end-1]),'%2.1f\n'))));
    title(['PD Elem ',num2str(iii)]);
end
figure;
for iii=15:28
    subplot(2,7,iii-14); boxplot(cat(1,M1v(iii,:),M2v(iii,:),M3v(iii,:),M4v(iii,:),M5v(iii,:),M6v(iii,:))','Positions',csdt([2,6:2:end]));
    xticklabels(strtrim(cellstr(num2str(csdt([2,6:2:end-1]),'%2.1f\n'))));
    title(['TD1 Elem ',num2str(iii-14)]);
end
figure;
for iii=29:42
    subplot(2,7,iii-28); boxplot(cat(1,M1v(iii,:),M2v(iii,:),M3v(iii,:),M4v(iii,:),M5v(iii,:),M6v(iii,:))','Positions',csdt([2,6:2:end]));
    xticklabels(strtrim(cellstr(num2str(csdt([2,6:2:end-1]),'%2.1f\n'))));
    title(['TD2 Elem ',num2str(iii-28)]);
end
% Normalized
figure;
for iii=1:14
    subplot(2,7,iii); boxplot(cat(1,M1v(iii,:)./nanmedian(M5v(iii,:)),M2v(iii,:)./nanmedian(M5v(iii,:)),M3v(iii,:)./nanmedian(M5v(iii,:)),M4v(iii,:)./nanmedian(M5v(iii,:)),M5v(iii,:)./nanmedian(M5v(iii,:)),M6v(iii,:)./nanmedian(M5v(iii,:)))','Positions',csdt([2,6:2:end]));%({'M1','M2','M3','M4','M5'});
    xticklabels(strtrim(cellstr(num2str(csdt([2,6:2:end-1]),'%2.1f\n'))));
    title(['PD Elem ',num2str(iii)]);
end
figure;
for iii=15:28
    subplot(2,7,iii-14); boxplot(cat(1,M1v(iii,:)./nanmedian(M5v(iii,:)),M2v(iii,:)./nanmedian(M5v(iii,:)),M3v(iii,:)./nanmedian(M5v(iii,:)),M4v(iii,:)./nanmedian(M5v(iii,:)),M5v(iii,:)./nanmedian(M5v(iii,:)),M6v(iii,:)./nanmedian(M5v(iii,:)))','Positions',csdt([2,6:2:end]));
    xticklabels(strtrim(cellstr(num2str(csdt([2,6:2:end-1]),'%2.1f\n'))));
    title(['TD1 Elem ',num2str(iii-14)]);
end
figure;
for iii=29:42
    subplot(2,7,iii-28); boxplot(cat(1,M1v(iii,:)./nanmedian(M5v(iii,:)),M2v(iii,:)./nanmedian(M5v(iii,:)),M3v(iii,:)./nanmedian(M5v(iii,:)),M4v(iii,:)./nanmedian(M5v(iii,:)),M5v(iii,:)./nanmedian(M5v(iii,:)),M6v(iii,:)./nanmedian(M5v(iii,:)))','Positions',csdt([2,6:2:end]));
    xticklabels(strtrim(cellstr(num2str(csdt([2,6:2:end-1]),'%2.1f\n'))));
    title(['TD2 Elem ',num2str(iii-28)]);
end

end