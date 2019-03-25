
% Conditional MI Plot Figs

scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
fileflag=1;

load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
load([scanArchivePath,'_manvalMI.mat']);
load([scanArchivePath,'_parampred_lores.mat']);
load([scanArchivePath,'_parampred.mat']);

MIman=sum(MI,2);
M0stdlr=cellfun(@std,M0vlr);
T1stdlr=cellfun(@std,T1vlr);
T2stdlr=cellfun(@std,T2vlr);
M0Hlr=log(sqrt(2*pi*exp(1))*M0stdlr);
T1Hlr=log(sqrt(2*pi*exp(1))*T1stdlr);
T2Hlr=log(sqrt(2*pi*exp(1))*T2stdlr);
M0std=std(M0v,[],2);
T1std=std(T1v,[],2);
T2std=std(T2v,[],2);

M0=nanmedian(M0v(15:end,:),2)';
T1=phantomProp.T1element.T30.T1./1000; T2=phantomProp.T1element.T30.T2./1000;
T1=[T1,phantomProp.T2element.T30.T1./1000]; T2=[T2,phantomProp.T2element.T30.T2./1000];

for iii=1:42
    htmp=M0v(iii,:); htmp(isnan(htmp))=[];
    M0MI(iii)=0.5.*size(M0v,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    M0H(iii)=entropy(htmp-nanmean(htmp)+.5);
    %     htmp=M0v_2(iii,:); htmp(isnan(htmp))=[];
%     M0MI2(iii)=0.5.*size(M0v_2,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    htmp=T1v(iii,:); htmp(isnan(htmp))=[];
    T1MI(iii)=0.5.*size(T1v,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    T1H(iii)=entropy(htmp-nanmean(htmp)+.5);
    %     htmp=T1v_2(iii,:); htmp(isnan(htmp))=[];
%     T1MI2(iii)=0.5.*size(T1v_2,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    htmp=T2v(iii,:); htmp(isnan(htmp))=[];
    T2MI(iii)=0.5.*size(T2v,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(real(htmp)-nanmean(real(htmp))+.5);
    T2H(iii)=entropy(real(htmp)-nanmean(real(htmp))+.5);
    %     htmp=T2v_2(iii,:); htmp(isnan(htmp))=[];
%     T2MI2(iii)=0.5.*size(T2v_2,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(real(htmp)-nanmean(real(htmp))+.5);
end











MI
[MI]=QALAS_synphan_manval_MIobjfun(1,[]);
MI
sum(MI,2)
figure;plot(ans)
load([scanArchivePath,'_parampred_lores.mat']);
size(M0pred)
size(M0v)
M0v
cellfun(@std,M0v)
figure;plot(ans)
cellfun(@std,T1v)
figure;plot(ans)
figure;plot(ans,sum(MI,2))
figure;plot(ans,sum(MI,2),'x')
figure;plot(ans(15:28),sum(MI(15:28,:),2),'x')
load([scanArchivePath,'_parampred.mat']);
figure;plot(std(T1v(15:28,:),[],2),sum(MI(15:28,:),2))
figure;plot(std(T1v(15:28,:),[],2),sum(MI(15:28,:),2),'x')
figure;plot(-sum(MI(15:28,:),2),std(T1v(15:28,:),[],2),'x')
figure;plot(sum(MI(15:28,:),2),std(T1v(15:28,:),[],2),'x')
figure;loglog(sum(MI(15:28,:),2),std(T1v(15:28,:),[],2),'x')
figure;plot(sum(MI(29:42,:),2),std(T1v(29:42,:),[],2),'x')
figure;plot(sum(MI(29:42,:),2),std(T2v(29:42,:),[],2),'x')
figure;plot(sum(MI(14:28,:),2),std(T2v(14:28,:),[],2),'x')