
% Conditional MI Plot Figs

%% Processing Flags
nind=[3:14,17:28,31:42];
gibbs_exclude=0;
t1correction=1;
t1correctionfull=1;
MIcalcflag=2;

%% Scan 1 Data
scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
fileflag=1;

% High Res
load([scanArchivePath,'_parampred.mat']);

if gibbs_exclude~=0
    Mv(~gibbs_mask)=nan;
    M0v(~gibbs_mask)=nan;
    T1v(~gibbs_mask)=nan;
    T2v(~gibbs_mask)=nan;
    M5v(~gibbs_mask)=nan;
end
if t1correction~=0
    t1manexclude(1:14)=[.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1];
    t1manexclude(15:28)=[0,0,0,0,0,.1,.1,.1,.1,0,0,0,.01,.01];
    t1manexclude(29:42)=[0,0,0,0,0,0,0,0,0,.1,.1,.12,.1,.12];
    for iii=1:42
        t1tmp=T1v(iii,:);
        %         t1tmp(t1tmp<t1manexclude(iii))=nan;
        T1v(iii,t1tmp<t1manexclude(iii))=nan;
        if t1correctionfull~=0
            %             m0tmp=M0v(iii,:);
            %             m0tmp(t1tmp<t1manexclude(iii))=nan;
            M0v(iii,t1tmp<t1manexclude(iii))=nan;
            %             t2tmp=T2v(iii,:);
            %             t2tmp(t1tmp<t1manexclude(iii))=nan;
            T2v(iii,t1tmp<t1manexclude(iii))=nan;
        end
    end
    %     T1v(T1v<.1)=nan;
end
M0v_1=M0v; T1v_1=T1v; T2v_1=T2v; %signuin1=nanvar(M5v(15:end,:),[],2)';
M0std_1=nanstd(M0v_1,[],2);
T1std_1=nanstd(T1v_1,[],2);
T2std_1=nanstd(T2v_1,[],2);

% Low Res
load([scanArchivePath,'_parampred_lores.mat']);
M0predlr_1=M0predlr; M0vlr_1=M0vlr; M1vlr_1=M1vlr; M2vlr_1=M2vlr; M3vlr_1=M3vlr; M4vlr_1=M4vlr; M5vlr_1=M5vlr;
T1predlr_1=T1predlr; T1vlr_1=T1vlr; T2predlr_1=T2predlr; T2vlr_1=T2vlr;

M0stdlr_1=cellfun(@nanstd,M0vlr_1);
T1stdlr_1=cellfun(@nanstd,T1vlr_1);
T2stdlr_1=cellfun(@nanstd,T2vlr_1);
M0Hlr_1=log(sqrt(2*pi*exp(1))*M0stdlr_1);
T1Hlr_1=log(sqrt(2*pi*exp(1))*T1stdlr_1);
T2Hlr_1=log(sqrt(2*pi*exp(1))*T2stdlr_1);

% Man Val
load([scanArchivePath,'_manvalMI.mat']);
MIman_1=sum(MI,2);

% load([scanArchivePath,'_parampred_lores.mat']);
% load([scanArchivePath,'_parampred.mat']);


%% Scan 2 Data
scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
fileflag=2;

% High Res
load([scanArchivePath,'_parampred.mat']);

if gibbs_exclude~=0
    Mv(~gibbs_mask)=nan;
    M0v(~gibbs_mask)=nan;
    T1v(~gibbs_mask)=nan;
    T2v(~gibbs_mask)=nan;
    M5v(~gibbs_mask)=nan;
end
if t1correction~=0
    t1manexclude=[0,0,0,0,0,.1,.1,.075,.1,.1,0,0,.01,.01];
    for iii=15:28
        t1tmp=T1v(iii,:);
        t1tmp(t1tmp<t1manexclude(iii-14))=nan;
        T1v(iii,:)=t1tmp;
        %         T1v(T1v<t1manexclude(iii-14))=nan;
    end
    %     T1v(T1v<.1)=nan;
end
M0v_2=M0v; T1v_2=T1v; T2v_2=T2v; %signuin2=nanvar(M5v(15:end,:),[],2)';
M0std_2=nanstd(M0v_2,[],2);
T1std_2=nanstd(T1v_2,[],2);
T2std_2=nanstd(T2v_2,[],2);

% Low Res
load([scanArchivePath,'_parampred_lores.mat']);
M0predlr_2=M0predlr; M0vlr_2=M0vlr; M1vlr_2=M1vlr; M2vlr_2=M2vlr; M3vlr_2=M3vlr; M4vlr_2=M4vlr; M5vlr_2=M5vlr;
T1predlr_2=T1predlr; T1vlr_2=T1vlr; T2predlr_2=T2predlr; T2vlr_2=T2vlr;

M0stdlr_2=cellfun(@nanstd,M0vlr_2);
T1stdlr_2=cellfun(@nanstd,T1vlr_2);
T2stdlr_2=cellfun(@nanstd,T2vlr_2);
M0Hlr_2=log(sqrt(2*pi*exp(1))*M0stdlr_2);
T1Hlr_2=log(sqrt(2*pi*exp(1))*T1stdlr_2);
T2Hlr_2=log(sqrt(2*pi*exp(1))*T2stdlr_2);


% Man Val
load([scanArchivePath,'_manvalMI.mat']);
MIman_2=sum(MI,2);


%% Phantom Properties
load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
% load([scanArchivePath,'_manvalMI.mat']);
% load([scanArchivePath,'_parampred_lores.mat']);
% load([scanArchivePath,'_parampred.mat']);

% MIman=sum(MI,2);
% M0stdlr=cellfun(@nanstd,M0vlr);
% T1stdlr=cellfun(@nanstd,T1vlr);
% T2stdlr=cellfun(@nanstd,T2vlr);
% M0Hlr=log(sqrt(2*pi*exp(1))*M0stdlr);
% T1Hlr=log(sqrt(2*pi*exp(1))*T1stdlr);
% T2Hlr=log(sqrt(2*pi*exp(1))*T2stdlr);
% M0std=nanstd(M0v,[],2);
% T1std=nanstd(T1v,[],2);
% T2std=nanstd(T2v,[],2);

M0=nanmedian(M0v(15:end,:),2)';
T1=phantomProp.T1element.T30.T1./1000; T2=phantomProp.T1element.T30.T2./1000;
T1=[T1,phantomProp.T2element.T30.T1./1000]; T2=[T2,phantomProp.T2element.T30.T2./1000];


%% High Res MI
for iii=1:42
    htmp=M0v_1(iii,:); htmp(isnan(htmp))=[];
    M0MI1(iii)=0.5.*size(M0v_1,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    M0H1(iii)=entropy(htmp-nanmean(htmp)+.5);
    
    htmp=M0v_2(iii,:); htmp(isnan(htmp))=[];
    M0MI2(iii)=0.5.*size(M0v_2,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    M0H2(iii)=entropy(htmp-nanmean(htmp)+.5);
    
    htmp=T1v_1(iii,:); htmp(isnan(htmp))=[];
    T1MI1(iii)=0.5.*size(T1v_1,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    T1H1(iii)=entropy(htmp-nanmean(htmp)+.5);
    
    htmp=T1v_2(iii,:); htmp(isnan(htmp))=[];
    T1MI2(iii)=0.5.*size(T1v_2,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(htmp-nanmean(htmp)+.5);
    T1H2(iii)=entropy(htmp-nanmean(htmp)+.5);
    
    htmp=T2v_1(iii,:); htmp(isnan(htmp))=[];
    T2MI1(iii)=0.5.*size(T2v_1,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(real(htmp)-nanmean(real(htmp))+.5);
    T2H1(iii)=entropy(real(htmp)-nanmean(real(htmp))+.5);
    
    htmp=T2v_2(iii,:); htmp(isnan(htmp))=[];
    T2MI2(iii)=0.5.*size(T2v_2,2).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)))-entropy(real(htmp)-nanmean(real(htmp))+.5);
    T2H2(iii)=entropy(real(htmp)-nanmean(real(htmp))+.5);    
end


%% Plot Figures


figure; plot(T1MI1(nind),T1std_1(nind),'bo');
hold on; plot(T1MI1(nind),T1stdlr_1(nind),'rx');

figure; plot(T1MI1(nind),T1stdlr_1(nind)-T1std_1(nind)','bo');

figure; plot(M0MI1(nind),M0stdlr_1(nind)-M0std_1(nind)','bo');
figure; plot(T2MI1(nind),T2stdlr_1(nind)-T2std_1(nind)','bo');



figure; plot(T1MI1(nind),T1H1(nind),'bo')
hold on; plot(T1MI1(nind),T1Hlr_1(nind),'rx');
legend('High Res','Low Res');


