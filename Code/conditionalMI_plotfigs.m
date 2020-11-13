
% Conditional MI Plot Figs

%% Processing Flags
nind=[3:14,17:28,31:42];
nind2=nind;
nind2(20:24)=[];
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
% M0Hlr_1=log((2*pi*exp(1)).^(12/2)*M0stdlr_1);
% T1Hlr_1=log((2*pi*exp(1)).^(12/2)*T1stdlr_1);
% T2Hlr_1=log((2*pi*exp(1)).^(12/2)*T2stdlr_1);

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
% M0Hlr_2=log((2*pi*exp(1)).^(12/2)*M0stdlr_2);
% T1Hlr_2=log((2*pi*exp(1)).^(12/2)*T1stdlr_2);
% T2Hlr_2=log((2*pi*exp(1)).^(12/2)*T2stdlr_2);

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


% figure; plot(T1MI1(nind),T1std_1(nind),'bo');
% hold on; plot(T1MI1(nind),T1stdlr_1(nind),'rx');
% 
% figure; plot(T1MI1(nind),T1stdlr_1(nind)-T1std_1(nind)','bo');
% figure; plot(T1MI2(nind),T1stdlr_2(nind)-T1std_2(nind)','bo');
% 
% figure; plot(M0MI1(nind),M0stdlr_1(nind)-M0std_1(nind)','bo');
% figure; plot(T2MI1(nind),T2stdlr_1(nind)-T2std_1(nind)','bo');

% figure; plot(T1MI1(nind),T1H1(nind),'bo')
% hold on; plot(T1MI1(nind),T1Hlr_1(nind),'rx');
% legend('High Res','Low Res');
% 
% figure; plot(T1MI1(nind),T1std_1(nind),'bo')
% hold on; plot(T1MI1(nind),T1stdlr_1(nind),'rx');
% legend('High Res','Low Res');

% figure; plot(T1MI1(nind),T1H1(nind)-T1Hlr_1(nind),'bo');

% Conditional MI
Hzmulr=0.5.*(5+length(T1vlr_1{1})).*log(2*pi*2.7183*(0.05*M0(3)*T1(3)*T2(3)));

% figure;
% plot(Hzmulr-M0Hlr_1(nind)-M0H1(nind),M0std_1(nind),'bo');
% % hold on;
% % plot(Hzmulr-M0Hlr_2(nind)-M0H2(nind),M0std_2(nind),'rx');
% xlabel('Conditional Mutual Information'); ylabel('M0 Standard Deviation');
% 
% figure;
% plot(Hzmulr-T1Hlr_1(nind)-T1H1(nind),T1std_1(nind),'bo');
% % hold on;
% % plot(Hzmulr-T1Hlr_2(nind)-T1H2(nind),T1std_2(nind),'rx');
% xlabel('Conditional Mutual Information'); ylabel('T1 Standard Deviation (s)');
% 
% figure;
% plot(Hzmulr-T2Hlr_1(nind)-T2H1(nind),T2std_1(nind),'bo');
% % hold on;
% % plot(Hzmulr-T2Hlr_2(nind)-T1H2(nind),T2std_2(nind),'rx');
% xlabel('Conditional Mutual Information'); ylabel('T2 Standard Deviation (s)');


figure; hold on;
x=Hzmulr+T1Hlr_1(nind2(1))-T1H1(nind2);
y=T1std_1(nind2)-T1stdlr_1(nind2(1));
plot(Hzmulr+T1Hlr_1(nind2(1))-T1H1(nind2),T1std_1(nind2)-T1stdlr_1(nind2(1)),'bo');
xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev.');
% title('Prediction from Element 3 Properties');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(0,-0.03,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
saveas(gcf,'Figures/Pred_el3','png');

figure; hold on;
x=Hzmulr+T1Hlr_1(nind2(11))-T1H1(nind2);
y=T1std_1(nind2)-T1stdlr_1(nind2(11));
plot(Hzmulr+T1Hlr_1(nind2(11))-T1H1(nind2),T1std_1(nind2)-T1stdlr_1(nind2(11)),'bo');
xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev.');
% title('Prediction from Element 13 Properties');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(0,-0.03,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
saveas(gcf,'Figures/Pred_el13','png');

figure; hold on;
x=Hzmulr+T1Hlr_1(nind2(21))-T1H1(nind2);
y=T1std_1(nind2)-T1stdlr_1(nind2(21));
plot(Hzmulr+T1Hlr_1(nind2(21))-T1H1(nind2),T1std_1(nind2)-T1stdlr_1(nind2(21)),'bo');
f=@(b,x) b(1).*exp(b(2).*x)+b(3);
[B,R,J,COVB,MSE]=nlinfit(x',y,@(b,x) b(1).*exp(b(2).*x)+b(3),[1;-1;0]);
[~,fitsortind]=sort(x);
plot(x(fitsortind),B(1).*exp(B(2).*x(fitsortind))+B(3),'color','black');

[lfit,gob]=fit(x',y,'exp1','StartPoint',[0,0.1]);
plot(x,lfit(1)*exp(lfit(2)*x));
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit);
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-0.7,-0.03,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev.');
% title('Prediction from Element 10 Properties');
saveas(gcf,'Figures/Pred_el10','png');

%%
nind2=[6:10,20:23,34:39];
% figure; hold on;
% [~,sortind]=sort(M0std_1(nind2));
% for iii=1:length(nind2)
%     plot(Hzmulr+M0Hlr_1(nind2(iii))-M0H1(nind2(sortind)),M0std_1(nind2(sortind))-M0stdlr_1(nind2(iii)),'-o');
% end
% xlabel('Conditional Mutual Information'); ylabel('Delta M0 Std. Dev.');
% % title('Prediction between Elements');
% saveas(gcf,'Figures/Pred_els_M0_lines','png');
figure; hold on;
[~,sortind]=sort(M0std_1(nind2));
for iii=1:length(nind2)
    plot(Hzmulr+M0Hlr_1(nind2(iii))-M0H1(nind2(sortind)),M0std_1(nind2(sortind))-M0stdlr_1(nind2(iii)),'bo');
end
xlabel('Conditional Mutual Information'); ylabel('Delta M0 Std. Dev.');

x=[];y=[];
for iii=1:length(nind2)
    x=[x,Hzmulr+M0Hlr_1(nind2(iii))-M0H1(nind2(sortind))];
    y=[y;M0std_1(nind2(sortind))-M0stdlr_1(nind2(iii))];
end
[lfit,gof]=robustfit(x',y);
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-2.8,-0.15,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));

% title('Prediction between Elements');
saveas(gcf,'Figures/Pred_els_M0','png');

% figure; hold on;
% [~,sortind]=sort(T1std_1(nind2));
% x=[]; y=[];
% for iii=1:length(nind2)
%     plot(Hzmulr+T1Hlr_1(nind2(iii))-T1H1(nind2(sortind)),T1std_1(nind2(sortind))-T1stdlr_1(nind2(iii)),'-o');
%     x=[x,Hzmulr+T1Hlr_1(nind2(iii))-T1H1(nind2(sortind))];
%     y=[y;T1std_1(nind2(sortind))-T1stdlr_1(nind2(iii))];
% end
% xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev.');
% % title('Prediction between Elements');
% saveas(gcf,'Figures/Pred_els_T1_lines','png');
figure; hold on;
[~,sortind]=sort(T1std_1(nind2));
for iii=1:length(nind2)
    plot(Hzmulr+T1Hlr_1(nind2(iii))-T1H1(nind2(sortind)),T1std_1(nind2(sortind))-T1stdlr_1(nind2(iii)),'bo');
end
xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev.');

x=[];y=[];
for iii=1:length(nind2)
    x=[x,Hzmulr+T1Hlr_1(nind2(iii))-T1H1(nind2(sortind))];
    y=[y;T1std_1(nind2(sortind))-T1stdlr_1(nind2(iii))];
end
[lfit,gof]=robustfit(x',y);
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-2.8,-0.3,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));

% title('Prediction between Elements');
saveas(gcf,'Figures/Pred_els_T1','png');

% figure; hold on;
% [~,sortind]=sort(T2std_1(nind2));
% for iii=1:length(nind2)
%     plot(Hzmulr+T2Hlr_1(nind2(iii))-T2H1(nind2(sortind(1:end-2))),T2std_1(nind2(sortind(1:end-2)))-T2stdlr_1(nind2(iii)),'-o');
% end
% xlabel('Conditional Mutual Information'); ylabel('Delta T2 Std. Dev.');
% % title('Prediction between Elements');
% saveas(gcf,'Figures/Pred_els_T2_lines','png');
figure; hold on;
[~,sortind]=sort(T2std_1(nind2));
for iii=1:length(nind2)
    plot(Hzmulr+T2Hlr_1(nind2(iii))-T2H1(nind2(sortind(1:end-2))),T2std_1(nind2(sortind(1:end-2)))-T2stdlr_1(nind2(iii)),'bo');
end
xlabel('Conditional Mutual Information'); ylabel('Delta T2 Std. Dev.');

x=[];y=[];
for iii=1:length(nind2)
    x=[x,Hzmulr+T2Hlr_1(nind2(iii))-T2H1(nind2(sortind))];
    y=[y;T2std_1(nind2(sortind))-T2stdlr_1(nind2(iii))];
end
[lfit,gof]=robustfit(x',y);
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-3.5,-0.2,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));

% title('Prediction between Elements');
saveas(gcf,'Figures/Pred_els_T2','png');


%% Intra Scan STD

x=Hzmulr+M0Hlr_1(nind2)-M0H2(nind2);
y=M0std_2(nind2)'-M0stdlr_1(nind2);
% [~,maxind]=max(y); x(maxind)=[]; y(maxind)=[];
% [~,minind]=min(y); x(minind)=[]; y(minind)=[];
figure; hold on;
plot(Hzmulr+M0Hlr_1(nind2)-M0H2(nind2),M0std_2(nind2)'-M0stdlr_1(nind2),'bo');
% h=lsline; set(h(1),'color','black');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-0.7,-0.03,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta M0 Std. Dev.');
% title('Prediction of Scan 1 High Res.');
saveas(gcf,'Figures/Pred_scan1_M0','png');

x=Hzmulr+M0Hlr_2(nind)-M0H1(nind);
y=M0std_1(nind)'-M0stdlr_2(nind);
figure; hold on;
plot(Hzmulr+M0Hlr_2(nind)-M0H1(nind),M0std_1(nind)'-M0stdlr_2(nind),'bo');
% h=lsline; set(h(1),'color','black');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-0.5,-0.1,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta M0 Std. Dev.');
% title('Prediction of Scan 2 High Res.');
saveas(gcf,'Figures/Pred_scan2_M0','png');

% figure; hold on;
% plot(Hzmulr+T1Hlr_1(nind)-T1H2(nind),T1std_2(nind)'-T1stdlr_1(nind),'bo');
% lsline;
% xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev. (s)');
% title('Prediction of Scan 1 High Res.');

figure; hold on;
y=T1std_2(nind([1:19,25:36]))'-T1stdlr_1(nind([1:19,25:36]));
x=Hzmulr+T1Hlr_1(nind([1:19,25:36]))-T1H2(nind([1:19,25:36]));
plot(Hzmulr+T1Hlr_1(nind([1:19,25:36]))-T1H2(nind([1:19,25:36])),T1std_2(nind([1:19,25:36]))'-T1stdlr_1(nind([1:19,25:36])),'bo');
% h=lsline; set(h(1),'color','black');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-0.4,-0.075,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev. (s)');
% title('Prediction of Scan 1 High Res.');
saveas(gcf,'Figures/Pred_scan1_T1','png');

x=Hzmulr+T1Hlr_2(nind([1:19,25:36]))-T1H1(nind([1:19,25:36]));
y=T1std_1(nind([1:19,25:36]))'-T1stdlr_2(nind([1:19,25:36]));
figure; hold on;
[~,sortind]=sort(Hzmulr+T1Hlr_2(nind)-T1H1(nind));
plot(Hzmulr+T1Hlr_2(nind([1:19,25:36]))-T1H1(nind([1:19,25:36])),T1std_1(nind([1:19,25:36]))'-T1stdlr_2(nind([1:19,25:36])),'bo');
% h=lsline; set(h(1),'color','black');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit);
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(-0.9,-0.06,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta T1 Std. Dev. (s)');
% title('Prediction of Scan 2 High Res.');
saveas(gcf,'Figures/Pred_scan2_T1','png');

x=Hzmulr+T2Hlr_1(nind(3:end))-T2H2(nind(3:end));
y=T2std_2(nind(3:end))'-T2stdlr_1(nind(3:end));
figure; hold on;
plot(Hzmulr+T2Hlr_1(nind(3:end))-T2H2(nind(3:end)),T2std_2(nind(3:end))'-T2stdlr_1(nind(3:end)),'bo');
% h=lsline; set(h(1),'color','black');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(0.60,0.08,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta T2 Std. Dev. (s)');
% title('Prediction of Scan 1 High Res.');
saveas(gcf,'Figures/Pred_scan1_T2','png');

x=Hzmulr+T2Hlr_2(nind(3:end))-T2H1(nind(3:end));
y=T2std_1(nind(3:end))'-T2stdlr_2(nind(3:end));
figure; hold on;
plot(Hzmulr+T2Hlr_2(nind(3:end))-T2H1(nind(3:end)),T2std_1(nind(3:end))'-T2stdlr_2(nind(3:end)),'bo');
% h=lsline; set(h(1),'color','black');
[lfit,gof]=robustfit(x',y');
% coef=coeffvalues(lfit)l;
xl=xlim;
plot(xl,lfit(2)*xl+lfit(1),'color','black');
text(0.60,.08,sprintf('y=%fx+%f\nr^2=%f\nRMSE=%f',lfit(2),lfit(1),gof.coeffcorr(2)^2,gof.s));
xlabel('Conditional Mutual Information'); ylabel('Delta T2 Std. Dev. (s)');
% title('Prediction of Scan 2 High Res.');
saveas(gcf,'Figures/Pred_scan2_T2','png');


%% Synthetic Patient Population Images
load('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/syntheticPatientPopulation.mat');
figure('pos',[10,10,1810,1210]);
for iii=1:10
    subplot(2,5,iii); imagesc(syntheticPatientPop{iii,1}(:,:,92)); axis off;
    title(sprintf('Model %i',iii));
end
saveas(gcf,'Figures/synthPatientPop_hires','png');
figure('pos',[10,10,1810,1210]);
for iii=1:10
    subplot(2,5,iii); imagesc(syntheticPatientPop{iii,2}(:,:,23)); axis off;
    title(sprintf('Model %i',iii));
end
saveas(gcf,'Figures/synthPatientPop_lores','png');

%% Phantom Images
scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
load([scanArchivePath,'.mat'],'realimg','bartrecon');

% for iii=1:5
%     tmpnii=load_untouch_nii([scanArchivePath,'_realimg.nii']);
%     tmpnii.img=realimg(:,:,:,iii);
%     save_untouch_nii(tmpnii,'tmpnii.nii.gz');
%     system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -interpolation NearestNeighbor -resample 50x50x50 -o resampleimg.nii.gz','tmpnii.nii.gz'));
%     tmprealimg = load_untouch_nii('resampleimg.nii.gz');
%     system('rm resampleimg.nii.gz');
%     system('rm tmpnii.nii.gz');
%     realimg_lores(:,:,:,iii) = tmprealimg.img;
% end


tmpim=fftshift(fftn(bartrecon(:,:,:,1)));
tmplo=zeros(size(tmpim));
tmplo(87:137,71:121,75:125)=tmpim(87:137,71:121,75:125);
tmplo=fftn(tmplo);

magnimg=abs(tmplo);
phaseimg=angle(bartrecon);
phasecorr=phaseimg-repmat(phaseimg(:,:,:,5),[1,1,1,5]);
realimg_lores=magnimg.*cos(phasecorr);
imagimg_lores=magnimg.*sin(phasecorr);


% tmpim=fftshift(fftn(realimg(:,:,:,1)));
% tmplo=zeros(size(tmpim));
% % kc=[size(tmplo)./2-25;size(tmplo)./2+25]';
% tmplo(87:137,71:121,75:125)=tmpim(87:137,71:121,75:125);
% tmplo=fftn(tmplo);

% figure;imagesc(realimg_lores(:,:,32,1)); axis off; colormap hot;
figure;imagesc(realimg(:,:,120,1)); axis off; colormap hot; caxis([-.5,2]);
saveas(gcf,'Figures/nistphantom_hires','png');

figure;imagesc(4.*realimg_lores(end:-1:1,end:-1:1,80)./max(realimg_lores(:))); axis off; colormap hot; caxis([-.5,2]);
saveas(gcf,'Figures/nistphantom_lores','png');

% scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
% load([scanArchivePath,'.mat'],'realimg');
