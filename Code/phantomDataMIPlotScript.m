
gibbs_exclude=0;
MIcalcflag=2;

scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
fileflag=1;
load([scanArchivePath,'_parampred.mat']);

if gibbs_exclude~=0
    Mv(~gibbs_mask)=nan;
    M0v(~gibbs_mask)=nan;
    T1v(~gibbs_mask)=nan;
    T2v(~gibbs_mask)=nan;
    M5v(~gibbs_mask)=nan;
end
M0v_1=M0v; T1v_1=T1v; T2v_1=T2v; %signuin1=nanvar(M5v(15:end,:),[],2)';
scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
fileflag=2;
load([scanArchivePath,'_parampred.mat']);
if gibbs_exclude~=0
    Mv(~gibbs_mask)=nan;
    M0v(~gibbs_mask)=nan;
    T1v(~gibbs_mask)=nan;
    T2v(~gibbs_mask)=nan;
    M5v(~gibbs_mask)=nan;
end
M0v_2=M0v; T1v_2=T1v; T2v_2=T2v; %signuin2=nanvar(M5v(15:end,:),[],2)';

load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
M0=nanmedian(M0v(15:end,:),2)';
T1=phantomProp.T1element.T30.T1./1000; T2=phantomProp.T1element.T30.T2./1000;
T1=[T1,phantomProp.T2element.T30.T1./1000]; T2=[T2,phantomProp.T2element.T30.T2./1000];

signu=3.4762E-4;

switch MIcalcflag
    case 1
        for iii=1:28
            %     [MItest1(iii,:)]=QALAS_synphan_MIcalc1p_objfun(1,[],[M0(iii),T1(iii),T2(iii),[M0(iii),T1(iii),T2(iii)].*.05,signu]);
            [MItest1(iii,:)]=QALAS_synphan_MIcalc1p1t_objfun(1,[],[M0(iii),T1(iii),T2(iii),[M0(iii),T1(iii),T2(iii)].*.05,signu],2);
            [MItest2(iii,:)]=QALAS_synphan_MIcalc1p1t_objfun(2,[],[M0(iii),T1(iii),T2(iii),[M0(iii),T1(iii),T2(iii)].*.05,signu],2);
        end
        MItest1=[zeros(14,5);MItest1];
        MItest2=[zeros(14,5);MItest2];
        
        figure; hold on;
        plot(nanvar(M0v_1(3:14,:),[],2)./nanmedian(M0v_1(3:14,:),2),'o');
        plot(nanvar(M0v_2(3:14,:),[],2)./nanmedian(M0v_2(3:14,:),2),'x');
        xlabel('M0 Element Number'); ylabel('M0 Relative Variance');
        title('M0 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(nanvar(M0v_1(17:28,:),[],2)./nanmedian(M0v_1(17:28,:),2),'o');
        plot(nanvar(M0v_2(17:28,:),[],2)./nanmedian(M0v_2(17:28,:),2),'x');
        xlabel('M0 Element Number'); ylabel('M0 Relative Variance');
        title('T1 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(nanvar(M0v_1(31:42,:),[],2)./nanmedian(M0v_1(31:42,:),2),'o');
        plot(nanvar(M0v_2(31:42,:),[],2)./nanmedian(M0v_2(31:42,:),2),'x');
        xlabel('T2 Element Number'); ylabel('M0 Relative Variance');
        title('T2 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(nanvar(T1v_1(17:28,:),[],2)./nanmedian(T1v_1(17:28,:),2),'o');
        plot(nanvar(T1v_2(17:28,:),[],2)./nanmedian(T1v_2(17:28,:),2),'x');
        xlabel('T1 Element Number'); ylabel('T1 Relative Variance');
        title('T1 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(nanvar(T1v_1(31:42,:),[],2)./nanmedian(T1v_1(31:42,:),2),'o');
        plot(nanvar(T1v_2(31:42,:),[],2)./nanmedian(T1v_2(31:42,:),2),'x');
        xlabel('T2 Element Number'); ylabel('T1 Relative Variance');
        title('T2 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(nanvar(T2v_1(17:28,:),[],2)./nanmedian(T2v_1(17:28,:),2),'o');
        plot(nanvar(T2v_2(17:28,:),[],2)./nanmedian(T2v_2(17:28,:),2),'x');
        xlabel('T2 Element Number'); ylabel('T2 Relative Variance');
        title('T1 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(nanvar(T2v_1(31:42,:),[],2)./nanmedian(T2v_1(31:42,:),2),'o');
        plot(nanvar(T2v_2(31:42,:),[],2)./nanmedian(T2v_2(31:42,:),2),'x');
        xlabel('T2 Element Number'); ylabel('T2 Relative Variance');
        title('T2 Elements');
        legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
        
        figure; hold on;
        plot(sum(MItest1(17:28,:),2),nanvar(M0v_1(17:28,:),[],2)./nanmedian(M0v_1(17:28,:),2),'o'); lsline;
        mdl1=fitlm(sum(MItest1(17:28,:),2),nanvar(M0v_1(17:28,:),[],2)./nanmedian(M0v_1(17:28,:),2));
        plot(sum(MItest1(31:42,:),2),nanvar(M0v_1(31:42,:),[],2)./nanmedian(M0v_1(31:42,:),2),'x'); lsline;
        mdl2=fitlm(sum(MItest1(31:42,:),2),nanvar(M0v_1(31:42,:),[],2)./nanmedian(M0v_1(31:42,:),2));
        plot(sum(MItest2(17:28,:),2),nanvar(M0v_2(17:28,:),[],2)./nanmedian(M0v_2(17:28,:),2),'*'); lsline;
        mdl3=fitlm(sum(MItest2(17:28,:),2),nanvar(M0v_2(17:28,:),[],2)./nanmedian(M0v_2(17:28,:),2));
        plot(sum(MItest2(31:42,:),2),nanvar(M0v_2(31:42,:),[],2)./nanmedian(M0v_2(31:42,:),2),'sq'); lsline;
        mdl4=fitlm(sum(MItest2(31:42,:),2),nanvar(M0v_2(31:42,:),[],2)./nanmedian(M0v_2(31:42,:),2));
        xlabel('Mutual Information'); ylabel('M0 Relative Variance');
        legend('T1 Elements, Scan 1',sprintf('R^2=%2.2f',mdl1.Rsquared.Ordinary),'T2 Elements, Scan 1',sprintf('R^2=%2.2f',mdl2.Rsquared.Ordinary),...
            'T1 Elements, Scan 2',sprintf('R^2=%2.2f',mdl3.Rsquared.Ordinary),'T2 Elements, Scan 2',sprintf('R^2=%2.2f',mdl4.Rsquared.Ordinary),...
            'Location','southwest');
        
        % figure; hold on;
        % plot(sum(MItest1(17:28,:),2),nanvar(T1v_1(17:28,:)./nanmean(T1v_1(17:28,:)),[],2),'o');
        % plot(sum(MItest2(17:28,:),2),nanvar(T1v_2(17:28,:)./nanmean(T1v_2(17:28,:)),[],2),'o');
        % plot(sum(MItest1(31:42,:),2),nanvar(T1v_1(31:42,:)./nanmean(T1v_1(31:42,:)),[],2),'o');
        % plot(sum(MItest2(31:42,:),2),nanvar(T1v_2(31:42,:)./nanmean(T1v_2(31:42,:)),[],2),'o');
        
        figure; hold on;
        plot(sum(MItest1(17:28,:),2),nanvar(T1v_1(17:28,:)./nanmean(T1v_1(17:28,:),2),[],2),'o'); lsline;
        mdl1=fitlm(sum(MItest1(17:28,:),2),nanvar(T1v_1(17:28,:)./nanmean(T1v_1(17:28,:),2),[],2));
        plot(sum(MItest1(31:42,:),2),nanvar(T1v_1(31:42,:)./nanmean(T1v_1(31:42,:),2),[],2),'x'); lsline;
        mdl2=fitlm(sum(MItest1(31:42,:),2),nanvar(T1v_1(31:42,:)./nanmean(T1v_1(31:42,:),2),[],2));
        plot(sum(MItest2(17:28,:),2),nanvar(T1v_2(17:28,:)./nanmean(T1v_2(17:28,:),2),[],2),'*'); lsline;
        mdl3=fitlm(sum(MItest2(17:28,:),2),nanvar(T1v_2(17:28,:)./nanmean(T1v_2(17:28,:),2),[],2));
        plot(sum(MItest2(31:42,:),2),nanvar(T1v_2(31:42,:)./nanmean(T1v_2(31:42,:),2),[],2),'sq'); lsline;
        mdl4=fitlm(sum(MItest2(31:42,:),2),nanvar(T1v_2(31:42,:)./nanmean(T1v_2(31:42,:),2),[],2));
        xlabel('Mutual Information'); ylabel('T1 Relative Variance');
        legend('T1 Elements, Scan 1',sprintf('R^2=%2.2f',mdl1.Rsquared.Ordinary),'T2 Elements, Scan 1',sprintf('R^2=%2.2f',mdl2.Rsquared.Ordinary),...
            'T1 Elements, Scan 2',sprintf('R^2=%2.2f',mdl3.Rsquared.Ordinary),'T2 Elements, Scan 2',sprintf('R^2=%2.2f',mdl4.Rsquared.Ordinary),...
            'Location','southwest');
        
        figure; hold on;
        plot(sum(MItest1(17:28,:),2),nanvar(T2v_1(17:28,:)./nanmedian(T2v_1(17:28,:),2),[],2),'o'); %lsline;
        plot(sum(MItest1(31:42,:),2),nanvar(T2v_1(31:42,:)./nanmedian(T2v_1(31:42,:),2),[],2),'x'); %lsline;
        plot(sum(MItest2(17:28,:),2),nanvar(T2v_2(17:28,:)./nanmedian(T2v_2(17:28,:),2),[],2),'*'); %lsline;
        plot(sum(MItest2(31:42,:),2),nanvar(T2v_2(31:42,:)./nanmedian(T2v_2(31:42,:),2),[],2),'sq'); %lsline;
        xlabel('Mutual Information'); ylabel('T2 Relative Variance');
        legend('T1 Elements, Scan 1','T2 Elements, Scan 1',...
            'T1 Elements, Scan 2','T2 Elements, Scan 2',...
            'Location','northwest');
        
    case 2
        for iii=1:42
            M0MI1(iii)=0.5.*log((2*pi*2.7183)^2.*(sqrt(M0(1).*.05)).^4)-entropy(M0v_1(iii,:)-mean(M0v_1(iii,:))+.5);
            M0MI2(iii)=0.5.*log((2*pi*2.7183)^2.*(sqrt(M0(1).*.05)).^4)-entropy(M0v_2(iii,:)-mean(M0v_2(iii,:))+.5);
            T1MI1(iii)=0.5.*log((2*pi*2.7183)^2.*(sqrt(T1(3).*.05)).^4)-entropy(T1v_1(iii,:)-mean(T1v_1(iii,:))+.5);
            T1MI2(iii)=0.5.*log((2*pi*2.7183)^2.*(sqrt(T1(3).*.05)).^4)-entropy(T1v_2(iii,:)-mean(T1v_2(iii,:))+.5);
            T2MI1(iii)=0.5.*log((2*pi*2.7183)^2.*(sqrt(T2(17).*.05)).^4)-entropy(real(T2v_1(iii,:))-mean(real(T2v_1(iii,:)))+.5);
            T2MI2(iii)=0.5.*log((2*pi*2.7183)^2.*(sqrt(T2(17).*.05)).^4)-entropy(real(T2v_2(iii,:))-mean(real(T2v_2(iii,:)))+.5);
        end
%         MItest1=[zeros(1,14),MItest1];
%         MItest2=[zeros(1,14),MItest2];

figure; 
subplot(1,3,1); hold on;
plot(M0MI1(3:14),nanstd(M0v_1(3:14,:),[],2),'o');
% yyaxis right;
plot(M0MI2(3:14),nanstd(M0v_2(3:14,:),[],2),'x');
xlabel('Mutual Information'); ylabel('M0 Standard Deviation');
legend('Scan 1','Scan 2');

subplot(1,3,2); hold on;
plot(T1MI1(17:28),nanstd(T1v_1(17:28,:),[],2),'o');
xlabel('Mutual Information'); ylabel('T1 Standard Deviation');
% yyaxis right;
plot(T1MI2(17:28),nanstd(T1v_2(17:28,:),[],2),'x');
xlabel('Mutual Information'); ylabel('T1 Standard Deviation');
legend('Scan 1','Scan 2');

subplot(1,3,3); hold on;
plot(T2MI1(31:42),nanstd(T2v_1(31:42,:),[],2),'o');
xlabel('Mutual Information'); ylabel('T2 Standard Deviation');
% yyaxis right;
plot(T2MI2(31:42),nanstd(T2v_2(31:42,:),[],2),'x');
xlabel('Mutual Information'); ylabel('T2 Standard Deviation');
legend('Scan 1','Scan 2');
        
end

