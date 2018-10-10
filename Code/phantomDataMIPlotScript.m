
scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
fileflag=1;
load([scanArchivePath,'_parampred.mat']);
M0v_1=M0v; T1v_1=T1v; T2v_1=T2v;
scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
fileflag=2;
load([scanArchivePath,'_parampred.mat']);
M0v_2=M0v; T1v_2=T1v; T2v_2=T2v;

load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/phantomProp.mat;
M0=median(M0v(15:end,:),2)';
T1=phantomProp.T1element.T30.T1./1000; T2=phantomProp.T1element.T30.T2./1000;
T1=[T1,phantomProp.T2element.T30.T1./1000]; T2=[T2,phantomProp.T2element.T30.T2./1000];

signuin=var(M5v(15:end,:),[],2)'; 

for iii=1:28
[MItest1(iii,:)]=QALAS_synphan_MIcalc1p_objfun(1,[],[M0(iii),T1(iii),T2(iii),[M0(iii),T1(iii),T2(iii)].*.05,signuin(iii)]);
 end
 
  for iii=1:28
[MItest2(iii,:)]=QALAS_synphan_MIcalc1p_objfun(2,[],[M0(iii),T1(iii),T2(iii),[M0(iii),T1(iii),T2(iii)].*.05,signuin(iii)]);
  end
  
  MItest1=[zeros(14,5);MItest1];
  MItest2=[zeros(14,5);MItest2];
  
figure; hold on;
plot(var(M0v_1(3:14,:),[],2)./median(M0v_1(3:14,:),2),'o');
plot(var(M0v_2(3:14,:),[],2)./median(M0v_2(3:14,:),2),'x');
xlabel('M0 Element Number'); ylabel('M0 Relative Variance');
title('M0 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');

figure; hold on;
plot(var(M0v_1(17:28,:),[],2)./median(M0v_1(17:28,:),2),'o');
plot(var(M0v_2(17:28,:),[],2)./median(M0v_2(17:28,:),2),'x');
xlabel('M0 Element Number'); ylabel('M0 Relative Variance');
title('T1 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');

figure; hold on;
plot(var(M0v_1(31:42,:),[],2)./median(M0v_1(31:42,:),2),'o');
plot(var(M0v_2(31:42,:),[],2)./median(M0v_2(31:42,:),2),'x');
xlabel('T2 Element Number'); ylabel('M0 Relative Variance');
title('T2 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');
  
figure; hold on;
plot(var(T1v_1(17:28,:),[],2)./median(T1v_1(17:28,:),2),'o');
plot(var(T1v_2(17:28,:),[],2)./median(T1v_2(17:28,:),2),'x');
xlabel('T1 Element Number'); ylabel('T1 Relative Variance');
title('T1 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');

figure; hold on;
plot(var(T1v_1(31:42,:),[],2)./median(T1v_1(31:42,:),2),'o');
plot(var(T1v_2(31:42,:),[],2)./median(T1v_2(31:42,:),2),'x');
xlabel('T2 Element Number'); ylabel('T1 Relative Variance');
title('T2 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');

figure; hold on;
plot(var(T2v_1(17:28,:),[],2)./median(T2v_1(17:28,:),2),'o');
plot(var(T2v_2(17:28,:),[],2)./median(T2v_2(17:28,:),2),'x');
xlabel('T2 Element Number'); ylabel('T2 Relative Variance');
title('T1 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');

figure; hold on;
plot(var(T2v_1(31:42,:),[],2)./median(T2v_1(31:42,:),2),'o');
plot(var(T2v_2(31:42,:),[],2)./median(T2v_2(31:42,:),2),'x');
xlabel('T2 Element Number'); ylabel('T2 Relative Variance');
title('T2 Elements');
legend('Scan 1 (MI=7.29E5)','Scan 2 (MI=7.49E5)');

figure; hold on;
plot(sum(MItest1(17:28,:),2),var(M0v_1(17:28,:),[],2)./median(M0v_1(17:28,:),2),'o'); lsline;
mdl1=fitlm(sum(MItest1(17:28,:),2),var(M0v_1(17:28,:),[],2)./median(M0v_1(17:28,:),2));
plot(sum(MItest1(31:42,:),2),var(M0v_1(31:42,:),[],2)./median(M0v_1(31:42,:),2),'x'); lsline;
mdl2=fitlm(sum(MItest1(31:42,:),2),var(M0v_1(31:42,:),[],2)./median(M0v_1(31:42,:),2));
plot(sum(MItest2(17:28,:),2),var(M0v_2(17:28,:),[],2)./median(M0v_2(17:28,:),2),'*'); lsline;
mdl3=fitlm(sum(MItest2(17:28,:),2),var(M0v_2(17:28,:),[],2)./median(M0v_2(17:28,:),2));
plot(sum(MItest2(31:42,:),2),var(M0v_2(31:42,:),[],2)./median(M0v_2(31:42,:),2),'sq'); lsline;
mdl4=fitlm(sum(MItest2(31:42,:),2),var(M0v_2(31:42,:),[],2)./median(M0v_2(31:42,:),2));
xlabel('Mutual Information'); ylabel('M0 Relative Variance');
legend('T1 Elements, Scan 1',sprintf('R^2=%2.2f',mdl1.Rsquared.Ordinary),'T2 Elements, Scan 1',sprintf('R^2=%2.2f',mdl2.Rsquared.Ordinary),...
    'T1 Elements, Scan 2',sprintf('R^2=%2.2f',mdl3.Rsquared.Ordinary),'T2 Elements, Scan 2',sprintf('R^2=%2.2f',mdl4.Rsquared.Ordinary),...
    'Location','southwest');

figure; hold on;
plot(sum(MItest1(17:28,:),2),var(T1v_1(17:28,:),[],2)./median(T1v_1(17:28,:),2),'o'); lsline;
mdl1=fitlm(sum(MItest1(17:28,:),2),var(T1v_1(17:28,:),[],2)./median(T1v_1(17:28,:),2));
plot(sum(MItest1(31:42,:),2),var(T1v_1(31:42,:),[],2)./median(T1v_1(31:42,:),2),'x'); lsline;
mdl2=fitlm(sum(MItest1(31:42,:),2),var(T1v_1(31:42,:),[],2)./median(T1v_1(31:42,:),2));
plot(sum(MItest2(17:28,:),2),var(T1v_2(17:28,:),[],2)./median(T1v_2(17:28,:),2),'*'); lsline;
mdl3=fitlm(sum(MItest2(17:28,:),2),var(T1v_2(17:28,:),[],2)./median(T1v_2(17:28,:),2));
plot(sum(MItest2(31:42,:),2),var(T1v_2(31:42,:),[],2)./median(T1v_2(31:42,:),2),'sq'); lsline;
mdl4=fitlm(sum(MItest2(31:42,:),2),var(T1v_2(31:42,:),[],2)./median(T1v_2(31:42,:),2));
xlabel('Mutual Information'); ylabel('T1 Relative Variance');
legend('T1 Elements, Scan 1',sprintf('R^2=%2.2f',mdl1.Rsquared.Ordinary),'T2 Elements, Scan 1',sprintf('R^2=%2.2f',mdl2.Rsquared.Ordinary),...
    'T1 Elements, Scan 2',sprintf('R^2=%2.2f',mdl3.Rsquared.Ordinary),'T2 Elements, Scan 2',sprintf('R^2=%2.2f',mdl4.Rsquared.Ordinary),...
    'Location','southwest');

figure; hold on;
plot(sum(MItest1(17:28,:),2),var(T2v_1(17:28,:),[],2)./median(T2v_1(17:28,:),2),'o');
plot(sum(MItest1(31:42,:),2),var(T2v_1(31:42,:),[],2)./median(T2v_1(31:42,:),2),'x');
plot(sum(MItest2(17:28,:),2),var(T2v_2(17:28,:),[],2)./median(T2v_2(17:28,:),2),'*');
plot(sum(MItest2(31:42,:),2),var(T2v_2(31:42,:),[],2)./median(T2v_2(31:42,:),2),'sq');