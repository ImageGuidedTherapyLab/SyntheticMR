
%% Make In Vivo Predictions
T1mean = [1200,  900, 900, 1200]./1000; % s
T1stdd = [ 100,  100,  100,  150]./1000; % s
T2mean = [ 100,   80, 80,  110]./1000; % s
T2stdd = [   5,    4,   4,   10]./1000; % s
M0mean = [ 0.9,  0.9,  0.9,  0.9];       % relative intensity
M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];
tisinputnovar=[M0mean;0.001*ones(size(M0stdd));T1mean;0.001*ones(size(T1stdd));T2mean;0.001*ones(size(T2stdd))];

flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.871;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
TD = [0.5,0.5,0.5,0.5];          % s
signu=3.4762E-4;
% signu=0.40;
acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD,signu];

lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
tmptissue = load_untouch_nii(lfname);
materialID=int32(tmptissue.img(15:165,20:200,92));

% Theoretical Optimal Case
popt=MI_QALAS_04242019_driverfun(tisinput,acqparam,materialID);
mipred_topt=MI_QALAS_objfun_04242019(popt,tisinput,acqparam,materialID);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(popt,tisinput,acqparam,1,500,1);
topt_predstd(1:3)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(popt,tisinput,acqparam,2,500,1);
topt_predstd(4:6)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(popt,tisinput,acqparam,3,500,1);
topt_predstd(7:9)=nanstd(post_samples(:,1:3),1);

[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(popt,tisinputnovar,acqparam,1,500,1);
topt_tminstd(1:3)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(popt,tisinputnovar,acqparam,2,500,1);
topt_tminstd(4:6)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(popt,tisinputnovar,acqparam,3,500,1);
topt_tminstd(7:9)=nanstd(post_samples(:,1:3),1);

% Scan 1
td_scan1=[0.0234,0.117448,0.117448,0.117448,0];
mipred_scan1=MI_QALAS_objfun_04242019(td_scan1,tisinput,acqparam,materialID);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan1,tisinput,acqparam,1,500,1);
scan1_predstd(1:3)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan1,tisinput,acqparam,2,500,1);
scan1_predstd(4:6)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan1,tisinput,acqparam,3,500,1);
scan1_predstd(7:9)=nanstd(post_samples(:,1:3),1);

[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan1,tisinputnovar,acqparam,1,500,1);
scan1_tminstd(1:3)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan1,tisinputnovar,acqparam,2,500,1);
scan1_tminstd(4:6)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan1,tisinputnovar,acqparam,3,500,1);
scan1_tminstd(7:9)=nanstd(post_samples(:,1:3),1);

% Scan 2
td_scan2=[0.500,0.150,0.0005,0.0005,0.500];
mipred_scan2=MI_QALAS_objfun_04242019(td_scan2,tisinput,acqparam,materialID);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan2,tisinput,acqparam,1,500,1);
scan2_predstd(1:3)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan2,tisinput,acqparam,2,500,1);
scan2_predstd(4:6)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan2,tisinput,acqparam,3,500,1);
scan2_predstd(7:9)=nanstd(post_samples(:,1:3),1);

[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan2,tisinputnovar,acqparam,1,500,1);
scan2_tminstd(1:3)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan2,tisinputnovar,acqparam,2,500,1);
scan2_tminstd(4:6)=nanstd(post_samples(:,1:3),1);
[logZ,nest_samples,post_samples,prior]=mn_qalasnp_04242019(td_scan2,tisinputnovar,acqparam,3,500,1);
scan2_tminstd(7:9)=nanstd(post_samples(:,1:3),1);


%% Process In Vivo Data

scanArchivePath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/ScanArchive_713792CMR1_20190424_160439142';
load([scanArchivePath,'_parampred.mat']);
M0pred1=M0pred; T1pred1=T1pred; T2pred1=real(T2pred);
realimg1_scan1=load_untouch_nii([scanArchivePath,'_realimg1.nii']);
realimg1_scan1=realimg1_scan1.img;
realimg5_scan1=load_untouch_nii([scanArchivePath,'_realimg.nii']);
realimg5_scan1=realimg5_scan1.img;

scanArchivePath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/ScanArchive_713792CMR1_20190424_161933252';
load([scanArchivePath,'_parampred.mat']);
M0pred2=M0pred; T1pred2=T1pred; T2pred2=real(T2pred);

clear M0pred T1pred T2pred;

segimgPath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/invivo_segimg.nii.gz';
segimg=load_untouch_nii(segimgPath);
segimg=segimg.img;

for iii=1:max(segimg(:))
    m0_roi_scan1{iii}=M0pred1(segimg==iii);
    m0_roi_scan2{iii}=M0pred2(segimg==iii);
    t1_roi_scan1{iii}=T1pred1(segimg==iii);
    t1_roi_scan2{iii}=T1pred2(segimg==iii);
    t2_roi_scan1{iii}=T2pred1(segimg==iii);
    t2_roi_scan2{iii}=T2pred2(segimg==iii);
    
    m0std_scan1(iii)=nanstd(m0_roi_scan1{iii});
    m0std_scan2(iii)=nanstd(m0_roi_scan2{iii});
    m0relstd_scan1(iii)=nanstd(m0_roi_scan1{iii})/abs(nanmean(m0_roi_scan1{iii}));
    m0relstd_scan2(iii)=nanstd(m0_roi_scan2{iii})/abs(nanmean(m0_roi_scan2{iii}));
    [h_m0(iii),p_m0(iii)]=vartest2(m0_roi_scan1{iii},m0_roi_scan2{iii});
    
    t1std_scan1(iii)=nanstd(t1_roi_scan1{iii});
    t1std_scan2(iii)=nanstd(t1_roi_scan2{iii});
    t1relstd_scan1(iii)=nanstd(t1_roi_scan1{iii})/abs(nanmean(t1_roi_scan1{iii}));
    t1relstd_scan2(iii)=nanstd(t1_roi_scan2{iii})/abs(nanmean(t1_roi_scan2{iii}));
    [h_t1(iii),p_t1(iii)]=vartest2(t1_roi_scan1{iii},t1_roi_scan2{iii});
    
    t2std_scan1(iii)=nanstd(t2_roi_scan1{iii});
    t2std_scan2(iii)=nanstd(t2_roi_scan2{iii});
    t2relstd_scan1(iii)=nanstd(t2_roi_scan1{iii})/abs(nanmean(t2_roi_scan1{iii}));
    t2relstd_scan2(iii)=nanstd(t2_roi_scan2{iii})/abs(nanmean(t2_roi_scan2{iii}));
    [h_t2(iii),p_t2(iii)]=vartest2(t2_roi_scan1{iii},t2_roi_scan2{iii});
end

save /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/invivo_results.mat -v7.3;


%% Make Figures

load('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/invivo_results.mat');

figure('Pos',[10,10,810,810]); bar([scan1_predstd',scan2_predstd',topt_predstd',topt_tminstd']);
xlabel('Tissue and Parametric Map');ylabel('Predicted Standard Deviation (s for relaxation times)');
xticklabels({'GM M0','GM T1','GM T2','WM M0','WM T1','WM T2','CSF M0','CSF T1','CSF T2'});
legend('Clinical Default','Clinical Optimized','Theoretical Optimal','Theoretical Minimum');
saveas(gcf,'Figures/invivo_predstd','png');

% figure; bar([m0std_scan1',m0std_scan2']);
% figure; bar(100*[m0relstd_scan1([1:12,25:28])',m0relstd_scan2([1:12,25:28])']);
% figure; bar([t1std_scan1',t1std_scan2']);
% figure; bar([t1relstd_scan1',t1relstd_scan2']);
% figure; bar([t2std_scan1',t2std_scan2']);
% figure; bar([t2relstd_scan1',t2relstd_scan2']);
% 
% figure; bar(t1std_scan1'-t1std_scan2');
% 
% figure; bar(100*[m0relstd_scan1'-m0relstd_scan2']);
% figure; bar(100*[t1relstd_scan1'-t1relstd_scan2']);
% figure; bar(100*[t2relstd_scan1'-t2relstd_scan2']);
% 
% figure; boxplot([t1_roi_scan1{1},t1_roi_scan2{1}]);


skstr130=load_untouch_nii('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/skullstrip130.nii.gz');
slicenum=130;

figure('pos',[10,10,1810,610]);
subplot(1,3,1); imagesc(M0pred1(30:235,25:160,slicenum).*~(skstr130.img(30:235,25:160)==1));
colormap hot; axis off; title('M0');
subplot(1,3,2); imagesc(T1pred1(30:235,25:160,slicenum).*~(skstr130.img(30:235,25:160)==1));
colormap hot; axis off; caxis([0,4.5]); title('T1');
subplot(1,3,3); imagesc(T2pred1(30:235,25:160,slicenum).*~(skstr130.img(30:235,25:160)==1));
colormap hot; axis off; caxis([0,.5]); title('T2');
saveas(gcf,'Figures/invivo_scan1_slice130','png');

figure('pos',[10,10,1810,610]);
subplot(1,3,1); imagesc(M0pred2(30:235,25:160,slicenum).*~(skstr130.img(30:235,25:160)==1));
colormap hot; axis off; title('M0');
subplot(1,3,2); imagesc(T1pred2(30:235,25:160,slicenum).*~(skstr130.img(30:235,25:160)==1));
colormap hot; axis off; caxis([0,4.5]); title('T1');
subplot(1,3,3); imagesc(T2pred2(30:235,25:160,slicenum).*~(skstr130.img(30:235,25:160)==1));
colormap hot; axis off; caxis([0,.5]); title('T2');
saveas(gcf,'Figures/invivo_scan2_slice130','png');

testind=[2,6,9,12,13,18,22,25,28];
m0ind=[1,2,4,6,8,11,18,25,26];
% m0ind=[2,6,8,11,15,18,23,24,25];
t1ind=[2,6,8,10,15,18,23,24,25];
t2ind=[2,4,7,9,13,14,19,20,21];
% t2ind=[3,4,7,9,13,15,19,21,22];
testind2=[3,4,7,9,13,16,19,21,22];
figure('pos',[10,10,810,810]);
bar(-[m0relstd_scan1(m0ind)'-m0relstd_scan2(m0ind)',t1relstd_scan1(t1ind)'-t1relstd_scan2(t1ind)',t2relstd_scan1(t2ind)'-t2relstd_scan2(t2ind)']);
% bar(-[m0relstd_scan1(testind)'-m0relstd_scan2(testind)',t1relstd_scan1(testind)'-t1relstd_scan2(testind)',t2relstd_scan1(testind2)'-t2relstd_scan2(testind2)']);
xlabel('ROI Location'); ylabel('Change in Std. Dev. after Optimization');
xticklabels({'Fr. WM','Fr. GM','Pa. WM','Pa. GM','Tm. WM','Tm. GM','Oc. WM','Oc. GM','CSF'});
legend('M0','T1','T2','location','northeast');
saveas(gcf,'Figures/invivo_deltastd','png');


%% Example Images

for iii=1:28
[h_m0(iii),p_m0(iii)]=vartest2(m0_roi_scan1{iii},m0_roi_scan2{iii}.*median(m0_roi_scan1{iii})./median(m0_roi_scan2{iii}));
[h_t1(iii),p_t1(iii)]=vartest2(t1_roi_scan1{iii},t1_roi_scan2{iii}.*median(t1_roi_scan1{iii})./median(t1_roi_scan2{iii}));
[h_t2(iii),p_t2(iii)]=vartest2(t2_roi_scan1{iii},t2_roi_scan2{iii}.*median(t2_roi_scan1{iii})./median(t2_roi_scan2{iii}));
end

% 134 135
m12=mean(t1_roi_scan1{2});
m16=mean(t1_roi_scan1{6});
m22=mean(t1_roi_scan2{2});
m26=mean(t1_roi_scan2{6});
s12=std(t1_roi_scan1{2});
s16=std(t1_roi_scan1{6});
s22=std(t1_roi_scan2{2});
s26=std(t1_roi_scan2{6});
figure;
subplot(1,2,1); imagesc(T1pred1(30:235,25:160,135)); colormap hot; caxis([m12-2*s12,m12+2*s12]);
subplot(1,2,2); imagesc(T1pred2(30:235,25:160,135).*m12./m22); colormap hot; caxis([m12-2*s12,m12+2*s12]);

% tmp=T1pred1(30:235,25:160,145);
% tmpnii=make_nii(tmp);
% save_nii(tmpnii,'tmpnii.nii');
% system('vglrun itksnap -g tmpnii.nii');

tmpseg=load_untouch_nii('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/invivoexampleseg.nii.gz');
tmpseg=tmpseg.img;

figure;
subplot(1,2,1); imagesc(T1pred1(30:235,25:160,145).*~tmpseg); colormap gray; caxis([m12-1.5*s12,m12+1.5*s12]); axis off;
subplot(1,2,2); imagesc(~tmpseg.*T1pred2(30:235,25:160,147).*m12./m22); colormap gray; caxis([m12-1.5*s12,m12+1.5*s12]); axis off;
saveas(gcf,'Figures/invivo_examplediff','png');


figure;
subplot(1,2,1); imagesc(T1pred1(30:235,25:160,134)); colormap hot; caxis([m16-2*s16,m16+2*s16]);
subplot(1,2,2); imagesc(T1pred2(30:235,25:160,134)*m16/m26); colormap hot; caxis([m16-2*s16,m16+2*s16]);

figure;
subplot(1,2,1); imagesc(T1pred1(30:235,25:160,145)); colormap hot; caxis([m16-1.5*s16,m16+1.5*s16]);
subplot(1,2,2); imagesc(T1pred2(30:235,25:160,147)*m16/m26); colormap hot; caxis([m16-1.5*s16,m16+1.5*s16]);