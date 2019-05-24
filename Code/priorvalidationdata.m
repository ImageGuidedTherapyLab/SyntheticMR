
%% Generate in vivo M0/T1/T2 correlation data

symriflag=0;
if symriflag==1

    load '/rsrch1/ip'/dmitchell2/github/SyntheticMR/Code/priorvalidationdata_pathnames.mat;
    
    m0=dicomread(filename1);
    t1=dicomread(filename2);
    t2=dicomread(filename3);
    
    segimg=load_untouch_nii([pathname,'tmpseg.nii.gz']);
    segimg=segimg.img;
    
    z_gm=m0(segimg==1);
    z_gm=[z_gm,t1(segimg==1)];
    z_gm=[z_gm,t2(segimg==1)];
    z_wm=m0(segimg==2);
    z_wm=[z_wm,t1(segimg==2)];
    z_wm=[z_wm,t2(segimg==2)];
    z_csf=m0(segimg==3);
    z_csf=[z_csf,t1(segimg==3)];
    z_csf=[z_csf,t2(segimg==3)];
    z=[z_gm(1:452,:),z_wm(1:452,:),z_csf];
    
else
    load('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/invivo_results.mat');

    % Handle outliers
    [~,ind]=min(t2_roi_scan1{28});
    t2_roi_scan1{28}(ind)=median(t2_roi_scan1{28});
    [~,ind]=max(t2_roi_scan1{6});
    t2_roi_scan1{6}(ind)=median(t2_roi_scan1{6});
    
    z=[m0_roi_scan1{6},t1_roi_scan1{6},t2_roi_scan1{6},...
        m0_roi_scan1{2},t1_roi_scan1{2},t2_roi_scan1{2},...
        m0_roi_scan1{28},t1_roi_scan1{28},t2_roi_scan1{28}];

end

zcov=cov(double(z));
zcorr=corr(double(z));

%% Plot in vivo M0/T1/T2 correlations

plotcorr(z,[0,13,0,4,0,1.5,0,13,0,4,0,1.5,0,13,0,4,0,1.5],[0,13,0,3.99,0,1.5,0,13,0,3.99,0,1.5,0,13,0,3.99,0,1.5]);
saveas(gcf,'Figures/priorval_zmeas','png');


%% MultiNest Prior Testing

T1mean = [1200,  900, 4000, 1200]./1000; % s
T1stdd = [ 100,  100,  200,  150]./1000; % s
T2mean = [ 100,   80, 1000,  110]./1000; % s
T2stdd = [   5,    4,   50,   10]./1000; % s
M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];

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

Nlive=500;
Nmcmc=5;
TDin=[.5,.5,.5,.5,.5];

[logZ,nest_samples,post_samples_unif_gm,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,1,Nlive,Nmcmc,'uniform');
[logZ,nest_samples,post_samples_unif_wm,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,2,Nlive,Nmcmc,'uniform');
[logZ,nest_samples,post_samples_unif_csf,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,3,Nlive,Nmcmc,'uniform');
psmin=min([size(post_samples_unif_gm,1),size(post_samples_unif_wm,1),size(post_samples_unif_csf,1)]);
z_unif=[post_samples_unif_gm(1:psmin,1:3),post_samples_unif_wm(1:psmin,1:3),post_samples_unif_csf(1:psmin,1:3)];

md=[tisinput(1,1:3),tisinput(3,1:3),tisinput(5,1:3)];
stdd=[tisinput(2,1:3),tisinput(4,1:3),tisinput(6,1:3)];
plotvec=zeros([1,2*length(md)]);
plotvec(1:2:end)=mean(z_unif,1)-3*stdd; plotvec(2:2:end)=mean(z_unif,1)+3*stdd;
plotcorr(z_unif,plotvec,[]);
saveas(gcf,'Figures/priorval_zunif','png');

mu_unif=[M0mean(1)-2*M0stdd(1) + 4*M0stdd(1).*rand(psmin,1),...
        M0mean(2)-2*M0stdd(2) + 4*M0stdd(2).*rand(psmin,1),...
        M0mean(3)-2*M0stdd(3) + 4*M0stdd(3).*rand(psmin,1),...
        T1mean(1)-2*T1stdd(1) + 4*T1stdd(1).*rand(psmin,1),...
        T1mean(2)-2*T1stdd(2) + 4*T1stdd(2).*rand(psmin,1),...
        T1mean(3)-2*T1stdd(3) + 4*T1stdd(3).*rand(psmin,1),...
        T2mean(1)-2*T2stdd(1) + 4*T2stdd(1).*rand(psmin,1),...
        T2mean(2)-2*T2stdd(2) + 4*T2stdd(2).*rand(psmin,1),...
        T2mean(3)-2*T2stdd(3) + 4*T2stdd(3).*rand(psmin,1)];

plotvec=zeros([1,2*length(md)]);
plotvec(1:2:end)=md-3*stdd; plotvec(2:2:end)=md+3*stdd;
plotcorr(mu_unif,plotvec,[]);
saveas(gcf,'Figures/priorval_muunif','png');

zunifcov=cov(z_unif);
zunifcorr=corr(z_unif);

[logZ,nest_samples,post_samples_gauss_gm,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,1,Nlive,Nmcmc,'gaussian');
[logZ,nest_samples,post_samples_gauss_wm,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,2,Nlive,Nmcmc,'gaussian');
[logZ,nest_samples,post_samples_gauss_csf,prior]=mn_qalasnp_04242019_covtest(TDin,tisinput,acqparam,3,Nlive,Nmcmc,'gaussian');
psmin=min([size(post_samples_gauss_gm,1),size(post_samples_gauss_wm,1),size(post_samples_gauss_csf,1)]);
z_gauss=[post_samples_gauss_gm(1:psmin,1:3),post_samples_gauss_wm(1:psmin,1:3),post_samples_gauss_csf(1:psmin,1:3)];

md=[tisinput(1,1:3),tisinput(3,1:3),tisinput(5,1:3)];
stdd=[tisinput(2,1:3),tisinput(4,1:3),tisinput(6,1:3)];
plotvec=zeros([1,2*length(md)]);
plotvec(1:2:end)=mean(z_gauss,1)-3*stdd; plotvec(2:2:end)=mean(z_gauss,1)+3*stdd;
plotcorr(z_gauss,plotvec,[]);
saveas(gcf,'Figures/priorval_zgauss','png');

mu_gauss=[normrnd(M0mean(1),M0stdd(1),[psmin,1]),...
        normrnd(M0mean(2),M0stdd(2),[psmin,1]),...
        normrnd(M0mean(3),M0stdd(3),[psmin,1]),...
        normrnd(T1mean(1),T1stdd(1),[psmin,1]),...
        normrnd(T1mean(2),T1stdd(2),[psmin,1]),...
        normrnd(T1mean(3),T1stdd(3),[psmin,1]),...
        normrnd(T2mean(1),T2stdd(1),[psmin,1]),...
        normrnd(T2mean(2),T2stdd(2),[psmin,1]),...
        normrnd(T2mean(3),T2stdd(3),[psmin,1])];

plotvec=zeros([1,2*length(md)]);
plotvec(1:2:end)=md-3*stdd; plotvec(2:2:end)=md+3*stdd;
plotcorr(mu_gauss,plotvec,[]);
saveas(gcf,'Figures/priorval_mugauss','png');

zgausscov=cov(z_gauss);
zgausscorr=corr(z_gauss);


%% Machine noise

segimgPath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/invivo_segimg.nii.gz';
segimg=load_untouch_nii(segimgPath);
segimg=segimg.img;

m0_roi_scan1{29}=M0pred1(segimg==29);
t1_roi_scan1{29}=T1pred1(segimg==29);
t2_roi_scan1{29}=T2pred1(segimg==29);

numeas=[m0_roi_scan1{29},t1_roi_scan1{29},t2_roi_scan1{29}];
numeascov=cov(numeas);
numeascorr=corr(numeas);




axthresh=[0,4,0,4,0,5];
plotthresh=[0,4,0,3.99,0,5];

nparam=3;
NumTicks = 2;
fontsz=6;
labels = {'M0_{background}','T1_{background}','T2_{background}'};
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nparam
    for j = 1:nparam
        if i == j
            subplot(nparam,nparam,nparam*(i-1)+j);

            tmp=numeas(:,i);
            histogram(tmp(and(tmp>=plotthresh(2*i-1),tmp<=plotthresh(2*i))),20);
            yl=get(gca,'YLim');
            axis([axthresh(2*i-1),axthresh(2*i),yl(1),yl(2)]);

            
            xl = get(gca,'XLim');
            set(gca,'XTick',linspace(xl(1),xl(2),NumTicks),'FontSize',fontsz);
            yl=get(gca,'YLim');
            set(gca,'YTick',linspace(yl(1),yl(2),NumTicks),'FontSize',fontsz);
        end
        
        if i < j
            subplot(nparam, nparam, nparam*(i-1)+j)%count)
            
            tmpi=numeas(:,i); tmpj=numeas(:,j);
            indi=and(tmpi>=plotthresh(2*i-1),tmpi<=plotthresh(2*i));
            indj=and(tmpj>=plotthresh(2*j-1),tmpj<=plotthresh(2*j));

            
            histogram2(tmpi(and(indi,indj)), tmpj(and(indi,indj)), [20,20], 'FaceColor','flat');
            view(30,60);
            
            yl=ylim; zl=zlim;
            axis([axthresh(2*i-1),axthresh(2*i),yl(1),yl(2),zl(1),zl(2)]);

            xl=xlim;
            axis([xl(1),xl(2),axthresh(2*j-1),axthresh(2*j),zl(1),zl(2)]);

            
            xl = get(gca,'XLim');
            set(gca,'XTick',linspace(xl(1),xl(2),NumTicks),'FontSize',fontsz);
            yl=get(gca,'YLim');
            set(gca,'YTick',linspace(yl(1),yl(2),NumTicks),'FontSize',fontsz);
            zl=get(gca,'ZLim');
            set(gca,'ZTick',linspace(zl(1),zl(2),NumTicks),'FontSize',fontsz);
        end
        
        if i > j
            subplot(nparam, nparam, nparam*(i-1)+j)%count)
                        
            tmpi=numeas(:,i); tmpj=numeas(:,j);
            indi=and(tmpi>=plotthresh(2*i-1),tmpi<=plotthresh(2*i));
            indj=and(tmpj>=plotthresh(2*j-1),tmpj<=plotthresh(2*j));
            
            scatter(tmpj(and(indi,indj)), tmpi(and(indi,indj)),'r.'); 
            
            yl=ylim;
            axis([axthresh(2*j-1),axthresh(2*j),yl(1),yl(2)]);

            xl=xlim;
            axis([xl(1),xl(2),axthresh(2*i-1),axthresh(2*i)]);

            xl = get(gca,'XLim');
            set(gca,'XTick',linspace(xl(1),xl(2),NumTicks),'FontSize',fontsz);
            yl=get(gca,'YLim');
            set(gca,'YTick',linspace(yl(1),yl(2),NumTicks),'FontSize',fontsz);
        end
        
        if j == 1
            ylabel(labels{i},'FontSize',fontsz+6);
        end
        if i == 3
            xlabel(labels{j},'FontSize',fontsz+6);
        end
    end
end
