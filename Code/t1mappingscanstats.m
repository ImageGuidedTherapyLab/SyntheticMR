
%% QALAS Data
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


%% T1 Mapping Data
for iii=1:70
    t1img(:,:,iii)=dicomread(sprintf('t1mappingscan/Image#0%02d.dcm',iii));
end

t1img=4.3217.*double(t1img)-80.504;

nii=make_nii(t1img);
save_nii(nii,'t1mappingscan/t1img.nii');



t1seg=load_untouch_nii('t1mappingscan/t1imgseg.nii.gz');
t1seg=t1seg.img;

for iii=1:42
    t1mean(iii)=nanmean(t1img(t1seg==iii));
    t1std(iii)=nanstd(t1img(t1seg==iii));
    t1n(iii)=sum(t1seg(:)==iii);
end

% figure; errorbar(t1mean(1:14),t1std(1:14),'o');
% xlabel('PD Element'); ylabel('T1 (ms)');
% figure; errorbar(t1mean(15:28),t1std(15:28),'o');
% xlabel('T1 Element'); ylabel('T1 (ms)');
% figure; errorbar(t1mean(29:42),t1std(29:42),'o');
% xlabel('T2 Element'); ylabel('T1 (ms)');

%% Manual T1 Values
manT1(1:14)=0;
manT1(15:28)=[1989,1454,984.1,706,496.7,351.5,247.13,175.3,125.9,89.0,62.7,44.54,30.84,21.719];
manT1(29:42)=[2480,2173,1907,1604,1332,1044,801.7,608.6,458.4,336.5,224.2,176.6,126.9,90.0];
manT1s(1:14)=0;
manT1s(15:28)=[1.0,2.5,0.33,1.5,0.41,0.91,0.086,0.11,0.33,0.17,0.13,0.09,0.016,0.004];
manT1s(29:42)=[10.8,14.7,10.3,7.2,0.8,3.2,1.7,1.03,0.33,0.18,0.09,0.09,0.03,0.05];

%% Plots
t1mean(t1mean<0)=0;
figure; hold on;
errorbar(t1mean(1:14),t1std(1:14),'bo');
errorbar(nanmean(T1v_1(1:14,:),2)*1000,T1std_1(1:14)*1000,'rx');
xlabel('PD Element'); ylabel('T1 (ms)');
legend('VFA T1','QALAS T1');
saveas(gcf,'Figures/t1scancompare_m0el','png');

figure; hold on;
errorbar(t1mean(15:28),t1std(15:28),'bo');
errorbar(nanmean(T1v_1(15:28,:),2)*1000,T1std_1(15:28)*1000,'rx');
errorbar(manT1(15:28)',manT1s(15:28)','gs');
xlabel('T1 Element'); ylabel('T1 (ms)');
legend('VFA T1','QALAS T1','True T1');
saveas(gcf,'Figures/t1scancompare_t1el','png');

figure; hold on;
errorbar(t1mean(29:42),t1std(29:42),'bo');
errorbar(nanmean(T1v_1(29:42,:),2)*1000,T1std_1(29:42)*1000,'rx');
errorbar(manT1(29:42)',manT1s(29:42)','gs');
xlabel('T2 Element'); ylabel('T1 (ms)');
legend('VFA T1','QALAS T1','True T1');
saveas(gcf,'Figures/t1scancompare_t2el','png');
