
function [] = MI_QALAS_plotfigs(tconoverride,ttotal,pdvval)

fmsym={'o','x','+','*','s','.','d','^'};
fmcol={'r','g','b','c','m','y','k','w'};

%% Load
load(sprintf('results/optresults_subsamp_%f_%f_%f.mat',tconoverride,ttotal,pdvval));
load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_goldstandards_%f_%f_%f.mat',tconoverride,ttotal,pdvval));
load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_%f.mat',tconoverride,ttotal,pdvval));
opt_hist_table=readtable(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/opt_history_%f_%f_%f.txt',tconoverride,ttotal,pdvval));

%% NAN Handling
M0save(isnan(M0save))=0; T1save(isnan(T1save))=0; T2save(isnan(T2save))=0;
goldstandardM0(isnan(goldstandardM0))=0; goldstandardT1(isnan(goldstandardT1))=0; goldstandardT2(isnan(goldstandardT2))=0;
synthdataM0(isnan(synthdataM0))=0; synthdataT1(isnan(synthdataT1))=0; synthdataT2(isnan(synthdataT2))=0;

%% Overall Error
for iii=1:size(M0save,3)
    M0OR(iii)=norm(M0save(:,:,iii)-goldstandardM0,2)/norm(goldstandardM0,2);
    T1OR(iii)=norm(T1save(:,:,iii)-goldstandardT1,2)/norm(goldstandardT1,2);
    T2OR(iii)=norm(T2save(:,:,iii)-goldstandardT2,2)/norm(goldstandardT2,2);
end

%% Tissue Overall Error
for jjj=1:3
    for iii=1:size(M0save,3)
        M0ORtis(iii,jjj)=norm((M0save(:,:,iii)-goldstandardM0).*(materialID==jjj),2)/norm(goldstandardM0.*(materialID==jjj),2);
        T1ORtis(iii,jjj)=norm((T1save(:,:,iii)-goldstandardT1).*(materialID==jjj),2)/norm(goldstandardT1.*(materialID==jjj),2);
        T2ORtis(iii,jjj)=norm((T2save(:,:,iii)-goldstandardT2).*(materialID==jjj),2)/norm(goldstandardT2.*(materialID==jjj),2);
    end
end

%% Relative Error
M0relerr=abs(M0save-repmat(goldstandardM0,[1,1,size(M0save,3)]))./abs(repmat(goldstandardM0,[1,1,size(M0save,3)]));
M0relerr(isnan(M0relerr))=0;
T1relerr=abs(T1save-repmat(goldstandardT1,[1,1,size(T1save,3)]))./abs(repmat(goldstandardT1,[1,1,size(T1save,3)]));
T1relerr(isnan(T1relerr))=0;
T2relerr=abs(T2save-repmat(goldstandardT2,[1,1,size(T2save,3)]))./abs(repmat(goldstandardT2,[1,1,size(T2save,3)]));
T2relerr(isnan(T2relerr))=0;
M0relerrn=abs(M0save-repmat(synthdataM0,[1,1,size(M0save,3)]))./abs(repmat(synthdataM0,[1,1,size(M0save,3)]));
M0relerrn(isnan(M0relerrn))=0;
T1relerrn=abs(T1save-repmat(synthdataT1,[1,1,size(T1save,3)]))./abs(repmat(synthdataT1,[1,1,size(T1save,3)]));
T1relerrn(isnan(T1relerrn))=0;
T2relerrn=abs(T2save-repmat(synthdataT2,[1,1,size(T2save,3)]))./abs(repmat(synthdataT2,[1,1,size(T2save,3)]));
T2relerrn(isnan(T2relerrn))=0;

% M0RE=M0relerr.*(repmat(materialID,[1,1,size(M0relerr,3)])==1);

for jjj=1:3
    for iii=1:size(M0relerr,3)
        tmp=M0relerr(:,:,iii).*(materialID==jjj);
%         tmp=tmp(~isnan(tmp));
        M0REmedian(iii,jjj)=median(tmp(tmp~=0));
        tmp=T1relerr(:,:,iii).*(materialID==jjj);
%         tmp=tmp(~isnan(tmp));
        T1REmedian(iii,jjj)=median(tmp(tmp~=0));
        tmp=T2relerr(:,:,iii).*(materialID==jjj);
%         tmp=tmp(~isnan(tmp));
        T2REmedian(iii,jjj)=median(tmp(tmp~=0));
        tmp=M0relerrn(:,:,iii).*(materialID==jjj);
%         tmp=tmp(~isnan(tmp));
        M0REnmedian(iii,jjj)=median(tmp(tmp~=0));
        tmp=T1relerrn(:,:,iii).*(materialID==jjj);
%         tmp=tmp(~isnan(tmp));
        T1REnmedian(iii,jjj)=median(tmp(tmp~=0));
        tmp=T2relerrn(:,:,iii).*(materialID==jjj);
%         tmp=tmp(~isnan(tmp));
        T2REnmedian(iii,jjj)=median(tmp(tmp~=0));
    end
end

for iii=1:3
    figure;plot(M0REmedian(:,iii),'o');axis([0,150,0,5]);
    figure;plot(T1REmedian(:,iii),'o');axis([0,150,0,5]);
    figure;plot(T2REmedian(:,iii),'o');axis([0,150,0,5]);
% figure;plot(M0REnmedian(:,1),'o');
end

% M0RE=M0RE(M0RE~=0);
% 
% M0RE=M0relerr(repmat(materialID,[1,1,size(M0relerr,3)])==1);

%% Plot Figures
% Parameter Maps and Error Maps
figure;
imagesc(M0save(:,:,1),[0,1]); colormap('gray'); colorbar;
title('M0 Recon First Iteration');
figure;
imagesc(M0save(:,:,end),[0,1]); colormap('gray'); colorbar;
title('M0 Recon Final Iteration');

figure;
imagesc(abs(M0save(:,:,1)-goldstandardM0)./abs(goldstandardM0),[0,1]); colormap('gray'); colorbar;
title('M0 Recon Relative Error First Iteration');
figure;
imagesc(abs(M0save(:,:,end)-goldstandardM0)./abs(goldstandardM0),[0,1]); colormap('gray'); colorbar;
title('M0 Recon Relative Error Final Iteration');

figure;
imagesc(T1save(:,:,1),[0,1]); colormap('gray'); colorbar;
title('T1 Recon First Iteration');
figure;
imagesc(T1save(:,:,end),[0,1]); colormap('gray'); colorbar;
title('T1 Recon Final Iteration');

figure;
imagesc(abs(T1save(:,:,1)-goldstandardT1)./abs(goldstandardT1),[0,1]); colormap('gray'); colorbar;
title('T1 Recon Relative Error First Iteration');
figure;
imagesc(abs(T1save(:,:,end)-goldstandardT1)./abs(goldstandardT1),[0,1]); colormap('gray'); colorbar;
title('T1 Recon Relative Error Final Iteration');

try
figure;
imagesc(T2save(:,:,1),[0,1]); colormap('gray'); colorbar;
title('T2 Recon First Iteration');
figure;
imagesc(T2save(:,:,end),[0,1]); colormap('gray'); colorbar;
title('T2 Recon Final Iteration');

figure;
imagesc(abs(T2save(:,:,1)-goldstandardT2)./abs(goldstandardT2),[0,1]); colormap('gray'); colorbar;
title('T2 Recon Relative Error First Iteration');
figure;
imagesc(abs(T2save(:,:,end)-goldstandardT2)./abs(goldstandardT2),[0,1]); colormap('gray'); colorbar;
title('T2 Recon Relative Error Final Iteration');
catch
end

% Reconstruction Statistics During Optimization
figure;
hold on;
for iii=1:3
    for jjj=1:3
        hold on; plot(squeeze(varsave(jjj,iii,:)),[fmsym{iii},fmcol{jjj}]);
    end
end
legend('M0 WM','M0 GM','M0 CSF','T1 WM','T1 GM','T1 CSF','T2 WM','T2 GM','T2 CSF');
xlabel('Mutual Information'); ylabel('Variance');

figure;
hold on;
for iii=1:3
    for jjj=1:3
        hold on; plot(squeeze(varsave(jjj,iii,:)),[fmsym{iii},fmcol{jjj}]);
    end
end
legend('M0 WM','M0 GM','M0 CSF','T1 WM','T1 GM','T1 CSF','T2 WM','T2 GM','T2 CSF'); axis([0 150 0 20]);
xlabel('Mutual Information'); ylabel('Variance');

figure;
hold on;
for iii=1:3
    for jjj=1:3
        hold on; plot(squeeze(meansave(jjj,iii,:)),[fmsym{iii},fmcol{jjj}]);
    end
end
legend('M0 WM','M0 GM','M0 CSF','T1 WM','T1 GM','T1 CSF','T2 WM','T2 GM','T2 CSF');
xlabel('Mutual Information'); ylabel('Mean');

figure;
hold on;
for iii=1:3
    for jjj=1:3
        hold on; plot(squeeze(meansave(jjj,iii,:)),[fmsym{iii},fmcol{jjj}]);
    end
end
legend('M0 WM','M0 GM','M0 CSF','T1 WM','T1 GM','T1 CSF','T2 WM','T2 GM','T2 CSF'); axis([0 150 0 20]);
xlabel('Mutual Information'); ylabel('Mean');

figure;
hold on;
for iii=1:3
    for jjj=1:3
        hold on; plot(squeeze(mediansave(jjj,iii,:)),[fmsym{iii},fmcol{jjj}]);
    end
end
legend('M0 WM','M0 GM','M0 CSF','T1 WM','T1 GM','T1 CSF','T2 WM','T2 GM','T2 CSF'); 
xlabel('Mutual Information'); ylabel('Median');

figure;
hold on;
for iii=1:3
    for jjj=1:3
        plot(squeeze(mediansave(jjj,iii,:)),[fmsym{iii},fmcol{jjj}]);
    end
end
legend('M0 WM','M0 GM','M0 CSF','T1 WM','T1 GM','T1 CSF','T2 WM','T2 GM','T2 CSF'); axis([0 150 0 20]);
xlabel('Mutual Information'); ylabel('Median');

% Parameter Performance During Optimization
figure;
plot(-squeeze(MIsave),'o');
xlabel('Iteration'); ylabel('Mutual Information');
title('Mutual Information Objective Function');

figure;
plot(squeeze(pspacesave(iii,1,:)),'o');
xlabel('Iteration'); ylabel('Flip Angle (deg)');
title('Flip Angle Optimization');

figure;
hold on;
for iii=2:5
    plot(squeeze(pspacesave(iii,1,:)),[fmsym{iii},fmcol{iii}]);
end
legend('TD2','TD3','TD4','TD5');
xlabel('Iteration'); ylabel('Delay Time (s)');
title('Delay Time Optimization');

figure;
plot(squeeze(pspacesave(2,1,:)+pspacesave(3,1,:)+pspacesave(4,1,:)+pspacesave(5,1,:)),'o');
xlabel('Iteration'); ylabel('Total Delay Time');
title('Total Delay Time');

end