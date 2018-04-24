
function [] = MI_QALAS_plotfigs(tconoverride,ttotal,pdvval)

fmsym={'o','x','+','*','s','.','d','^'};
fmcol={'r','g','b','c','m','y','k','w'};

%% Load
load(sprintf('results/optresults_subsamp_%f_%f_%f.mat',tconoverride,ttotal,pdvval));
load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_goldstandards_%f_%f_%f.mat',tconoverride,ttotal,pdvval));
load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_%f.mat',tconoverride,ttotal,pdvval));
opt_hist_table=readtable(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/opt_history_%f_%f_%f.txt',tconoverride,ttotal,pdvval));

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