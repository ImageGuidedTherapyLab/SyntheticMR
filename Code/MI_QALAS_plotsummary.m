function [plotvarinit,MIstats,pspacestats,tisrecon,tisreconerr,tisreconrelerr,gssave,sdsave] = MI_QALAS_plotsummary(tconoverridein,ttotalin,pdvvalin)

dirname=sprintf('MIQALASplotsummary_%i_%2.2f-%2.2f_%2.2f-%2.2f',tconoverridein,ttotalin(1),ttotalin(end),pdvvalin(1),pdvvalin(end));
resultspathname=['/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/',dirname,'/'];
if exist(resultspathname,'dir')~=1
    mkdir(resultspathname);
end

fmsym={'o','x','+','*','s','.','d','^'};
fmcol={'r','g','b','c','m','y','k','w'};

for iii=1:length(ttotalin)
    for jjj=1:length(pdvvalin)
        %% Load
        load(sprintf('results/optresults_subsamp_%f_%f_%f.mat',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_goldstandards_%f_%f_%f.mat',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_%f.mat',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        opt_hist_table=readtable(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/opt_history_%f_%f_%f.txt',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        
        %% Total Scan Time
        if pdarg(3)==-1
            subsmplconstrain=ones([size(materialID,1),size(materialID,2)]);
        else
            subsmplconstrain=bart(sprintf('poisson -Y %i -Z %i -y %f -z %f -V %f',size(materialID,1),size(materialID,2),pdarg(1),pdarg(2),pdarg(3)));
        end
        %             tconstrain=ttotal*60/ceil(sum(subsmplconstrain(:))/100)-TE_T2prep-TDpT2-nacq*Tacq-TDinv; % seconds
        subpct{iii,jjj}=sum(subsmplconstrain(:))/numel(subsmplconstrain);
        tmpdimvec=1:size(pspacesave,1);
        tdtime{iii,jjj}(:)=squeeze(sum(pspacesave(tmpdimvec(strncmp('TD',pspacelabels,2)),:,:),1));
        scantime{iii,jjj}(:)=(tdtime{iii,jjj}(:)+TE_T2prep+TDpT2+nacq*Tacq+TDinv)*ceil(sum(subsmplconstrain(:))/100);
        
        %% NAN Handling
        M0save(isnan(M0save))=0; T1save(isnan(T1save))=0; T2save(isnan(T2save))=0; T2save=real(T2save);
%         M0save(M0save>10)=0; T1save(T1save>10)=0; T2save(T2save>10)=0;
        goldstandardM0(isnan(goldstandardM0))=0; goldstandardT1(isnan(goldstandardT1))=0; goldstandardT2(isnan(goldstandardT2))=0;
        synthdataM0(isnan(synthdataM0))=0; synthdataT1(isnan(synthdataT1))=0; synthdataT2(isnan(synthdataT2))=0;
        
        %% Overall Error
        for kkk=1:size(M0save,3)
            M0OR{iii,jjj}(kkk)=norm(M0save(:,:,kkk)-goldstandardM0,2)/norm(goldstandardM0,2);
            T1OR{iii,jjj}(kkk)=norm(T1save(:,:,kkk)-goldstandardT1,2)/norm(goldstandardT1,2);
            T2OR{iii,jjj}(kkk)=norm(T2save(:,:,kkk)-goldstandardT2,2)/norm(goldstandardT2,2);
        end
        
        %% Tissue Overall Error
        for lll=1:3
            for kkk=1:size(M0save,3)
                M0ORtis{iii,jjj}(kkk,lll)=norm((M0save(:,:,kkk)-goldstandardM0).*(materialID==lll),2)/norm(goldstandardM0.*(materialID==lll),2);
                T1ORtis{iii,jjj}(kkk,lll)=norm((T1save(:,:,kkk)-goldstandardT1).*(materialID==lll),2)/norm(goldstandardT1.*(materialID==lll),2);
                T2ORtis{iii,jjj}(kkk,lll)=norm((T2save(:,:,kkk)-goldstandardT2).*(materialID==lll),2)/norm(goldstandardT2.*(materialID==lll),2);
            end
        end
        
        %% Example Images
        tisrecon{iii,jjj,1}=cat(3,M0save(:,:,1),M0save(:,:,end));
        tisrecon{iii,jjj,2}=cat(3,T1save(:,:,1),T1save(:,:,end));
        tisrecon{iii,jjj,3}=cat(3,T2save(:,:,1),T2save(:,:,end));
        tisreconerr{iii,jjj,1}=cat(3,abs(M0save(:,:,1)-goldstandardM0),abs(M0save(:,:,end)-goldstandardM0));
        tisreconerr{iii,jjj,2}=cat(3,abs(T1save(:,:,1)-goldstandardT1),abs(T1save(:,:,end)-goldstandardT1));
        tisreconerr{iii,jjj,3}=cat(3,abs(T2save(:,:,1)-goldstandardT2),abs(T2save(:,:,end)-goldstandardT2));
        tisreconrelerr{iii,jjj,1}=cat(3,abs(M0save(:,:,1)-goldstandardM0)./goldstandardM0,abs(M0save(:,:,end)-goldstandardM0)./goldstandardM0);
        tisreconrelerr{iii,jjj,2}=cat(3,abs(T1save(:,:,1)-goldstandardT1)./goldstandardT1,abs(T1save(:,:,end)-goldstandardT1)./goldstandardT1);
        tisreconrelerr{iii,jjj,3}=cat(3,abs(T2save(:,:,1)-goldstandardT2)./goldstandardT2,abs(T2save(:,:,end)-goldstandardT2)./goldstandardT2);
        tisreconlabel{iii,jjj,1}=sprintf('%2.2f s - %2.2f pdv - M0',ttotalin(iii),pdvvalin(jjj));
        tisreconlabel{iii,jjj,2}=sprintf('%2.2f s - %2.2f pdv - T1',ttotalin(iii),pdvvalin(jjj));
        tisreconlabel{iii,jjj,3}=sprintf('%2.2f s - %2.2f pdv - T2',ttotalin(iii),pdvvalin(jjj));
        
        %% Gold Standards
        gssave{iii,jjj,1}=goldstandardM0; gssave{iii,jjj,2}=goldstandardT1; gssave{iii,jjj,3}=goldstandardT2;
        sdsave{iii,jjj,1}=synthdataM0; sdsave{iii,jjj,2}=synthdataT1; sdsave{iii,jjj,3}=synthdataT2;
        
        %% Optimization Stats
        MIstats{iii,jjj}=squeeze(MIsave);
        pspacestats{iii,jjj}=squeeze(pspacesave);
    end
end

%% Optimization Stats
titlename='Mutual Information';
figure; hold on;
for jjj=1:length(pdvvalin)
    for iii=1:length(ttotalin)
        plot(MIstats{iii,jjj},'-o');
    end
end
xlabel('Iteration'); ylabel('Mutual Information'); title(titlename);
saveas(gcf,[resultspathname,titlename],'png');

for kkk=1:size(pspacestats{1},1)
    titlename=pspacelabels{kkk};
    figure; hold on;
    for jjj=1:length(pdvvalin)
        for iii=1:length(ttotalin);
            plot(pspacestats{iii,jjj}(kkk,:),'o');
        end
    end
    xlabel('Iteration'); ylabel(pspacelabels{kkk}); title(titlename);
    saveas(gcf,[resultspathname,titlename],'png');
end

% titlename='Flip Angle';
% figure; hold on;
% for jjj=1:length(pdvvalin)
%     for iii=1:length(ttotalin);
%         plot(pspacestats{iii,jjj}(1,:),'o');
%     end
% end
% xlabel('Iteration'); ylabel('Flip Angle'); title(titlename);
% saveas(gcf,[resultspathname,titlename],'png');
% 
% titlename='TD1';
% figure; hold on;
% for jjj=1:length(pdvvalin)
%     for iii=1:length(ttotalin);
%         plot(pspacestats{iii,jjj}(2,:),'o');
%     end
% end
% xlabel('Iteration'); ylabel('TD1'); title(titlename);
% saveas(gcf,[resultspathname,titlename],'png');
% 
% titlename='TD2';
% figure; hold on;
% for jjj=1:length(pdvvalin)
%     for iii=1:length(ttotalin);
%         plot(pspacestats{iii,jjj}(3,:),'o');
%     end
% end
% xlabel('Iteration'); ylabel('TD2'); title(titlename);
% saveas(gcf,[resultspathname,titlename],'png');
% 
% titlename='TD3';
% figure; hold on;
% for jjj=1:length(pdvvalin)
%     for iii=1:length(ttotalin);
%         plot(pspacestats{iii,jjj}(4,:),'o');
%     end
% end
% xlabel('Iteration'); ylabel('TD3'); title(titlename);
% saveas(gcf,[resultspathname,titlename],'png');
% 
% titlename='TD4';
% figure; hold on;
% for jjj=1:length(pdvvalin)
%     for iii=1:length(ttotalin);
%         plot(pspacestats{iii,jjj}(5,:),'o');
%     end
% end
% xlabel('Iteration'); ylabel('TD4'); title(titlename);
% saveas(gcf,[resultspathname,titlename],'png');

%% Example Images
tismax=[1.5,5,1];
for kkk=1:3
    for jjj=1:length(pdvvalin)
        for iii=1:length(ttotalin)
            titlename=['Recon Init ',tisreconlabel{iii,jjj,kkk}];
            figure; imagesc(tisrecon{iii,jjj,kkk}(:,:,1),[0,tismax(kkk)]); colormap('gray'); colorbar;
            title(titlename);
            saveas(gcf,[resultspathname,titlename,'.png']);
            close;
            titlename=['Recon Opt ',tisreconlabel{iii,jjj,kkk}];
            figure; imagesc(tisrecon{iii,jjj,kkk}(:,:,2),[0,tismax(kkk)]); colormap('gray'); colorbar;
            title(titlename);
            saveas(gcf,[resultspathname,titlename,'.png']);
            close;
            titlename=['Recon Error Init ',tisreconlabel{iii,jjj,kkk}];
            figure; imagesc(tisreconerr{iii,jjj,kkk}(:,:,1),[0,tismax(kkk)]); colormap('gray'); colorbar;
            title(titlename);
            saveas(gcf,[resultspathname,titlename,'.png']);
            close;
            titlename=['Recon Error Opt ',tisreconlabel{iii,jjj,kkk}];
            figure; imagesc(tisreconerr{iii,jjj,kkk}(:,:,2),[0,tismax(kkk)]); colormap('gray'); colorbar;
            title(titlename);
            saveas(gcf,[resultspathname,titlename,'.png']);
            close;
            titlename=['Recon Rel Error Init ',tisreconlabel{iii,jjj,kkk}];
            figure; imagesc(tisreconrelerr{iii,jjj,kkk}(:,:,1),[0,0.5]); colormap('gray'); colorbar;
            title(titlename);
            saveas(gcf,[resultspathname,titlename,'.png']);
            close;
            titlename=['Recon Rel Error Opt ',tisreconlabel{iii,jjj,kkk}];
            figure; imagesc(tisreconrelerr{iii,jjj,kkk}(:,:,2),[0,0.5]); colormap('gray'); colorbar;
            title(titlename);
            saveas(gcf,[resultspathname,titlename,'.png']);
            close;
        end
    end
end

%% Main Results Figure
plotvar=[scantime{1,1}(1),T1OR{1,1}(1)];
for iii=1:size(M0OR,1)
    for jjj=1:size(M0OR,2)
        plotvar=[plotvar;scantime{iii,jjj}(end),T1OR{iii,jjj}(end)];
        %         plot(scantime{iii,jjj}(1),T1OR{iii,jjj}(1));
        %         plot(scantime{iii,jjj}(end),T1OR{iii,jjj}(end));
    end
end

%% Main Results Tissue Figure
% White Matter (N=2)
plotvarinit=[]; plotvaropt=[];
for iii=1:size(M0ORtis,1)
    for jjj=1:size(M0ORtis,2)
        plotvarinit=[plotvarinit;subpct{iii,jjj},tdtime{iii,jjj}(1),tdtime{iii,jjj}(end),scantime{iii,jjj}(1),scantime{iii,jjj}(end),...
            M0ORtis{iii,jjj}(1,2),M0ORtis{iii,jjj}(end,2),T1ORtis{iii,jjj}(1,2),T1ORtis{iii,jjj}(end,2),T2ORtis{iii,jjj}(1,2),T2ORtis{iii,jjj}(end,2)];
        %         plotvaropt=[plotvaropt;scantime{iii,jjj}(end,2),M0ORtis{iii,jjj}(end,2),T1ORtis{iii,jjj}(end,2),T2ORtis{iii,jjj}(end,2)];
    end
end

%% Plot Figures
plotsym={'-o','--x',':+','-.s','-*','--d',':v','-.^','-<','-->',':p','-.h'};
% plotsym={'b-o','b--x','b:x','b-.s','b-*','r-o','r--x','r:x','r-.s','r-*'};
szgs=size(M0ORtis);

titlename='White Matter T1 Reconstruction';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,3),plotvarinit(iii:szgs(2):end,9),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarinit(iii:szgs(2):end,8),['r',plotsym{iii}]);
end
% plot(plotvarinit(szgs(1):szgs(2):end,2),plotvarinit(szgs(1):szgs(2):end,8),['r',plotsym{1}]);
% plot(plotvarinit(1:4:end,3),plotvarinit(1:4:end,9),'-o');
% plot(plotvarinit(2:4:end,3),plotvarinit(2:4:end,9),'--x');
% plot(plotvarinit(3:4:end,3),plotvarinit(3:4:end,9),':+');
% plot(plotvarinit(4:4:end,3),plotvarinit(4:4:end,9),'-.s');
% plot(plotvarinit(4:4:end,2),plotvarinit(4:4:end,8),'r-o');
xlabel('TD Sum (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter T1 Reconstruction - Acq Time';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,5),plotvarinit(iii:szgs(2):end,9),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,4),plotvarinit(iii:szgs(2):end,8),['r',plotsym{iii}]);
end
% plot(plotvarinit(1:4:end,5),plotvarinit(1:4:end,9),'b-o');
% plot(plotvarinit(2:4:end,5),plotvarinit(2:4:end,9),'b--x');
% plot(plotvarinit(3:4:end,5),plotvarinit(3:4:end,9),'b:+');
% plot(plotvarinit(4:4:end,5),plotvarinit(4:4:end,9),'b-.s');
% plot(plotvarinit(1:4:end,4),plotvarinit(1:4:end,8),'r-o');
% plot(plotvarinit(2:4:end,4),plotvarinit(2:4:end,8),'r--x');
% plot(plotvarinit(3:4:end,4),plotvarinit(3:4:end,8),'r:+');
% plot(plotvarinit(4:4:end,4),plotvarinit(4:4:end,8),'r-.s');
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title(titlename);
% legendkey={'100% Opt','70% Opt','50% Opt','25% Opt','100% Init','70% Init','50% Init','25% Init'};
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter M0 Reconstruction';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,3),plotvarinit(iii:szgs(2):end,7),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarinit(iii:szgs(2):end,6),['r',plotsym{iii}]);
end
% plot(plotvarinit(szgs(1):szgs(2):end,2),plotvarinit(szgs(1):szgs(2):end,6),['r',plotsym{1}]);
% plot(plotvarinit(1:4:end,3),plotvarinit(1:4:end,7),'-o');
% plot(plotvarinit(2:4:end,3),plotvarinit(2:4:end,7),'--x');
% plot(plotvarinit(3:4:end,3),plotvarinit(3:4:end,7),':+');
% plot(plotvarinit(4:4:end,3),plotvarinit(4:4:end,7),'-.s');
% plot(plotvarinit(4:4:end,2),plotvarinit(4:4:end,6),'r-o');
xlabel('TD Sum (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter M0 Reconstruction - Acq Time';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,5),plotvarinit(iii:szgs(2):end,7),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,4),plotvarinit(iii:szgs(2):end,6),['r',plotsym{iii}]);
end
% plot(plotvarinit(1:4:end,5),plotvarinit(1:4:end,7),'b-o');
% plot(plotvarinit(2:4:end,5),plotvarinit(2:4:end,7),'b--x');
% plot(plotvarinit(3:4:end,5),plotvarinit(3:4:end,7),'b:+');
% plot(plotvarinit(4:4:end,5),plotvarinit(4:4:end,7),'b-.s');
% plot(plotvarinit(1:4:end,4),plotvarinit(1:4:end,6),'r-o');
% plot(plotvarinit(2:4:end,4),plotvarinit(2:4:end,6),'r--x');
% plot(plotvarinit(3:4:end,4),plotvarinit(3:4:end,6),'r:+');
% plot(plotvarinit(4:4:end,4),plotvarinit(4:4:end,6),'r-.s');
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title(titlename);
% legendkey={'100% Opt','70% Opt','50% Opt','25% Opt','100% Init','70% Init','50% Init','25% Init'};
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter Mmeas';
szms=size(Mmeassave);
MmeasWM=Mmeassave;
MmeasWM(isnan(MmeasWM))=0;
MmeasWM=MmeasWM.*repmat(materialID==2,[1,1,szms(3),szms(4),szms(5)]);
MmeasWM=reshape(MmeasWM,[prod(szms(1:3)),szms(4),szms(5)]);
figure; hold on;
for jjj=1:size(MmeasWM,3)
    for iii=1:size(MmeasWM,2)
        plotpts=MmeasWM(:,iii,jjj);
        plot(iii*ones(size(plotpts)),plotpts,'o');
    end
end
xlabel('Acquisition Number'); ylabel('Mmeas'); title(titlename);
saveas(gcf,[resultspathname,titlename],'png');

%% Data Matrix
% Data array: MI, TD1,TD2,TD3,TD4, scan time, wm m0, wm t1, wm t2
dataarray=[];
dataarrayfull=[];
for iii=1:size(MIstats,1)
    for jjj=1:size(MIstats,2)
        dataarray=[dataarray;[MIstats{iii,jjj}([1,end]),pspacestats{iii,jjj}(:,[1,end])',...
            scantime{iii,jjj}([1,end])',M0ORtis{iii,jjj}([1,end],2),T1ORtis{iii,jjj}([1,end],2),...
            T2ORtis{iii,jjj}([1,end],2)]];
        dataarrayfull=[dataarrayfull;[MIstats{iii,jjj},pspacestats{iii,jjj}',...
            scantime{iii,jjj}',M0ORtis{iii,jjj}(:,2),T1ORtis{iii,jjj}(:,2),T2ORtis{iii,jjj}(:,2)]];
    end
end

titlename='Parameter Trends';
% tight_subplot(9,9,[.01,.03],[.1,.01],[.01,.01]);
for jjj=1:size(dataarray,2)
    for iii=1:size(dataarray,2)
        subplot(size(dataarray,2),size(dataarray,2),size(dataarray,2)*(jjj-1)+iii);
        if iii==jjj
            hist(dataarray(:,iii));
        elseif iii>jjj
%             histogram2(dataarray(:,iii),dataarray(:,jjj));
            plot(dataarrayfull(:,iii),dataarrayfull(:,jjj),'g+');
        else
            hold on;
            plot(dataarray(1:2:end,iii),dataarray(1:2:end,jjj),'bo');
            plot(dataarray(2:2:end,iii),dataarray(2:2:end,jjj),'rx');
        end
    end
end
saveas(gcf,[resultspathname,titlename],'png');

% nparam = 3;
%
% count = 1;
%
% while count <= nparam*nparam
%     for i = 1:nparam
%         for j = 1:nparam
%             if i == j
%                 subplot(nparam, nparam, count)
%                 histogram(post_samples(:,i), 50)
%
% %                 subplot(nparam, nparam, count)
% %                 scatter(post_samples(:,i), post_samples(:,nparam+1), '.')
%             end
%
%             if i < j
%                 subplot(nparam, nparam, count)
%                 histogram2(post_samples(:,i), post_samples(:,j), [50 50], 'FaceColor','flat');
%             end
%
%             if i > j
%                 subplot(nparam, nparam, count)
%                 scatter3(post_samples(:,i), post_samples(:,j), post_samples(:,nparam+1), '.')
%             end
%
%             count = count + 1;
%         end
%     end
% end
%

end