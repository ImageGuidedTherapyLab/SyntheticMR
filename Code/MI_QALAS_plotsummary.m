function [plotvarinit] = MI_QALAS_plotsummary(tconoverridein,ttotalin,pdvvalin)

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
        tdtime{iii,jjj}(:)=squeeze(sum(pspacesave(2:5,:,:),1));
        scantime{iii,jjj}(:)=(tdtime{iii,jjj}(:)+TE_T2prep+TDpT2+nacq*Tacq+TDinv)*ceil(sum(subsmplconstrain(:))/100);
        
        %% NAN Handling
        M0save(isnan(M0save))=0; T1save(isnan(T1save))=0; T2save(isnan(T2save))=0;
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
    end
end

%% Main Results Figure
figure; hold on;
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

figure; plot(plotvarinit(1:4:end,3),plotvarinit(1:4:end,9),'o');
hold on; plot(plotvarinit(2:4:end,3),plotvarinit(2:4:end,9),'x');
hold on; plot(plotvarinit(3:4:end,3),plotvarinit(3:4:end,9),'+');
hold on; plot(plotvarinit(4:4:end,3),plotvarinit(4:4:end,9),'s');
hold on; plot(plotvarinit(4:4:end,2),plotvarinit(4:4:end,8),'ro');
xlabel('TD Sum (s)'); ylabel('Overall Error'); title('White Matter T1 Reconstruction');
legend('100% Opt','70% Opt','50% Opt','25% Opt','Init','Location','NorthWest');

figure; plot(plotvarinit(1:4:end,5),plotvarinit(1:4:end,9),'bo');
hold on; plot(plotvarinit(2:4:end,5),plotvarinit(2:4:end,9),'bx');
hold on; plot(plotvarinit(3:4:end,5),plotvarinit(3:4:end,9),'b+');
hold on; plot(plotvarinit(4:4:end,5),plotvarinit(4:4:end,9),'bs');
hold on; plot(plotvarinit(1:4:end,4),plotvarinit(1:4:end,8),'ro');
hold on; plot(plotvarinit(2:4:end,4),plotvarinit(2:4:end,8),'rx');
hold on; plot(plotvarinit(3:4:end,4),plotvarinit(3:4:end,8),'r+');
hold on; plot(plotvarinit(4:4:end,4),plotvarinit(4:4:end,8),'rs');
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title('White Matter T1 Reconstruction');
legend('100% Opt','70% Opt','50% Opt','25% Opt','100% Init','70% Init','50% Init','25% Init','Location','NorthEast');

figure; plot(plotvarinit(1:4:end,3),plotvarinit(1:4:end,7),'o');
hold on; plot(plotvarinit(2:4:end,3),plotvarinit(2:4:end,7),'x');
hold on; plot(plotvarinit(3:4:end,3),plotvarinit(3:4:end,7),'+');
hold on; plot(plotvarinit(4:4:end,3),plotvarinit(4:4:end,7),'s');
hold on; plot(plotvarinit(4:4:end,2),plotvarinit(4:4:end,6),'ro');
xlabel('TD Sum (s)'); ylabel('Overall Error'); title('White Matter M0 Reconstruction');
legend('100% Opt','70% Opt','50% Opt','25% Opt','Init','Location','NorthWest');

figure; plot(plotvarinit(1:4:end,5),plotvarinit(1:4:end,7),'bo');
hold on; plot(plotvarinit(2:4:end,5),plotvarinit(2:4:end,7),'bx');
hold on; plot(plotvarinit(3:4:end,5),plotvarinit(3:4:end,7),'b+');
hold on; plot(plotvarinit(4:4:end,5),plotvarinit(4:4:end,7),'bs');
hold on; plot(plotvarinit(1:4:end,4),plotvarinit(1:4:end,6),'ro');
hold on; plot(plotvarinit(2:4:end,4),plotvarinit(2:4:end,6),'rx');
hold on; plot(plotvarinit(3:4:end,4),plotvarinit(3:4:end,6),'r+');
hold on; plot(plotvarinit(4:4:end,4),plotvarinit(4:4:end,6),'rs');
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title('White Matter M0 Reconstruction');
legend('100% Opt','70% Opt','50% Opt','25% Opt','100% Init','70% Init','50% Init','25% Init','Location','NorthEast');

%% Data Matrix


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