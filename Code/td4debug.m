function [tduniform,plotvarinit,plotmeangm,plotmeanwm,plotmedgm,plotmedwm,plotvargm,plotvarwm,Mmeasfullsave,Mmeassave,M0save,T1save,T2save,M0ORtis,T1ORtis,T2ORtis,meancell,mediancell,varcell] = td4debug(TD4in,subsampin,bartrecon)

tduniform=rand(TD4in,4);

flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
plotvarinit=[]; plotmeangm=[]; plotmeanwm=[]; plotmedgm=[]; plotmedwm=[]; plotvargm=[]; plotvarwm=[];
for iii=1:TD4in %length(TD4in)  %TD4 times
    for jjj=1:length(subsampin)   %subsample pct
        
        %% Tissue Properties
        load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/synthphantom_goldstandards.mat;
        %         goldstandardM0(isnan(goldstandardM0))=0; goldstandardT1(isnan(goldstandardT1))=0; goldstandardT2(isnan(goldstandardT2))=0;
        %         synthdataM0(isnan(synthdataM0))=0; synthdataT1(isnan(synthdataT1))=0; synthdataT2(isnan(synthdataT2))=0;
        %
        
        %% Total Scan Time
        if subsampin(jjj)==-1
            subsmplconstrain=ones([size(materialID,1),size(materialID,2)]);
        else
            subsmplconstrain=bart(sprintf('poisson -Y %i -Z %i -y %f -z %f -V %f',size(materialID,1),size(materialID,2),1,1,subsampin(jjj)));
        end
        %             tconstrain=ttotal*60/ceil(sum(subsmplconstrain(:))/100)-TE_T2prep-TDpT2-nacq*Tacq-TDinv; % seconds
        subpct{iii,jjj}=sum(subsmplconstrain(:))/numel(subsmplconstrain);
        %     tmpdimvec=1:size(pspacesave,1);
        %     tdtime{iii,jjj}(:)=squeeze(sum(pspacesave(tmpdimvec(strncmp('TD',pspacelabels,2)),:,:),1));
        %     scantime{iii,jjj}(:)=(tdtime{iii,jjj}(:)+TE_T2prep+TDpT2+nacq*Tacq+TDinv)*ceil(sum(subsmplconstrain(:))/100);
        
        %% Create synthetic QALAS measurements
        TD=tduniform(iii,:); %TD=[0.5,0.001,0.1,TD4in(iii)];
        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [Mmeasfull,Mmeas]=qalas(goldstandardM0,goldstandardM0,goldstandardT1,goldstandardT2,TR,TE_T2prep,flipAngle,nacq,dt);
        signu=.001;
        %         stdmapmeas=normrnd(0,signu,size(materialID));
        %         Mmeas=Mmeas+stdmapmeas;
        
        tdtime{iii,jjj}=sum(TD);
        scantime{iii,jjj}=(tdtime{iii,jjj}+TE_T2prep+TDpT2+nacq*Tacq+TDinv)*ceil(sum(subsmplconstrain(:))/100);
        
        %% Subsample synthetic measurements
        %             Mmeas(isnan(Mmeas))=0;
%         bartrecon=0;
        Mmeassub=Mmeas;
        if subsampin(jjj)~=-1
            subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i -C 10',size(materialID,1),size(materialID,2),1,1,subsampin(jjj)));
            Mmeassub(isnan(Mmeassub))=0;
            kmeas=bart('fft 3',Mmeassub);
            subsample=squeeze(subsmplmask(1,:,:));
            if bartrecon==1
                for reconind=1:size(kmeas,4)
                    Mmeassub(:,:,:,reconind)=bart('pics -l1 -r0.01 -S -w1',kmeas(:,:,:,reconind).*subsample,ones(size(kmeas(:,:,:,reconind))));
                end
            else
                subsample=repmat(subsample,[1,1,size(kmeas,3),size(kmeas,4)]);
                Mmeassub=bart('fft -i 3',kmeas.*subsample)*size(kmeas,4)/numel(kmeas);
            end
            Mmeassub=double(real(Mmeassub));
            Mmeassub(repmat(materialID,[1,1,size(Mmeassub,3),size(Mmeassub,4)])==0)=nan;
        end
        
        %% Reconstruct synthetic QALAS measurements
        % Optimization solution for M0 and T1 prediction
        smeas=size(Mmeassub);
        Mmeasvec=reshape(Mmeassub,[prod(smeas(1:3)),smeas(4:end)]);
        mmvsize=size(Mmeasvec,1);
        TR=TR; TE_T2prep=TE_T2prep; flipAngle=flipAngle; nacq=nacq; dt=dt;
        parfor zzz=1:mmvsize
            if sum(isnan(squeeze(Mmeasvec(zzz,:))))>0
                M0predvec(zzz)=0;%nan;
                T1predvec(zzz)=0;%nan;
                T2predvec(zzz)=0;%nan;
            else
                [M0predvec(zzz),T1predvec(zzz),T2predvec(zzz)]=qalasrecon(squeeze(Mmeasvec(zzz,:)),TR,TE_T2prep,flipAngle,nacq,dt);
            end
        end
        M0pred(:,:,:)=reshape(M0predvec,smeas(1:3));
        T1pred(:,:,:)=reshape(T1predvec,smeas(1:3));
        T2pred(:,:,:)=reshape(T2predvec,smeas(1:3));
        
        %         dtplot=cumsum(dt)';
        %         mmvec=reshape(Mmeassub,[size(Mmeassub,1)*size(Mmeassub,2)*size(Mmeassub,3),size(Mmeassub,4)]);
        %         mmvecfull=reshape(Mmeasfull,[size(Mmeasfull,1)*size(Mmeasfull,2)*size(Mmeasfull,3),size(Mmeasfull,4)]);
        %         for yyy=1:size(Mmeassub,4)
        %             gmmmplot(:,yyy)=mmvec(materialID(:)==1,yyy);
        %             wmmmplot(:,yyy)=mmvec(materialID(:)==2,yyy);
        %         end
        
        %         titlename=sprintf('Gray Matter, TD4=%2.2f, SS=%2.2f',TD4in(iii),subsampin(jjj));
        %         figure; hold on;
        %         for yyy=1:size(gmmmplot,1)
        %             plot(dtplot([2,6:2:end-1]),squeeze(gmmmplot(yyy,:)),'-o');
        %         end
        %         plot(dtplot,squeeze(sind(flipAngle)*Mmeasfull(16,58,1,:)),'r--s');
        %         xlabel('Time (s)'); ylabel('M'); title(titlename);
        %
        %         titlename=sprintf('White Matter, TD4=%2.2f, SS=%2.2f',TD4in(iii),subsampin(jjj));
        %         figure; hold on;
        %         for yyy=1:size(wmmmplot,1)
        %             plot(dtplot([2,6:2:end-1]),squeeze(wmmmplot(yyy,:)),'-o');
        %         end
        %         plot(dtplot,squeeze(sind(flipAngle)*Mmeasfull(49,36,1,:)),'r--s');
        %         xlabel('Time (s)'); ylabel('M'); title(titlename);
        
        Mmeasfullsave{iii,jjj}=Mmeasfull;
        Mmeassave{iii,jjj}=Mmeassub;
        M0save{iii,jjj}=M0pred;
        T1save{iii,jjj}=T1pred;
        T2save{iii,jjj}=T2pred;
        
        M0p1set=M0pred(:).*(materialID(:)==1); M0p1set=M0p1set(M0p1set~=0);
        M0p2set=M0pred(:).*(materialID(:)==2); M0p2set=M0p2set(M0p2set~=0);
        M0p3set=M0pred(:).*(materialID(:)==3); M0p3set=M0p3set(M0p3set~=0);
        T1p1set=T1pred(:).*(materialID(:)==1); T1p1set=T1p1set(T1p1set~=0);
        T1p2set=T1pred(:).*(materialID(:)==2); T1p2set=T1p2set(T1p2set~=0);
        T1p3set=T1pred(:).*(materialID(:)==3); T1p3set=T1p3set(T1p3set~=0);
        T2p1set=T2pred(:).*(materialID(:)==1); T2p1set=T2p1set(T2p1set~=0);
        T2p2set=T2pred(:).*(materialID(:)==2); T2p2set=T2p2set(T2p2set~=0);
        T2p3set=T2pred(:).*(materialID(:)==3); T2p3set=T2p3set(T2p3set~=0);
        
        M0varmeas=[quantile(M0p1set,0.975)-quantile(M0p1set,0.025)...
            quantile(M0p2set,0.975)-quantile(M0p2set,0.025)...
            quantile(M0p3set,0.975)-quantile(M0p3set,0.025)];
        T1varmeas=[quantile(T1p1set,0.975)-quantile(T1p1set,0.025)...
            quantile(T1p2set,0.975)-quantile(T1p2set,0.025)...
            quantile(T1p3set,0.975)-quantile(T1p3set,0.025)];
        T2varmeas=[quantile(T2p1set,0.975)-quantile(T2p1set,0.025)...
            quantile(T2p2set,0.975)-quantile(T2p2set,0.025)...
            quantile(T2p3set,0.975)-quantile(T2p3set,0.025)];
        varstats=[M0varmeas;T1varmeas;T2varmeas];
        
        meanstats=[mean(M0p1set),mean(M0p2set),mean(M0p3set);
            mean(T1p1set),mean(T1p2set),mean(T1p3set);
            mean(T2p1set),mean(T2p2set),mean(T2p3set)];
        
        medianstats=[median(M0p1set),median(M0p2set),median(M0p3set);
            median(T1p1set),median(T1p2set),median(T1p3set);
            median(T2p1set),median(T2p2set),median(T2p3set)];
        
        meancell{iii,jjj}=meanstats;
        mediancell{iii,jjj}=medianstats;
        varcell{iii,jjj}=varstats;
        
        goldstandardM0(isnan(goldstandardM0))=0; goldstandardT1(isnan(goldstandardT1))=0; goldstandardT2(isnan(goldstandardT2))=0;
        synthdataM0(isnan(synthdataM0))=0; synthdataT1(isnan(synthdataT1))=0; synthdataT2(isnan(synthdataT2))=0;
        
        %% Overall Error
        M0OR{iii,jjj}=norm(M0pred(:,:)-goldstandardM0,2)/norm(goldstandardM0,2);
        T1OR{iii,jjj}=norm(T1pred(:,:)-goldstandardT1,2)/norm(goldstandardT1,2);
        T2OR{iii,jjj}=norm(T2pred(:,:)-goldstandardT2,2)/norm(goldstandardT2,2);
        
        %% Tissue Overall Error
        for lll=1:3
            M0ORtis{iii,jjj}(lll)=norm((M0pred(:,:)-goldstandardM0).*(materialID==lll),2)/norm(goldstandardM0.*(materialID==lll),2);
            T1ORtis{iii,jjj}(lll)=norm((T1pred(:,:)-goldstandardT1).*(materialID==lll),2)/norm(goldstandardT1.*(materialID==lll),2);
            T2ORtis{iii,jjj}(lll)=norm((T2pred(:,:)-goldstandardT2).*(materialID==lll),2)/norm(goldstandardT2.*(materialID==lll),2);
        end
%         save td4debugresults.mat -v7.3;
        save(sprintf('td4debugresults%i.mat',bartrecon),'-v7.3');
    end
end

for iii=1:length(TD4in)
    for jjj=1:length(subsampin)
        plotvarinit=[plotvarinit;subpct{iii,jjj},tdtime{iii,jjj}(1),tdtime{iii,jjj}(end),scantime{iii,jjj}(1),scantime{iii,jjj}(end),...
            M0ORtis{iii,jjj}(1,2),M0ORtis{iii,jjj}(end,2),T1ORtis{iii,jjj}(1,2),T1ORtis{iii,jjj}(end,2),T2ORtis{iii,jjj}(1,2),T2ORtis{iii,jjj}(end,2)];
        plotmeangm=[plotmeangm;meancell{iii,jjj}(1,1),meancell{iii,jjj}(1,1),meancell{iii,jjj}(2,1),meancell{iii,jjj}(2,1),meancell{iii,jjj}(3,1),meancell{iii,jjj}(3,1)];
        plotmeanwm=[plotmeanwm;meancell{iii,jjj}(1,2),meancell{iii,jjj}(1,2),meancell{iii,jjj}(2,2),meancell{iii,jjj}(2,2),meancell{iii,jjj}(3,2),meancell{iii,jjj}(3,2)];
        plotmedgm=[plotmedgm;mediancell{iii,jjj}(1,1),mediancell{iii,jjj}(1,1),mediancell{iii,jjj}(2,1),mediancell{iii,jjj}(2,1),mediancell{iii,jjj}(3,1),mediancell{iii,jjj}(3,1)];
        plotmedwm=[plotmedwm;mediancell{iii,jjj}(1,2),mediancell{iii,jjj}(1,2),mediancell{iii,jjj}(2,2),mediancell{iii,jjj}(2,2),mediancell{iii,jjj}(3,2),mediancell{iii,jjj}(3,2)];
        plotvargm=[plotvargm;varcell{iii,jjj}(1,1),varcell{iii,jjj}(1,1),varcell{iii,jjj}(2,1),varcell{iii,jjj}(2,1),varcell{iii,jjj}(3,1),varcell{iii,jjj}(3,1)];
        plotvarwm=[plotvarwm;varcell{iii,jjj}(1,2),varcell{iii,jjj}(1,2),varcell{iii,jjj}(2,2),varcell{iii,jjj}(2,2),varcell{iii,jjj}(3,2),varcell{iii,jjj}(3,2)];
        %         plotvaropt=[plotvaropt;scantime{iii,jjj}(end,2),M0ORtis{iii,jjj}(end,2),T1ORtis{iii,jjj}(end,2),T2ORtis{iii,jjj}(end,2)];
    end
end
% save td4debugresults.mat -v7.3;
save(sprintf('td4debugresults%i.mat',bartrecon),'-v7.3');

%% Plot Figures
plotsym={'-o','--x',':+','-.s','-*','--d',':v','-.^','-<','-->',':p','-.h'};
% plotsym={'b-o','b--x','b:x','b-.s','b-*','r-o','r--x','r:x','r-.s','r-*'};
szgs=size(M0ORtis);

titlename='GM M0 Recon Mean and Median';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,2),plotmeangm(iii:szgs(2):end,1),['r',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotmedgm(iii:szgs(2):end,1),['b',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Mean or Median'); title(titlename);
legendkey={'100% Mean','100% Median','70% Mean','70% Median','50% Mean','50% Median','25% Mean','25% Median'};
legend(legendkey,'Location','NorthEast');

titlename='GM T1 Recon Mean and Median';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,2),plotmeangm(iii:szgs(2):end,3),['r',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotmedgm(iii:szgs(2):end,3),['b',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Mean or Median'); title(titlename);
legendkey={'100% Mean','100% Median','70% Mean','70% Median','50% Mean','50% Median','25% Mean','25% Median'};
legend(legendkey,'Location','NorthEast');

titlename='GM T2 Recon Mean and Median';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,2),plotmeangm(iii:szgs(2):end,5),['r',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotmedgm(iii:szgs(2):end,5),['b',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Mean or Median'); title(titlename);
legendkey={'100% Mean','100% Median','70% Mean','70% Median','50% Mean','50% Median','25% Mean','25% Median'};
legend(legendkey,'Location','NorthEast');

titlename='WM M0 Recon Mean and Median';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,2),plotmeanwm(iii:szgs(2):end,1),['r',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotmedwm(iii:szgs(2):end,1),['b',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Mean or Median'); title(titlename);
legendkey={'100% Mean','100% Median','70% Mean','70% Median','50% Mean','50% Median','25% Mean','25% Median'};
legend(legendkey,'Location','NorthEast');

titlename='WM T1 Recon Mean and Median';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,2),plotmeanwm(iii:szgs(2):end,3),['r',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotmedwm(iii:szgs(2):end,3),['b',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Mean or Median'); title(titlename);
legendkey={'100% Mean','100% Median','70% Mean','70% Median','50% Mean','50% Median','25% Mean','25% Median'};
legend(legendkey,'Location','NorthEast');

titlename='WM T2 Recon Mean and Median';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,2),plotmeanwm(iii:szgs(2):end,5),['r',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotmedwm(iii:szgs(2):end,5),['b',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Mean or Median'); title(titlename);
legendkey={'100% Mean','100% Median','70% Mean','70% Median','50% Mean','50% Median','25% Mean','25% Median'};
legend(legendkey,'Location','NorthEast');

titlename='GM M0 Recon Var';
figure; hold on;
for iii=1:szgs(2)
    %     plot(plotvar1(iii:szgs(2):end,3),plotvar1(iii:szgs(2):end,7),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvargm(iii:szgs(2):end,1),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Variance'); title(titlename);
legendkey={'100% Init','70% Init','50% Init','25% Init'};
legend(legendkey{1:szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='GM T1 Recon Var';
figure; hold on;
for iii=1:szgs(2)
    %     plot(plotvar1(iii:szgs(2):end,3),plotvar1(iii:szgs(2):end,9),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvargm(iii:szgs(2):end,3),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Variance'); title(titlename);
legendkey={'100% Init','70% Init','50% Init','25% Init'};
legend(legendkey{1:szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='GM T2 Recon Var';
figure; hold on;
for iii=1:szgs(2)
    %     plot(plotvar1(iii:szgs(2):end,3),plotvar1(iii:szgs(2):end,11),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvargm(iii:szgs(2):end,5),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Variance'); title(titlename);
legendkey={'100% Init','70% Init','50% Init','25% Init'};
legend(legendkey{1:szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='WM M0 Recon Var';
figure; hold on;
for iii=1:szgs(2)
    %     plot(plotvar2(iii:szgs(2):end,3),plotvar2(iii:szgs(2):end,7),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarwm(iii:szgs(2):end,1),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Variance'); title(titlename);
legendkey={'100% Init','70% Init','50% Init','25% Init'};
legend(legendkey{1:szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='WM T1 Recon Var';
figure; hold on;
for iii=1:szgs(2)
    %     plot(plotvar2(iii:szgs(2):end,3),plotvar2(iii:szgs(2):end,9),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarwm(iii:szgs(2):end,3),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Variance'); title(titlename);
legendkey={'100% Init','70% Init','50% Init','25% Init'};
legend(legendkey{1:szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='WM T2 Recon Var';
figure; hold on;
for iii=1:szgs(2)
    %     plot(plotvar2(iii:szgs(2):end,3),plotvar2(iii:szgs(2):end,11),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarwm(iii:szgs(2):end,5),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Variance'); title(titlename);
legendkey={'100% Init','70% Init','50% Init','25% Init'};
legend(legendkey{1:szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');



titlename='White Matter M0 Reconstruction';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,3),plotvarinit(iii:szgs(2):end,7),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarinit(iii:szgs(2):end,6),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter M0 Reconstruction - Acq Time';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,5),plotvarinit(iii:szgs(2):end,7),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,4),plotvarinit(iii:szgs(2):end,6),['r',plotsym{iii}]);
end
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter T1 Reconstruction';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,3),plotvarinit(iii:szgs(2):end,9),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarinit(iii:szgs(2):end,8),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter T1 Reconstruction - Acq Time';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,5),plotvarinit(iii:szgs(2):end,9),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,4),plotvarinit(iii:szgs(2):end,8),['r',plotsym{iii}]);
end
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter T2 Reconstruction';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,3),plotvarinit(iii:szgs(2):end,11),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,2),plotvarinit(iii:szgs(2):end,10),['r',plotsym{iii}]);
end
xlabel('TD Sum (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

titlename='White Matter T2 Reconstruction - Acq Time';
figure; hold on;
for iii=1:szgs(2)
    plot(plotvarinit(iii:szgs(2):end,5),plotvarinit(iii:szgs(2):end,11),['b',plotsym{iii}]);
    plot(plotvarinit(iii:szgs(2):end,4),plotvarinit(iii:szgs(2):end,10),['r',plotsym{iii}]);
end
xlabel('Acquisition Time (s)'); ylabel('Overall Error'); title(titlename);
legendkey={'100% Opt','100% Init','70% Opt','70% Init','50% Opt','50% Init','25% Opt','25% Init'};
legend(legendkey{1:2*szgs(2)},'Location','NorthEast');
% saveas(gcf,[resultspathname,titlename],'png');

crashplots=0;
if crashplots==1
    for jjj=1:length(subsampin)
        for iii=1:length(TD4in)
            clear gmmmplot wmmmplot;
            TD=[0.5,0.001,0.1,TD4in(iii)];
            dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
            dtplot=cumsum(dt)';
            mmvec=reshape(Mmeassave{iii,jjj},[size(Mmeassave{iii,jjj},1)*size(Mmeassave{iii,jjj},2)*size(Mmeassave{iii,jjj},3),size(Mmeassave{iii,jjj},4)]);
            mmvecfull=reshape(Mmeasfullsave{iii,jjj},[size(Mmeasfullsave{iii,jjj},1)*size(Mmeasfullsave{iii,jjj},2)*size(Mmeasfullsave{iii,jjj},3),size(Mmeasfullsave{iii,jjj},4)]);
            for yyy=1:size(Mmeassave{iii,jjj},4)
                gmmmplot(:,yyy)=mmvec(materialID(:)==1,yyy);
                wmmmplot(:,yyy)=mmvec(materialID(:)==2,yyy);
            end
            
            titlename=sprintf('Gray Matter, TD4=%2.2f, SS=%2.2f',TD4in(iii),subsampin(jjj));
            figure; hold on;
            for yyy=1:size(gmmmplot,1)
                plot(dtplot([2,6:2:end-1]),squeeze(gmmmplot(yyy,:)),'-o');
            end
            plot(dtplot,squeeze(sind(flipAngle)*Mmeasfullsave{iii,jjj}(16,58,1,:)),'r--s');
            xlabel('Time (s)'); ylabel('M'); title(titlename);
            
            titlename=sprintf('White Matter, TD4=%2.2f, SS=%2.2f',TD4in(iii),subsampin(jjj));
            figure; hold on;
            for yyy=1:size(wmmmplot,1)
                plot(dtplot([2,6:2:end-1]),squeeze(wmmmplot(yyy,:)),'-o');
            end
            plot(dtplot,squeeze(sind(flipAngle)*Mmeasfullsave{iii,jjj}(49,36,1,:)),'r--s');
            xlabel('Time (s)'); ylabel('M'); title(titlename);
        end
    end
end

end