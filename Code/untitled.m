flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
plotvarinit=[]; plotmeangm=[]; plotmeanwm=[]; plotmedgm=[]; plotmedwm=[]; plotvargm=[]; plotvarwm=[];
plotsym={'-o','--x',':+','-.s','-*','--d',':v','-.^','-<','-->',':p','-.h'};
figure; hold on;
for jjj=1:length(subsampin)   %subsample pct
    plotx=[];ploty=[];ploty2=[];
    for iii=1:length(TD4in)  %TD4 times
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
        TD=[0.5,0.001,0.1,TD4in(iii)];
        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [Mmeasfull,Mmeas]=qalas(goldstandardM0,goldstandardM0,goldstandardT1,goldstandardT2,TR,TE_T2prep,flipAngle,nacq,dt);
        % stdmapmeas=normrnd(0,signu,size(materialID));
        % Mmeas=Mmeas+stdmapmeas;
        plotx=[plotx,TD4in(iii)];
        ploty=[ploty,Mmeasfull(16,58,1,1)];
        ploty2=[ploty2,Mmeasfull(49,36,1,1)];
%             Mmeasfullsave{iii,jjj}(16,58,1,:)
%             Mmeasfullsave{iii,jjj}(49,36,1,:)
    end
        figure(1); hold on;
        plot(plotx,ploty,plotsym{iii});
        figure(2); hold on;
        plot(plotx,ploty2,plotsym{iii});
end
