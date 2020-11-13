function [plotvar1,plotvar2,plotvar3] = MIQALASrecondebug(tconoverridein,ttotalin,pdvvalin)

plotvar1=[]; plotvar2=[]; plotvar3=[];
% meancell=cell(length(ttotalin),length(pdvvalin),3,3,2);
% mediancell=cell(length(ttotalin),length(pdvvalin),3,3,2);
% varcell=cell(length(ttotalin),length(pdvvalin),3,3,2);
for iii=1:length(ttotalin)
    for jjj=1:length(pdvvalin)
        
        %% Load
        % Load single tissue optimization data
        load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_test.mat',tconoverridein,ttotalin(iii)));
        pspace1opt=pspacesave;
        
        % Load full optimization data
        load(sprintf('results/optresults_subsamp_%f_%f_%f.mat',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_goldstandards_%f_%f_%f.mat',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        load(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_%f.mat',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        opt_hist_table=readtable(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/opt_history_%f_%f_%f.txt',tconoverridein,ttotalin(iii),pdvvalin(jjj)));
        
        %% Tissue Properties
        load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/synthphantom_goldstandards.mat;
        
        for kkk=1:2
            %% Create synthetic QALAS measurements
            if kkk==1
                TD=squeeze(pspace1opt(:,:,1));
            else
                TD=squeeze(pspace1opt(:,:,end));
            end
            dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
            [~,Mmeas]=qalas(goldstandardM0,goldstandardM0,goldstandardT1,goldstandardT2,TR,TE_T2prep,flipAngle,nacq,dt);
            % stdmapmeas=normrnd(0,signu,size(materialID));
            % Mmeas=Mmeas+stdmapmeas;
            
            %% Subsample synthetic measurements
            Mmeassub=Mmeas;
            if pdarg(3)~=-1
                subsmplmask=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -V %i',size(materialID,1),size(materialID,2),pdarg(1),pdarg(2),pdarg(3)));
                Mmeassub(isnan(Mmeassub))=0;
                kmeas=bart('fft 3',Mmeassub);
                subsample=squeeze(subsmplmask(1,:,:));
                subsample=repmat(subsample,[1,1,size(kmeas,3),size(kmeas,4)]);
                %             Mmeassub=bart('pics',kmeas.*subsample,ones(size(kmeas)));
                Mmeassub=bart('fft -i 3',kmeas.*subsample)*size(kmeas,4)/numel(kmeas);
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
            
            for mmm=1:3
                for lll=1:3
                    meancell{iii,jjj,lll,mmm,kkk}=squeeze(meanstats(lll,mmm,:));
                    mediancell{iii,jjj,lll,mmm,kkk}=squeeze(medianstats(lll,mmm,:));
                    varcell{iii,jjj,lll,mmm,kkk}=squeeze(varstats(lll,mmm,:));
                end
            end
        end
    end
end

for iii=1:length(ttotalin)
    for jjj=1:length(pdvvalin)
        plotvar1=[plotvar1;meancell{iii,jjj,1,1,1},meancell{iii,jjj,1,1,2},meancell{iii,jjj,2,1,1},meancell{iii,jjj,2,1,2},meancell{iii,jjj,3,1,1},meancell{iii,jjj,3,1,2}];
        plotvar2=[plotvar2;mediancell{iii,jjj,1,2,1},mediancell{iii,jjj,1,2,2},mediancell{iii,jjj,2,2,1},mediancell{iii,jjj,2,2,2},mediancell{iii,jjj,3,2,1},mediancell{iii,jjj,3,2,2}];
        plotvar3=[plotvar3;varcell{iii,jjj,1,3,1},varcell{iii,jjj,1,3,2},varcell{iii,jjj,2,3,1},varcell{iii,jjj,2,3,2},varcell{iii,jjj,3,3,1},varcell{iii,jjj,3,3,2}];
%         plotvaropt=[plotvaropt;scantime{iii,jjj}(end,2),M0ORtis{iii,jjj}(end,2),T1ORtis{iii,jjj}(end,2),T2ORtis{iii,jjj}(end,2)];
    end
end



end