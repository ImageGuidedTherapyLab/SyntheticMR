
TR=0.001:0.0002:0.004;
flipAngle=1:1:10;
for trind=1:length(TR)
    for faind=1:length(flipAngle)
        
        trind
        faind
        % loop multinest
        [logZ, nest_samples, post_samples] = multinest_qalas(TR(trind), flipAngle(faind))
        
        [npost,xedges,yedges]=histcounts2(post_samples(:,1),post_samples(:,2),[50 50]);
        normpost=npost./sum(npost(:));
        normprior=ones(size(normpost))./numel(normpost);
        [xbin,~]=discretize(post_samples(:,1),50);
        [ybin,~]=discretize(post_samples(:,2),50);
        
        logLbinmean=zeros(50);
        for iii=1:50
            for jjj=1:50
                logLbin=post_samples(intersect(find(xbin==iii),find(ybin==jjj)),3);
                lllb=length(logLbin);
                if lllb==1
                    logLbinmean(iii,jjj)=logLbin;
                elseif lllb>1
                    logLbinsum=logLbin(1);
                    for kkk=2:lllb
                        logLbinsum=logplus(logLbinsum,logLbin(kkk));
                    end
                    logLbinmean(iii,jjj)=logLbinsum-log(lllb);
                end
            end
        end
        
        Hprior=normprior.*log(normprior);
        Hprior=-sum(Hprior(:));
        Hpost=normpost.*normprior.*logLbinmean;
        Hpost=sum(Hpost(:));
        MI(trind,faind)=Hprior+Hpost;
        
    end
end