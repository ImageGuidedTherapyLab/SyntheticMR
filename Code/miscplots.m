
load results_hi.mat
load results_lo.mat

xlabels={'M0','T1 (s)','T2 (s)','M0','T1 (s)','T2 (s)','M0','T1 (s)','T2 (s)'};
ylabels={'PDF','PDF','PDF','PDF','PDF','PDF','PDF','PDF','PDF'};
titles={'Gray Matter','Gray Matter','Gray Matter','White Matter','White Matter','White Matter','Cerebrospinal Fluid','Cerebrospinal Fluid','Cerebrospinal Fluid'};

figure('pos',[10,10,1210,1210]);

for iii=1:9
    switch iii
           
        case {2,5}
            subplot(3,3,iii); hold on;
            posteriors_nofig(post_samples_lo,iii,{prior_lo{iii,1}}); posteriors_nofig(post_samples_hi,iii,{prior_hi{iii,1}});
            xlabel(xlabels{iii}); ylabel(ylabels{iii}); title(titles{iii});
            legend('Low Information Recon','High Information Recon');
            
        otherwise
            subplot(3,3,iii); hold on;
            posteriors_nofig(post_samples_hi,iii,{prior_hi{iii,1}}); posteriors_nofig(post_samples_lo,iii,{prior_lo{iii,1}});
            xlabel(xlabels{iii}); ylabel(ylabels{iii}); title(titles{iii});
            legend('Low Information Recon','High Information Recon');
    end
end

figure('pos',[10,10,1210,1210]);

for iii=1:9
    switch iii
        case {4,5}
            subplot(3,3,iii); hold on;
            posteriors_nofig(post_samples_hi,iii,{prior_hi{iii,1}}); posteriors_nofig(post_samples_lo,iii,{prior_lo{iii,1}});
            xlabel(xlabels{iii}); ylabel(ylabels{iii}); title(titles{iii});
            legend('Low Information Recon','High Information Recon');
            
        otherwise
            subplot(3,3,iii); hold on;
            posteriors_nofig(post_samples_lo,iii,{prior_lo{iii,1}}); posteriors_nofig(post_samples_hi,iii,{prior_hi{iii,1}});
            xlabel(xlabels{iii}); ylabel(ylabels{iii}); title(titles{iii});
            legend('Low Information Recon','High Information Recon');
    end
end
