        figure('pos',[10,10,1810,610]);
        subplot(1,3,1); hold on;
        plot(M0MI1([3:14,17:28,31:42]),nanstd(M0v_1([3:14,17:28,31:42],:),[],2),'bo');
        % yyaxis right;
        plot(M0MI2([3:14,17:28,31:42]),nanstd(M0v_2([3:14,17:28,31:42],:),[],2),'bo');
%         xlabel('Mutual Information'); ylabel('M0 Standard Deviation');
%         legend('Scan 1','Scan 2');
        
        subplot(1,3,2); hold on;
        plot(T1MI1([3:14,17:28,31:42]),nanstd(T1v_1([3:14,17:28,31:42],:),[],2),'bo');
%         xlabel('Mutual Information'); ylabel('T1 Standard Deviation');
        % yyaxis right;
        plot(T1MI2([3:14,17:28,31:42]),nanstd(T1v_2([3:14,17:28,31:42],:),[],2),'bo');
%         xlabel('Mutual Information'); ylabel('T1 Standard Deviation');
%         legend('Scan 1','Scan 2');
        
        subplot(1,3,3); hold on;
        plot(T2MI1([3:14,17:28,31:42]),nanstd(T2v_1([3:14,17:28,31:42],:),[],2),'bo');
        % yyaxis right;
        plot(T2MI2([3:14,17:28,31:42]),nanstd(T2v_2([3:14,17:28,31:42],:),[],2),'bo');
        axis([382,390,0,1]);
%         legend('Scan 1','Scan 2');
        
        saveas(gcf,'Figures/mivarresultsfull','tiff');
