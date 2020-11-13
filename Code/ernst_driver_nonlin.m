

T1=.6:.01:1.4;
TR=.1:.005:.5;

% T1=.6:.1:1.4;
% TR=.1:.05:.5;

% T1=[.6,1,1.4]; TR=[.1,.5];

h=waitbar(0);

for iii=1:length(T1)
    for jjj=1:length(TR)
        
        waitbar(((iii-1)*length(TR)+jjj)/(length(T1)*length(TR)),h);
        
        tisinput=[T1(iii),.1,0,.01];
%         tisinput=[T1(iii),0.01,.1,0.01];
        acqparam=[1,1,TR(jjj),.1];
        signu=3.4762E-3;
    
        tic;
%         [popt(iii,jjj)]=fminsearch(@(x) MI_objfun_ernst_nonlin(x,tisinput,acqparam,signu),45);
%         [popt(iii,jjj)]=fminunc(@(x) MI_objfun_ernst_nonlin(x,tisinput,acqparam,signu),45);
        [popt(iii,jjj)]=fmincon(@(x) MI_objfun_ernst_nonlin(x,tisinput,acqparam,signu),...
            45,[],[],[],[],0,90,[],...
            optimset('FinDiffRel',.0001,'TolX',1E-20,'TolFun',1E-20,'Display','iter-detailed'));
        %optimset('FinDiffRelStep',findiffrelstep,'TolX',tolx,'TolFun',tolfun,'MaxIter',maxiter,'Display','iter-detailed','OutputFcn',@outfun));
        toc;
        
        ernst_angle(iii,jjj)=acosd(exp(-TR(jjj)/T1(iii)));
        
    end
end

flipAngle=1:90;
for iii=1:90
    [MIobjfun(iii)]=MI_objfun_ernst_nonlin(flipAngle(iii),[1,.1,.1,.01],[1,1,.5,.1],signu);
end
for iii=1:90
    [svalue(iii)]=spoiledgre(1,1,flipAngle(iii),.5,.1,1,.1);
end
for iii=1:90
    [svalue2(iii)]=spoiledgre(1,1,flipAngle(iii),.5,.1,4,.1);
end


save ernst_angle_globalsearch_results_nl.mat T1 TR popt ernst_angle flipAngle MIobjfun svalue svalue2 -v7.3;

%% Plots

load ernst_angle_globalsearch_results_nl.mat;

xratio=0:.1:5;
figure;plot(xratio,acosd(exp(-xratio)));
xticks(0:.5:5);
xlabel('TR/T1'); ylabel('Ernst Angle (^\circ)'); %title('Theoretical Ernst Angle');
saveas(gcf,'Figures/ernstangleargs_nl','png');

figure; plot(flipAngle,svalue);
hold on; plot(flipAngle,svalue2);
xlabel('Flip Angle (^\circ)'); ylabel('Signal Intensity'); %title('Ernst Angle Optimization');
legend('White matter','Cerebrospinal fluid','location','east');
saveas(gcf,'Figures/ernstangle2tissue_nl','png');

figure;plot(flipAngle,svalue/.05);
xlabel('Flip Angle (^\circ)'); ylabel('SNR'); %title('Flip Angle Optimization');
yyaxis right;
plot(flipAngle,-MIobjfun);
ylabel('Mutual Information');
legend('Signal Model SNR','Information Model MI','location','east');
saveas(gcf,'Figures/mi_globsrch_ernstangles_nl','png');

figure; contourf(TR,T1,popt,15); cbar=colorbar;
ylabel(cbar,'Flip Angle (^\circ)');
xlabel('TR (s)'); ylabel('T1 (s)'); %title('Mutual Information-Optimized Flip Angle');
saveas(gcf,'Figures/mi_ernstangles_nl','png');

figure; contourf(TR,T1,ernst_angle,15); cbar=colorbar;
ylabel(cbar,'Flip Angle (^\circ)');
xlabel('TR (s)'); ylabel('T1 (s)'); %title('Theoretical Ernst Angle');
saveas(gcf,'Figures/thr_ernstangles_nl','png');

figure; contourf(TR,T1,100*(popt-ernst_angle)./ernst_angle,15,'LineStyle','none'); cbar=colorbar;
caxis([-.5,.5]);
ylabel(cbar,'Relative Error (%)');
xlabel('TR (s)'); ylabel('T1 (s)'); %title('Relative Error');
saveas(gcf,'Figures/err_ernstangles_nl','png');

% Bland-Altman plots
figure; plot(ernst_angle(:),popt(:),'r.');
xlabel('Ernst Angle (^\circ)'); ylabel('MI-Optimized Flip Angle (^\circ)');
h=lsline; set(h(1),'color','black');
[lfit,gof]=fit(ernst_angle(:),popt(:),'poly1');
coef=coeffvalues(lfit);
text(23,60,sprintf('y=%fx+%f\nAdjusted r^2=%f\nSSE=%f\nRMSE=%f\nn=%i',coef(1),coef(2),gof.adjrsquare,gof.sse,gof.rmse,numel(ernst_angle)));
saveas(gcf,'Figures/ernst_blandalt1_nl','png');

figure; plot((ernst_angle(:)+popt(:))/2,-ernst_angle(:)+popt(:),'r.');
diffpop=popt(:)-ernst_angle(:);
axis([20,70,-.06,.06]);
hold on; plot([20,70],[mean(diffpop),mean(diffpop)],'black-');
plot([20,70],[mean(diffpop)-1.96*std(diffpop),mean(diffpop)-1.96*std(diffpop)],'black-.');
plot([20,70],[mean(diffpop)+1.96*std(diffpop),mean(diffpop)+1.96*std(diffpop)],'black-.');
xlabel('Mean Flip Angle (^\circ)'); ylabel('Flip Angle Difference (^\circ)');
text(37.5,-0.015,sprintf('d+1.96s=%f\nd=%f\nd-1.96s=%f',mean(diffpop)+1.96*std(diffpop),mean(diffpop),mean(diffpop)-1.96*std(diffpop)));
saveas(gcf,'Figures/ernst_blandalt2_nl','png');

figure; histogram(ernst_angle(:)-popt(:));
