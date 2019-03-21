

T1=.6:.01:1.4;
TR=.1:.005:.5;

for iii=1:length(T1)
    for jjj=1:length(TR)
        
        tisinput=[T1(iii),.1,.1,.01];
        acqparam=[1,1,TR(jjj),.1];
        
        tic;
        [popt(iii,jjj)]=fmincon(@(x) MI_objfun_ernst(x,tisinput,acqparam),...
            45,[],[],[],[],0,90,[],...
            optimset('Display','iter-detailed'));
        %optimset('FinDiffRelStep',findiffrelstep,'TolX',tolx,'TolFun',tolfun,'MaxIter',maxiter,'Display','iter-detailed','OutputFcn',@outfun));
        toc;
        
        ernst_angle(iii,jjj)=acosd(exp(-TR(jjj)/T1(iii)));
        
    end
end

flipAngle=1:90;
for iii=1:90
    [MIobjfun(iii)]=MI_objfun_ernst(flipAngle(iii),[1,.1,.1,.01],[1,1,.5,.1]);
end
for iii=1:90
    [svalue(iii)]=spoiledgre(1,1,flipAngle(iii),.5,.1,1,.1);
end
for iii=1:90
    [svalue2(iii)]=spoiledgre(1,1,flipAngle(iii),.5,.1,4,.1);
end


save ernst_angle_globalsearch_results.mat T1 TR popt ernst_angle flipAngle MIobjfun svalue svalue2 -v7.3;

%% Plots

load ernst_angle_globalsearch_results.mat;

xratio=0:.1:5;
figure;plot(xratio,acosd(exp(-xratio)));
xticks(0:.5:5);
xlabel('TR/T1'); ylabel('Ernst Angle (^\circ)'); title('Theoretical Ernst Angle');
saveas(gcf,'Figures/ernstangleargs','png');

figure; plot(flipAngle,svalue);
hold on; plot(flipAngle,svalue2);
xlabel('Flip Angle (^\circ)'); ylabel('Signal Intensity'); title('Ernst Angle Optimization');
legend('White matter','Cerebrospinal fluid','location','east');
saveas(gcf,'Figures/ernstangle2tissue','png');

figure;plot(flipAngle,svalue/.05);
xlabel('Flip Angle (^\circ)'); ylabel('SNR'); title('Flip Angle Optimization');
yyaxis right;
plot(flipAngle,-MIobjfun);
ylabel('Mutual Information');
legend('Signal Model SNR','Information Model MI','location','east');
saveas(gcf,'Figures/mi_globsrch_ernstangles','png');

figure; contourf(TR,T1,popt,15); cbar=colorbar;
ylabel(cbar,'Flip Angle (^\circ)');
xlabel('TR (s)'); ylabel('T1 (s)'); title('Mutual Information-Optimized Flip Angle');
saveas(gcf,'Figures/mi_ernstangles','png');

figure; contourf(TR,T1,ernst_angle,15); cbar=colorbar;
ylabel(cbar,'Flip Angle (^\circ)');
xlabel('TR (s)'); ylabel('T1 (s)'); title('Theoretical Ernst Angle');
saveas(gcf,'Figures/thr_ernstangles','png');

figure; contourf(TR,T1,100*(popt-ernst_angle)./ernst_angle,15); cbar=colorbar;
ylabel(cbar,'Relative Error (%)');
xlabel('TR (s)'); ylabel('T1 (s)'); title('Relative Error');
saveas(gcf,'Figures/err_ernstangles','png');
