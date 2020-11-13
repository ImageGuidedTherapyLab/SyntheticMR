
%% MI_QALAS_popt_driver

for optcase=3:3
    
%     try
        
        close all; clearvars -except optcase;
        
        %% Optimization Space Acquisition Parameters
        % optcase=1;
        geometrycase=2.5;
        B1inhomflag=1;
        pspacelabels={'flipAngle','TD(1)','TD(2)','TD(3)','TD(4)'};
        pinit=[4,1,1,1,1];
        pmin=[0,0,0,0,0];
        pmax=[180,2,2,2,2];
        findiffrelstep=1.e-5;
        tolx=1.e-5;
        tolfun=1.e-5;
        maxiter=500;
        
        %% Driver Function
%         tic;
        [popt,fval,exitflag,output,lambda,grad,hessian] = MI_QALAS_popt...
            (optcase,geometrycase,B1inhomflag,pspacelabels,pinit,pmin,pmax,findiffrelstep,tolx,tolfun,maxiter);
%         toc;
        
        %% Save
        figure(1)
        saveas(gcf,sprintf('Figures/MIopt%i.png',optcase));
        figure(2)
        saveas(gcf,sprintf('Figures/Paramopt%i.png',optcase));
        save(sprintf('results/optresults%i.mat',optcase),'-v7.3');
        
%     catch
%     end
    
end