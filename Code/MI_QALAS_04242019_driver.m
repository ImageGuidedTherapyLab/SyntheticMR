
% function [] = MI_QALAS_04242019_driver()

%% MI_QALAS_popt_driver
        
        %% Optimization Space Acquisition Parameters
        lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
        tmptissue = load_untouch_nii(lfname);
        materialID=int32(tmptissue.img(15:165,20:200,92));

        
        %% Tissue Properties
        % GM/WM/CSF M0/T1/T2 Values
        %              GM    WM   CSF  Tumor
        T1mean = [1200,  900, 4000, 1200]./1000; % s
            T1stdd = [ 100,  100,  200,  150]./1000; % s

        T2mean = [ 100,   80, 1000,  110]./1000; % s
            T2stdd = [   5,    4,   50,   10]./1000; % s

        M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
            M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity

        overwritecsf=1;
        if overwritecsf==1
            T1mean = [1200,  900, 900, 1200]./1000; % s
            T1stdd = [ 100,  100,  100,  150]./1000; % s
            
            T2mean = [ 100,   80, 80,  110]./1000; % s
            T2stdd = [   5,    4,   4,   10]./1000; % s
            
            M0mean = [ 0.9,  0.9,  0.9,  0.9];       % relative intensity
            M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
        end
        
        tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];
               
        %% Default Acquisition Parameters
        flipAngle = 4;           % deg
        TR = 0.005;              % s
        TE_T2prep = 0.100;       % s
        Tacq = 0.871;            % s
        TDpT2 = 0.4;             % s
        TDinv = 0.1;            % s
        nacq = 5;
        TD = [0.5,0.5,0.5,0.5];          % s
        
        acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];
               
%             pinit=[4,1,1,1,1]';
%             pinit=tconstrain/4*ones([4,1]);
            pinit=.5*ones([5,1]);
            pAeq=zeros(length(pinit));
            pAeq(2,:)=[1,1,1,1,1];
            pbeq=zeros(size(pinit));
            pbeq(2)=3;
            pmin=[0,0,0,0,0]';
            pmax=5*[1,1,1,1,1]';
            findiffrelstep=1.e-6;
            tolx=1.e-3;%1.e-5;
            tolfun=1.e-3;%1.e-5;
            maxiter=500;
                
        if exist('opt_history.txt','file')==2
            delete('opt_history.txt');
        end
%         filenametmp=sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/MI_QALAS_subsample_poptreconstmp%i.mat',runnumber);
%         if exist(filenametmp,'file')==2
%             delete(filenametmp);
%         end
        
        %% Driver Function
        %         [popt,fval,exitflag,output,lambda,grad,hessian] = MI_QALAS_subsample_popt...
        %             (pdarg,tisinput,materialID,synthdataT1,synthdataT2,synthdataM0,B1inhomflag,pspacelabels,subsmpllabels,acqparam,pinit,pAeq,pbeq,pmin,pmax,findiffrelstep,tolx,tolfun,maxiter,filenametmp);
        tic;
        [popt,fval,exitflag,output,lambda,grad,hessian]=fmincon(@(x) MI_QALAS_objfun_04242019(x,tisinput,acqparam,materialID),...
            pinit,pAeq,pbeq,[],[],pmin,pmax,[],...
            optimset('TolX',tolx,'TolFun',tolfun,'MaxIter',maxiter,'Display','iter-detailed','OutputFcn',@outfun));
%             optimset('FinDiffRelStep',findiffrelstep,'TolX',tolx,'TolFun',tolfun,'MaxIter',maxiter,'Display','iter-detailed','OutputFcn',@outfun));
        toc;
        
        %% Save
%             figure(1)
%             saveas(gcf,sprintf('Figures/MIopt_subsamp_%f_%f_%f.png',tconoverride,ttotal,pdarg(3)));
%             figure(2)
%             saveas(gcf,sprintf('Figures/Paramopt_subsamp_%f_%f_%f.png',tconoverride,ttotal,pdarg(3)));
%             save(sprintf('results/optresults_subsamp_%f_%f_%f.mat',tconoverride,ttotal,pdarg(3)),'-v7.3');
%             save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_goldstandards_%f_%f_%f.mat',tconoverride,ttotal,pdarg(3)),'synthdataM0','synthdataT1','synthdataT2','goldstandardM0','goldstandardT1','goldstandardT2','-v7.3');
% 
%             system(sprintf('mv opt_history.txt results/opt_history_%f_%f_%f.txt',tconoverride,ttotal,pdarg(3)));
%             system(sprintf('mv %s /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_%f.mat',filenametmp,tconoverride,ttotal,pdarg(3)));


% ttotal=(sum(popt)+TE_T2prep+nacq*Tacq+TDinv);%*numel(materialID)/100;