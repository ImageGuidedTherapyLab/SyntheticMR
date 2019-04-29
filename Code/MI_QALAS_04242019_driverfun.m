
function [popt] = MI_QALAS_04242019_driverfun(tisinput,acqparam,materialID)

%% MI_QALAS_popt_driver
                       
%% Assign Acquisition Parameters
% Default Parameters
flipAngle=acqparam(1);
TR=acqparam(2);
TE_T2prep=acqparam(3);
Tacq=acqparam(4);
TDpT2=acqparam(5);
TDinv=acqparam(6);
nacq=acqparam(7);
TD=acqparam(8:6+nacq);
signu=acqparam(end);

% TD=xopt;

               
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