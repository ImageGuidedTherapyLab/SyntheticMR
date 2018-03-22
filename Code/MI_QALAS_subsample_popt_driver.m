
%% MI_QALAS_popt_driver

filename='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/MI_QALAS_subsample_poptrecons.mat';
if exist(filename,'file')==2
    delete(filename);
end

for ttotal=[5,10]
    for optcase=4:7;
        %     try
        
        close all; clearvars -except optcase ttotal;
        
        %% Optimization Space Acquisition Parameters
        % optcase=1;
%         optcaseacc=[1,2,3,4,100];
        optcaseacc=[1,1,1,1];
        optcasevar=[1,3,10,20];
        geometrycase=2;
        lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
        switch geometrycase
            case 1
                tmpmatid = int32(1);
            case 1.5
                tmpmatid = int32([1,2,3]);
            case 1.75
                system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
                tmptissue = load_untouch_nii('resampleimg.nii.gz');
                system('rm resampleimg.nii.gz');
                tmpmatid=int32(tmptissue.img(:,:,ceil(size(tmptissue.img,3)/2)));
            case 2
                tmptissue = load_untouch_nii(lfname);
                tmpmatid=int32(tmptissue.img(15:165,20:200,92));
            case 2.5
                system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
                tmptissue = load_untouch_nii('resampleimg.nii.gz');
                system('rm resampleimg.nii.gz');
                tmpmatid = int32(tmptissue.img);
            case 0
                system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
                tmptissue = load_untouch_nii('resampleimg.nii.gz');
                system('rm resampleimg.nii.gz');
                tmpmatid = int32(tmptissue.img);
            otherwise
                tmptissue = load_untouch_nii(lfname);
                tmpmatid = int32(tmptissue.img);
        end
        
        %% Default Acquisition Parameters
        flipAngle = 4;           % deg
        TR = 0.005;              % s
        TE_T2prep = 0.100;       % s
        Tacq = 0.500;            % s
        TDpT2 = 0.03;            % s
        TDinv = 0;               % s
        nacq = 5;
        TD = [1,1,1,1];          % s
        
        acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];
%             dt=[0,0,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];

        subsmplconstrain=bart(sprintf('poisson -Y %i -Z %i -y %f -z %f -V %f',size(tmpmatid,1),size(tmpmatid,2),optcaseacc(optcase-3),optcaseacc(optcase-3),optcasevar(optcase-3)));
        tconstrain=ttotal*60/ceil(sum(subsmplconstrain(:))/100)-TE_T2prep-TDpT2-nacq*Tacq-TDinv; % seconds
        
        B1inhomflag=1;
        pspacelabels={'flipAngle','TD(1)','TD(2)','TD(3)'};%,'TD(4)'};
        %         pspacelabels={'flipAngle'};
        subsmpllabels={};
        %         pspacelabels={};
        %         subsmpllabels={'variance'};
        pinit=[4,1,1,1]';
        pAeq=zeros(length(pinit));
        pAeq(2,2:end)=[1,1,1];
        pbeq=zeros(size(pinit));
        pbeq(2)=tconstrain;
        pmin=[0,0,0,0]';
        pmax=[180,2,2,2]';
        %         pinit=[10];
        %         pmin=[0];
        %         pmax=[200];
        findiffrelstep=1.e-5;
        tolx=1.e-5;
        tolfun=1.e-5;
        maxiter=500;
        
        %% Driver Function
        %         tic;
        [popt,fval,exitflag,output,lambda,grad,hessian] = MI_QALAS_subsample_popt...
            (optcase,geometrycase,B1inhomflag,pspacelabels,subsmpllabels,acqparam,pinit,pAeq,pbeq,pmin,pmax,findiffrelstep,tolx,tolfun,maxiter);
        %         toc;
        
        %% Save
        figure(1)
        saveas(gcf,sprintf('Figures/MIopt_subsamp_%f_%f.png',ttotal,optcase));
        figure(2)
        saveas(gcf,sprintf('Figures/Paramopt_subsamp_%f_%f.png',ttotal,optcase));
        save(sprintf('results/optresults_subsamp_%f_%f.mat',ttotal,optcase),'-v7.3');
        system(sprintf('mv opt_history.txt results/opt_history_%f_%f.txt',ttotal,optcase));
        system(sprintf('mv /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/MI_QALAS_subsample_poptrecons.mat /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f.mat',ttotal,optcase));
        %     catch
        %     end
        
    end
end