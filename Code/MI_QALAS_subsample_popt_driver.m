
function [] = MI_QALAS_subsample_popt_driver(tconoverride,acqtimes,pdvval,runnumber)

%% MI_QALAS_popt_driver

filename='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/MI_QALAS_subsample_poptrecons.mat';
if exist(filename,'file')==2
    delete(filename);
end

for ttotal=acqtimes%5:10;%[5,10]
    for pdv=pdvval%[20,25];%4:7;
        
        close all; clearvars -except pdv tconoverride ttotal acqtimes pdvval runnumber;
        pdxaccel=1;
        pdyaccel=1;
        pdarg=[pdxaccel,pdyaccel,pdv];
        
        %% Optimization Space Acquisition Parameters
        geometrycase=2;
        lfname = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
        switch geometrycase
            case 1
                materialID = int32(1);
            case 1.5
                materialID = int32([1,2,3]);
            case 1.75
                system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
                tmptissue = load_untouch_nii('resampleimg.nii.gz');
                system('rm resampleimg.nii.gz');
                materialID=int32(tmptissue.img(:,:,ceil(size(tmptissue.img,3)/2)));
            case 2
                tmptissue = load_untouch_nii(lfname);
                materialID=int32(tmptissue.img(15:165,20:200,92));
            case 2.5
                system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',lfname));
                tmptissue = load_untouch_nii('resampleimg.nii.gz');
                system('rm resampleimg.nii.gz');
                materialID = int32(tmptissue.img);
            case 0
                system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 1x1x3vox -interpolation NearestNeighbor -resample 256x192x64 -o resampleimg.nii.gz',lfname));
                tmptissue = load_untouch_nii('resampleimg.nii.gz');
                system('rm resampleimg.nii.gz');
                materialID = int32(tmptissue.img);
            otherwise
                tmptissue = load_untouch_nii(lfname);
                materialID = int32(tmptissue.img);
        end 
        
        %% Tissue Properties
% M0/T1/T2 Variance Flags
M0varflag = 1;
T1varflag = 1;
T2varflag = 1;
% CSFvarflag

% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
T1mean = [1200,  900, 4000, 1200]./1000; % s
if T1varflag~=0
    T1stdd = [ 100,  100,  200,  150]./1000; % s
else
    T1stdd = [   0,    0,    0,    0];
end
T2mean = [ 100,   80, 1000,  110]./1000; % s
if T2varflag~=0
    T2stdd = [   5,    4,   50,   10]./1000; % s
else
    T2stdd = [   0,    0,    0,    0];
end
M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
if M0varflag~=0
    M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity
else
    M0stdd = [   0,    0,    0,    0];
end

tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];
        
            %% Create synthetic data
        stdmapT1=normrnd(0,1,size(materialID));
        stdmapT2=normrnd(0,1,size(materialID));
        stdmapM0=normrnd(0,1,size(materialID));
        
        synthdataT1=(materialID==1).*(tisinput(3,1)+tisinput(4,1)*stdmapT1)+(materialID==2).*(tisinput(3,2)+tisinput(4,2)*tisinput(4,2)*stdmapT1)+(materialID==3).*(tisinput(3,3)+tisinput(4,3)*stdmapT1);
        synthdataT2=(materialID==1).*(tisinput(5,1)+tisinput(6,1)*stdmapT2)+(materialID==2).*(tisinput(5,2)+tisinput(6,2)*tisinput(6,2)*stdmapT2)+(materialID==3).*(tisinput(5,3)+tisinput(6,3)*stdmapT2);
        synthdataM0=(materialID==1).*(tisinput(1,1)+tisinput(2,1)*stdmapM0)+(materialID==2).*(tisinput(1,2)+tisinput(2,2)*tisinput(2,2)*stdmapM0)+(materialID==3).*(tisinput(1,3)+tisinput(2,3)*stdmapM0);
        synthdataT1(synthdataT1==0)=nan; synthdataT2(synthdataT2==0)=nan; synthdataM0(synthdataM0==0)=nan;
 
        goldstandardT1=(materialID==1).*tisinput(3,1)+(materialID==2).*tisinput(3,2)+(materialID==3).*tisinput(3,3);
        goldstandardT2=(materialID==1).*tisinput(5,1)+(materialID==2).*tisinput(5,2)+(materialID==3).*tisinput(5,3);
        goldstandardM0=(materialID==1).*tisinput(1,1)+(materialID==2).*tisinput(1,2)+(materialID==3).*tisinput(1,3);
        goldstandardT1(goldstandardT1==0)=nan; goldstandardT2(goldstandardT2==0)=nan; goldstandardM0(goldstandardM0==0)=nan;
        
        %% Default Acquisition Parameters
        flipAngle = 4;           % deg
        TR = 0.005;              % s
        TE_T2prep = 0.100;       % s
        Tacq = 0.500;            % s
        TDpT2 = 0.4;             % s
        TDinv = 0.03;            % s
        nacq = 5;
        TD = [1,1,1,1];          % s
        
        acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];
        
        if tconoverride==1
            tconstrain=ttotal;
        else
            if pdarg(3)==-1
                subsmplconstrain=ones([size(materialID,1),size(materialID,2)]);
            else
                subsmplconstrain=bart(sprintf('poisson -Y %i -Z %i -y %f -z %f -V %f',size(materialID,1),size(materialID,2),pdarg(1),pdarg(2),pdarg(3)));
            end
            tconstrain=ttotal*60/ceil(sum(subsmplconstrain(:))/100)-TE_T2prep-TDpT2-nacq*Tacq-TDinv; % seconds
        end
        % 80% elliptical sampling?
        
        B1inhomflag=1;
        pspacelabels={'flipAngle','TD(1)','TD(2)','TD(3)','TD(4)'};
        subsmpllabels={};   %{'variance'};
        pinit=[4,1,1,1,1]';
        pAeq=zeros(length(pinit));
        pAeq(2,2:end)=[1,1,1,1];
        pbeq=zeros(size(pinit));
        pbeq(2)=tconstrain;
        pmin=[0,0,0,0,0]';
        pmax=[180,2,2,2,2]';
        findiffrelstep=1.e-6;
        tolx=1.e-5;
        tolfun=1.e-5;
        maxiter=500;
        
        filenametmp=sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/MI_QALAS_subsample_poptreconstmp%i.mat',runnumber);
        if exist(filenametmp,'file')==2
            delete(filenametmp);
        end
        
        %% Driver Function
        [popt,fval,exitflag,output,lambda,grad,hessian] = MI_QALAS_subsample_popt...
            (pdarg,tisinput,materialID,synthdataT1,synthdataT2,synthdataM0,B1inhomflag,pspacelabels,subsmpllabels,acqparam,pinit,pAeq,pbeq,pmin,pmax,findiffrelstep,tolx,tolfun,maxiter,filenametmp);
        
        %% Save
        figure(1)
        saveas(gcf,sprintf('Figures/MIopt_subsamp_%f_%f_%f.png',tconoverride,ttotal,pdarg(3)));
        figure(2)
        saveas(gcf,sprintf('Figures/Paramopt_subsamp_%f_%f_%f.png',tconoverride,ttotal,pdarg(3)));
        save(sprintf('results/optresults_subsamp_%f_%f_%f.mat',tconoverride,ttotal,pdarg(3)),'-v7.3');
        try
        save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_goldstandards_%f_%f_%f.mat',tconoverride,ttotal,pdarg(3)),'synthdataM0','synthdataT1','synthdataT2','goldstandardM0','goldstandardT1','goldstandardT2','-v7.3');
        catch
        end
        system(sprintf('mv opt_history.txt results/opt_history_%f_%f_%f.txt',tconoverride,ttotal,pdarg(3)));
        system(sprintf('mv %s /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/results/MI_QALAS_subsample_poptrecons_%f_%f_%f.mat',filenametmp,tconoverride,ttotal,pdarg(3)));
        
    end
end