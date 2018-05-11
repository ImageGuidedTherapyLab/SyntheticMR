function [M0pred,T1pred,T2pred,varstats,meanstats,medianstats] = MIQALASrecondebug()

%% Reconstruct to Compute Variances
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
        
        synthdataT1=(materialID==1).*(tisinput(3,1)+tisinput(4,1)*stdmapT1)+(materialID==2).*(tisinput(3,2)+tisinput(4,2)*stdmapT1)+(materialID==3).*(tisinput(3,3)+tisinput(4,3)*stdmapT1);
        synthdataT2=(materialID==1).*(tisinput(5,1)+tisinput(6,1)*stdmapT2)+(materialID==2).*(tisinput(5,2)+tisinput(6,2)*stdmapT2)+(materialID==3).*(tisinput(5,3)+tisinput(6,3)*stdmapT2);
        synthdataM0=(materialID==1).*(tisinput(1,1)+tisinput(2,1)*stdmapM0)+(materialID==2).*(tisinput(1,2)+tisinput(2,2)*stdmapM0)+(materialID==3).*(tisinput(1,3)+tisinput(2,3)*stdmapM0);
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
        TD = [0.5,0.5,0.5,0.5];          % s
        
        acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];
        
%% Create synthetic QALAS measurements
dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
[~,Mmeas]=qalas(synthdataM0,synthdataM0,synthdataT1,synthdataT2,TR,TE_T2prep,flipAngle,nacq,dt);
stdmapmeas=normrnd(0,signu,size(materialID));
Mmeas=Mmeas+stdmapmeas;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE FOLLOWING SECTION ONLY WORKS FOR 2D!
% FIX IT FOR N-D!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subsample synthetic measurements
Mmeassub=Mmeas;
if pdarg(3)~=-1
    Mmeassub(isnan(Mmeassub))=0;
    kmeas=bart('fft 3',Mmeassub);
    subsample=squeeze(subsmplmask(1,:,:));
    subsample=repmat(subsample,[1,1,size(kmeas,3),size(kmeas,4)]);
    %         Mmeassub=bart('pics',kmeas.*subsample,ones(size(kmeas)));
    Mmeassub=bart('fft -i 3',kmeas.*subsample)*size(kmeas,4)/numel(kmeas);
    Mmeassub=real(Mmeassub);
    Mmeassub(repmat(materialID,[1,1,size(Mmeassub,3),size(Mmeassub,4)])==0)=nan;
end
%     end

%% Reconstruct synthetic QALAS measurements
% Optimization solution for M0 and T1 prediction
xinit=[mean(tisinput(1,1:3)),mean(tisinput(3,1:3))];%,mean(tisinput(5,1:3))];
smeas=size(Mmeassub);
Mmeasvec=reshape(Mmeassub,[prod(smeas(1:3)),smeas(4:end)]);
mmvsize=size(Mmeasvec,1);
parfor iii=1:size(Mmeasvec,1)
    if sum(isnan(squeeze(Mmeasvec(iii,:))))>0
        M0predvec(iii)=0;%nan;
        T1predvec(iii)=0;%nan;
        T2predvec(iii)=0;%nan;
    else
        xm=fminsearch(@(x) qalasobjfun(x,squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt),xinit);
        M0predvec(iii)=xm(1);
        T1predvec(iii)=xm(2);
        T2predvec(iii)=qalasT2calc(M0predvec(iii),T1predvec(iii),squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt);
    end
end
M0pred(:,:,:)=reshape(M0predvec,smeas(1:3));
T1pred(:,:,:)=reshape(T1predvec,smeas(1:3));
T2pred(:,:,:)=reshape(T2predvec,smeas(1:3));

M0p1set=M0pred(:).*(materialID(:)==1); M0p1set=M0p1set(M0p1set~=0);
M0p2set=M0pred(:).*(materialID(:)==2); M0p2set=M0p2set(M0p2set~=0);
M0p3set=M0pred(:).*(materialID(:)==3); M0p3set=M0p3set(M0p3set~=0);
T1p1set=T1pred(:).*(materialID(:)==1); T1p1set=T1p1set(T1p1set~=0);
T1p2set=T1pred(:).*(materialID(:)==2); T1p2set=T1p2set(T1p2set~=0);
T1p3set=T1pred(:).*(materialID(:)==3); T1p3set=T1p3set(T1p3set~=0);
T2p1set=T2pred(:).*(materialID(:)==1); T2p1set=T2p1set(T2p1set~=0);
T2p2set=T2pred(:).*(materialID(:)==2); T2p2set=T2p2set(T2p2set~=0);
T2p3set=T2pred(:).*(materialID(:)==3); T2p3set=T2p3set(T2p3set~=0);

M0varmeas=[quantile(M0p1set,0.975)-quantile(M0p1set,0.025)...
    quantile(M0p2set,0.975)-quantile(M0p2set,0.025)...
    quantile(M0p3set,0.975)-quantile(M0p3set,0.025)];
T1varmeas=[quantile(T1p1set,0.975)-quantile(T1p1set,0.025)...
    quantile(T1p2set,0.975)-quantile(T1p2set,0.025)...
    quantile(T1p3set,0.975)-quantile(T1p3set,0.025)];
T2varmeas=[quantile(T2p1set,0.975)-quantile(T2p1set,0.025)...
    quantile(T2p2set,0.975)-quantile(T2p2set,0.025)...
    quantile(T2p3set,0.975)-quantile(T2p3set,0.025)];
varstats=[M0varmeas;T1varmeas;T2varmeas];

meanstats=[mean(M0p1set),mean(M0p2set),mean(M0p3set);
    mean(T1p1set),mean(T1p2set),mean(T1p3set);
    mean(T2p1set),mean(T2p2set),mean(T2p3set)];

medianstats=[median(M0p1set),median(M0p2set),median(M0p3set);
    median(T1p1set),median(T1p2set),median(T1p3set);
    median(T2p1set),median(T2p2set),median(T2p3set)];

end