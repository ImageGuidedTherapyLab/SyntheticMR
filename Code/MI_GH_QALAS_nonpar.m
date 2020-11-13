%% Measurement model
% clear all
% close all
% format shortg

function [kspace,MI,Hz]=MI_GH_QALAS_nonpar()

%% Input files
labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation

%% signal model parameters
%%          GM       WM        CSF     Tumor
T1mean = [1400   , 1000    ,  4000    , 250  ]./1000; % s
T1stdd = [100    ,  100    ,   200    ,  50  ]./1000; % s
T2mean = [100    ,   75    ,   600    , 250  ]./1000; % s
T2stdd = [  5    ,    5    ,    30    ,  50  ]./1000; % s
M0mean = [  0.9  ,    0.8  ,     1.0  ,   1.2];       % relative intensity
M0stdd = [   .05 ,     .05 ,      .05 ,    .4];       % relative intensity
TR = 0.0026;              % s
TE_T2prep = 0.100;        % s
flipAngle = 5;            % degrees
T2TE = 0.100;             % s
T1TE = [0:50:1000]./1000; % s
%[0:200:10000]./1000;
T1TD = [0:10:200]./1000;  % s
%[0:7:140]./1000;
nacq=2;
% nacq = length(dt)/2-6;

% repetitiontime  = 4000; % ms
% delaytime       = [0:200:10000 ] ; % ms
% echotime        = [0:7:140 ] ; % ms

%% Quadrature Points
NumGP=4;
[x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumGP,[mean(T1mean(1:2))],[mean(T1stdd(1:2))]);
[x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumGP,[mean(T2mean(1:2))],[mean(T2stdd(1:2))]);
% [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumGP,[mean(T1mean(1:2)),T1mean(3)],[mean(T1stdd(1:2)),T1stdd(3)]);
% [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumGP,[mean(T2mean(1:2)),T2mean(3)],[mean(T2stdd(1:2)),T2stdd(3)]);

%% Loading tissue types
disp('loading tissue types');

try
    system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
    tissuelabel = load_untouch_nii('resampleimg.nii.gz');
    system('rm resampleimg.nii.gz');
catch
    tissuelabel = load_untouch_nii(labelfilename);
end
% tissuelabel  = load_untouch_nii(labelfilename);

%% TODO - downsample materialID by factor of 2
% materialID = imresize3(int32(tissuelabel.img),0.25);

materialID = int32(tissuelabel.img);
% materialID=[1,2,3];
% materialID=1;

%% loop over delay time and echo time
lt1te=length(T1TE);
lt1td=length(T1TD);
[x1,x2]=ndgrid(xm_t1,xm_t2);
x1=x1(:);
x2=x2(:);

M0map = M0mean(1) * (materialID == 1) + M0mean(2) * (materialID == 2) + M0mean(3) * (materialID == 3) ;
for gp2 = 1:NumGP
    for gp1 = 1:NumGP
        % Generate parameters maps from mean
        T1map = eps * (materialID==0) + xm_t1(gp1,1) * or(materialID == 1 , materialID == 2) + xm_t1(gp1,1) * (materialID == 3) ;
        T2map = eps * (materialID==0) + xm_t2(gp2,1) * or(materialID == 1 , materialID == 2) + xm_t2(gp2,1) * (materialID == 3) ;
        
        for jjj = 1:length(T1TD)
            for iii = 1:length(T1TE)
                dt = T2TE*[0,0,1,0,0,0,0,0]+T1TE(iii)*[0,0,1,0,0,0,1,0]+T1TD(jjj)*[0,0,0,0,0,1,0,0];
                
                disp(sprintf('Model eval: %d of %d, %d of %d, %d of %d, %d of %d',iii,length(T1TE),jjj,length(T1TD),gp1,NumGP,gp2,NumGP))
                %                 waitbar(((iii-1)*length(T1TD)*NumGP*NumGP+(jjj-1)*NumGP*NumGP+(gp1-1)*NumGP+gp2)/(length(T1TE)*length(T1TD)*NumGP*NumGP),h,...
                %                     ['Model Evaluation: ',sprintf('%d of %d, %d of %d, %d of %d, %d of %d',iii,length(T1TE),jjj,length(T1TD),gp1,NumGP,gp2,NumGP)]);
                
                [~,Mmodel]=qalas(M0map,M0map,T1map,T2map,TR,TE_T2prep,flipAngle,nacq,dt);
                
                % tranform to fourier space. note this is a function of TD and TE
                for mmm=1:size(Mmodel,4)
                    kspace(:,:,:,mmm,iii,jjj,gp1,gp2) = fftshift(fftn(Mmodel(:,:,:,mmm)));
                end
            end
        end
    end
end

save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas_mi_results_nq%d.mat',NumGP),'kspace','-v7.3');

ksize=size(kspace);
ksize=ksize(1:6);
ksr=real(kspace);
ksi=imag(kspace);

wn_k=kron(wn_t1,wn_t2');

Ezr=zeros(ksize);
Ezi=zeros(ksize);
Sigrr=zeros(ksize);
Sigii=zeros(ksize);
Sigri=zeros(ksize);
for jjj=1:size(wn_k,2)
    for iii=1:size(wn_k,1)
        disp(sprintf('Quadrature: %d of %d, %d of %d',iii,NumGP,jjj,NumGP))
        Ezr=Ezr+wn_k(iii,jjj).*ksr(:,:,:,:,:,:,iii,jjj);
        Ezi=Ezi+wn_k(iii,jjj).*ksi(:,:,:,:,:,:,iii,jjj);
        Sigrr=Sigrr+wn_k(iii,jjj).*ksr(:,:,:,:,:,:,iii,jjj).^2;
        Sigii=Sigii+wn_k(iii,jjj).*ksi(:,:,:,:,:,:,iii,jjj).^2;
        Sigri=Sigri+wn_k(iii,jjj).*ksr(:,:,:,:,:,:,iii,jjj).*ksi(:,:,:,:,:,:,iii,jjj);
    end
end

N=2;
signu=0.1;
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas_mi_results_nq%d.mat',NumGP),'kspace','MI','Hz','-v7.3');
