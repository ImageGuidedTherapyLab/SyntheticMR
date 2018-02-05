%% Measurement model
% clear all
% close all
% format shortg

function [Kspace,KspaceCenter,MI,H_z]=MagicSignalModel1p(materialID,T1stdflag,T2stdflag,NumGP)

NumGPz=NumGP;
if T1stdflag~=0; T1stdflag=1; end;
if T2stdflag~=0; T2stdflag=1; end;

%% Input files
% labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation

%% signal model parameters
%%         GM      WM      CSF     Tumor
T1mean = [1400  , 1000   ,  4000   , 250  ]; % ms
T1stdd = [100   ,  100   ,  200    ,  50  ]*T1stdflag; % ms
% T1stdd = [0,0,0,0];
T2mean = [100   , 75   , 600   , 250  ]; % ms
T2stdd = [5    ,  5   ,  30   ,  50  ]*T2stdflag; % ms
% T2stdd=[0,0,0,0];
M0mean = [0.9   , 0.8   , 1.0   , 1.2  ]; % relative intensity
M0stdd = [.05    ,  .05   ,  .05   ,  .4  ]; % relative intensity
noisestdd = 0.01;
flipangle = 120/180.*pi; % radian
excitationangle =  90/180.*pi; % radian
repetitiontime  = 4000; % ms
delaytime       = [0:200:4000 ] ; % ms
echotime        = [0:7:140 ] ; % ms

%% Quadrature Points
% NumGP=4;
% [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumGP,[mean(T1mean(1:2)),T1mean(3)],[mean(T1stdd(1:2)),T1stdd(3)]);
% [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumGP,[mean(T2mean(1:2)),T2mean(3)],[mean(T2stdd(1:2)),T2stdd(3)]);
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumGP,[mean(T1mean(1:2)),T1mean(3),mean(T2mean(1:2)),T2mean(3)],[mean(T1stdd(1:2)),T1stdd(3),mean(T2stdd(1:2)),T2stdd(3)]);
[xz,xnz,xmz,wz,wnz]=GaussHermiteNDGauss(NumGPz,[0,0],[noisestdd,noisestdd]);

%% Loading tissue types
disp('loading tissue types');

% system(sprintf('c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
% tissuelabel = load_untouch_nii('resampleimg.nii.gz');
% system('rm resampleimg.nii.gz');

% tissuelabel  = load_untouch_nii(labelfilename);

%% TODO - downsample materialID by factor of 2
% materialID = imresize3(int32(tissuelabel.img),0.25);

% materialID = int32(tissuelabel.img);
% materialID=[1,2,3];
% materialID=1;

%% loop over delay time and echo time
for iii = 1:length(delaytime)
    for jjj = 1:length(echotime)
        for gp1 = 1:NumGP
            for gp2 = 1:NumGP
                for gp3 = 1:NumGP
                    for gp4 = 1:NumGP
                        %% Generate parameters maps from mean
                        %% FIXME - @dmitchell412 need quadrature points here
                        %                         T1map = eps * (materialID==0) + T1gmwm(gp1) * or(materialID == 1 , materialID == 2) + T1mean(gp2) * (materialID == 3) ;
                        T1map = eps * (materialID==0) + xm(gp1,1) * or(materialID == 1 , materialID == 2) + xm(gp2,2) * (materialID == 3) ;
                        %                         T2map = eps * (materialID==0) + T2gmwm(gp3) * or(materialID == 1 , materialID == 2) + T2mean(gp4) * (materialID == 3) ;
                        T2map = eps * (materialID==0) + xm(gp3,3) * or(materialID == 1 , materialID == 2) + xm(gp4,4) * (materialID == 3) ;
                        M0map = M0mean(1) * (materialID == 1) + M0mean(2) * (materialID == 2) + M0mean(3) * (materialID == 3) ;
                        
                        disp( sprintf( '%d of %d, %d of %d',iii ,length(delaytime ) ,jjj,length(echotime        )))
                        eTD = exp(-delaytime(iii)./T1map );
                        eTR = exp(-repetitiontime./T1map );
                        eTE = exp(-echotime(jjj)./T2map );
                        % compute magnetization image in image domain
                        Mmodel = M0map.*( 1 - (1 - cos(flipangle)) * eTD  - cos(flipangle) * eTR ).*((1- cos(flipangle)*eTR*cos(excitationangle) ).^(-1) ).*eTE ;
                        % tranform to fourier space. note this is a function of TD and TE
                        Kspace{iii,jjj,gp1,gp2,gp3,gp4} = fftshift(fftn(Mmodel));
                        ksc=floor(size(Kspace)/2+1);
%                         KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4) = Kspace(ksc(1),ksc(2),ksc(3));
                        KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4)=Kspace{iii,jjj,gp1,gp2,gp3,gp4}(1);
                    end
                end
            end
        end
    end
end

% save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/magicsignalmodelresults_nq%d.mat',NumGP),'KspaceCenter','-v7.3');
save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/magicsignalmodel1p_results_%d%d%d%d.mat',materialID,T1stdflag,T2stdflag,NumGP),'Kspace','KspaceCenter','-v7.3');

a=(2*pi*noisestdd)^(-1);
c=2*noisestdd^2;
for iii = 1:length(delaytime)
    iii
    for jjj = 1:length(echotime)
        Hz_i=0;
        for gp1 = 1:NumGP
            for gp2 = 1:NumGP
                for gp3 = 1:NumGP
                    for gp4 = 1:NumGP
                        Hz_j=0;
                        for jr=1:NumGPz
                            for ji=1:NumGPz
                                b_r=xmz(jr,1)*ones(size(wn))+real(squeeze(KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4)));
                                b_i=xmz(ji,2)*ones(size(wn))+imag(squeeze(KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4)));
                                Hz_ktmp=0; Hz_k=0;
                                Hz_ktmp=a*wn.*exp(-((b_r-real(squeeze(KspaceCenter(iii,jjj,:,:,:,:)))).^2+(b_i-imag(squeeze(KspaceCenter(iii,jjj,:,:,:,:)))).^2)/c);
                                Hz_k=sum(Hz_ktmp(:));
                                Hz_j=Hz_j+wnz(jr,ji)*log(Hz_k);
                            end
                        end
                        Hz_i=Hz_i+wn(gp1,gp2,gp3,gp4)*Hz_j;
                    end
                end
            end
        end
        H_z(iii,jjj)=-Hz_i;
        MI(iii,jjj)=H_z(iii,jjj)-0.5*log((2*pi*exp(1))^2*noisestdd^4);
    end
end

save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/magicsignalmodel1p_results_%d%d%d%d.mat',materialID,T1stdflag,T2stdflag,NumGP),'Kspace','KspaceCenter','MI','H_z','-v7.3');

% figure;imagesc(abs(KspaceCenter(end:-1:1,:,1,1,1,1)));colorbar;
figure;imagesc(MI(end:-1:1,:));colorbar;
set(gca,'XTickLabel',[28:35:140],'YTickLabel',[3600:-400:0]);
xlabel('Echo Time (ms)'); ylabel('Delay Time (ms)');
title('Mutual Information between Signal Model and Relaxation Parameters');

% for i=1:length(wn_mu)
%     Hz_j=0;
%     for j=1:length(wn_z)
%         b_r=quadpoints.z.xp(j,1)*ones(size(zs(:,:,:,i)))+real(zs(:,:,:,i));
%         b_i=quadpoints.z.xp(j,2)*ones(size(zs(:,:,:,i)))+imag(zs(:,:,:,i));
%         b_r=xmz(jr,1)*ones(size(squeeze(KspaceCenter(iii,jjj,:,:,:,:))))+real(squeeze(K
%         Hz_ktmp=0; Hz_k=0;
% %         for k=1:length(wn_mu)
%             Hz_ktmp=a*wn.*exp(-((b_r-real(squeeze(KspaceCenter(iii,jjj,:,:,:,:)))).^2+(b_i-imag(squeeze(KspaceCenter(iii,jjj,:,:,:,:)))).^2)/c);
%             Hz_k=sum(Hz_ktmp(:));
% %             Hz_k=Hz_k+wn_mu(k)*a*exp(-((b_r-real(zs(:,:,:,k))).^2+(b_i-imag(zs(:,:,:,k))).^2)/c);
% %         end
%         Hz_j=Hz_j+wn_z(j)*log(Hz_k);
%     end
%     Hz_i=Hz_i+wn_mu(i)*Hz_j;
% end
% H_z=-Hz_i;
% MI = H_z-0.5*log((2*pi*exp(1))^2*noisestdd^4);


% for iii = 1:length(delaytime)
%     for jjj = 1:length(echotime)
%         Ez1(iii,jjj)=0;
%         for gp1 = 1:NumGP
%             Ez2(iii,jjj)=0;
%             for gp2 = 1:NumGP
%                 Ez3(iii,jjj)=0;
%                 for gp3 = 1:NumGP
%                     Ez4(iii,jjj)=0;
%                     for gp4 = 1:NumGP
%                         Ez4(iii,jjj)=Ez4(iii,jjj)+wn_t2(gp4,2)*KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4);
%                     end
%                     Ez3(iii,jjj)=Ez3(iii,jjj)+wn_t2(gp3,1)*Ez4(iii,jjj);
%                 end
%                 Ez2(iii,jjj)=Ez2(iii,jjj)+wn_t1(gp2,2)*Ez3(iii,jjj);
%             end
%             Ez1(iii,jjj)=Ez1(iii,jjj)+wn_t1(gp1,1)*Ez2(iii,jjj);
%         end
%     end
% end
% 
% for iii = 1:length(delaytime)
%     for jjj = 1:length(echotime)
%         Hz1(iii,jjj)=0;
%         for gp1 = 1:NumGP
%             Hz2(iii,jjj)=0;
%             for gp2 = 1:NumGP
%                 Hz3(iii,jjj)=0;
%                 for gp3 = 1:NumGP
%                     Hz4(iii,jjj)=0;
%                     for gp4 = 1:NumGP
%                         Hz4(iii,jjj)=Hz4(iii,jjj)+wn_t2(gp4,2)*KspaceCenter(iii,jjj,gp1,gp2,gp3,gp4)^2;
%                     end
%                     Hz3(iii,jjj)=Hz3(iii,jjj)+wn_t2(gp3,1)*Hz4(iii,jjj);
%                 end
%                 Hz2(iii,jjj)=Hz2(iii,jjj)+wn_t1(gp2,2)*Hz3(iii,jjj);
%             end
%             Hz1(iii,jjj)=Hz1(iii,jjj)+wn_t1(gp1,1)*Hz2(iii,jjj);
%         end
%     end
% end
% Hz1=0.5*log(2*pi*2.7183*(Hz1-Ez1.^2));

