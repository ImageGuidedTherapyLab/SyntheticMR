%% Measurement model
clear all
close all
format shortg

%% Input files
labelfilename = 'ICBM_grey_white_csf.nii.gz'; % population tissue segmentation

%% signal model parameters
T1mean = [100   , 120   , 320   , 250  ]; % ms
T1stdd = [20    ,  20   ,  20   ,  50  ]; % ms
T2mean = [100   , 120   , 320   , 250  ]; % ms
T2stdd = [20    ,  20   ,  20   ,  50  ]; % ms
M0mean = [0.7   , 0.7   , 1.0   , 1.2  ]; % relative intensity
M0stdd = [.2    ,  .2   ,  .3   ,  .4  ]; % ms
flipangle = 120/180.*pi; % radian
excitationangle =  90/180.*pi; % radian
repetitiontime  = 4.0; % seconds
delaytime       = [0:.4:4 ] ; % seconds  
echotime        = [0:14:140 ] ; % ms  

%% Loading tissue types
disp('loading tissue types');
tissuelabel  = load_untouch_nii(labelfilename );
materialID = int32(tissuelabel.img);

%% Generate parameters maps from mean
%% FIXME - @dmitchell412 need quadrature points here
T1map = T1mean(1) * (materialID == 1) + T1mean(2) * (materialID == 2) + T1mean(3) * (materialID == 3) ; 
T2map = T2mean(1) * (materialID == 1) + T2mean(2) * (materialID == 2) + T2mean(3) * (materialID == 3) ; 
M0map = M0mean(1) * (materialID == 1) + M0mean(2) * (materialID == 2) + M0mean(3) * (materialID == 3) ; 


%% loop over delay time and echo time
for iii = 1:length(delaytime )
  for jjj = 1:length(echotime)
    disp( sprintf( '%d of %d, %d of %d',iii ,length(delaytime ) ,jjj,length(echotime        )))
    eTD = exp(-delaytime(iii)/T1map );
    eTR = exp(-repetitiontime/T1map );
    eTE = exp(-echotime(jjj)/T2map );
    % compute magnetization image in image domain
    Mmodel = M0map.*( 1 - (1 - cos(flipangle)) * eTD  - cos(flipangle) * eTD ).*((1- cos(flipangle)*eTR*cos(excitationangle) ).^(-1) ).*eTE ;
    % tranform to fourier space. note this is a function of TD and TE
    Kspace = fft(Mmodel) 
    % TODO - @dmitchell412 use k-space center at every quadrature point to compute MI
  end
end
