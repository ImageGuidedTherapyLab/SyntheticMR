
% TDin=0:.25:5;
TDin=[0.28,0.48];
% TDin=[.4,.5;.44,.5;.3,.5;.38,.5;.34,.5;.36,.5;.48,.5;.38,.5;.38,.5;.28,.5];
snrin=[1,10,20,100];

load /rsrch1/ip/dmitchell2/github/SyntheticMR/Code/synthphantom_goldstandards.mat;
flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
dtgld=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,.5,Tacq,.5,Tacq,.5,Tacq,.5];
[Mmeasfullgld,Mmeasgld]=qalas(goldstandardM0,goldstandardM0,goldstandardT1,goldstandardT2,TR,TE_T2prep,flipAngle,nacq,dtgld);
noisein=max(Mmeasgld(:))./snrin;

load('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/syntheticPatientPopulation.mat');

for patID=1:10
    patientmatID_hires=syntheticPatientPop{patID,1}(15:165,20:200,92);
    patientmatID_lores=syntheticPatientPop{patID,2}(4:42,5:50,23);
    patient_tisinput=syntheticPatientPop{patID,3};
%     for iii=4:length(snrin)
        for jjj=1:length(TDin)
            [logZ{patID,jjj},nest_samples{patID,jjj},post_samples{patID,jjj},prior{patID,jjj}]=cond_mn_qalasnp(patientmatID_lores,patient_tisinput,TDin(patID,jjj),noisein(4),100,1);
        end
%     end
save conditionalMIresults3.mat logZ nest_samples post_samples prior -v7.3;
end