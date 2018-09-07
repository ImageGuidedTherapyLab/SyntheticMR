
TDin=0:.25:5;
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

for iii=4:length(snrin)
    for jjj=1:length(TDin)
        [logZ{iii,jjj},nest_samples{iii,jjj},post_samples{iii,jjj},prior{iii,jjj}]=mn_qalasnp(TDin(jjj),noisein(iii),500);
    end
end