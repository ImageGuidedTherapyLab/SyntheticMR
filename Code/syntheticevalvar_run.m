%% implement full parameter space; vectorize up front; run forloop over syntheticevalvar_np; reshape after

TDpT2=0:.1:3;
TDpT1=0:.1:3;
tislabel=1;
vardecay=0;

for iii=1:length(TDpT2)
    for jjj=1:length(TDpT1)
        iii
        jjj
        [varmap(iii,jjj,:),~,~,~] = syntheticevalvar (tislabel,TDpT2(iii),TDpT1(jjj),vardecay);
    end
end

figure; pcolor(TDpT2,TDpT1,varmap(:,:,1)); %shading interp; 
colorbar; xlabel('Acq. Spacing (s)'); ylabel('TD (s)'); title('M0 95% Confidence Interval Width');
figure; pcolor(TDpT2,TDpT1,varmap(:,:,2)); %shading interp; 
colorbar; xlabel('Acq. Spacing (s)'); ylabel('TD (s)'); title('T1 95% Confidence Interval Width');
figure; pcolor(TDpT2,TDpT1,varmap(:,:,3)); %shading interp; 
colorbar; xlabel('Acq. Spacing (s)'); ylabel('TD (s)'); title('T2 95% Confidence Interval Width');
