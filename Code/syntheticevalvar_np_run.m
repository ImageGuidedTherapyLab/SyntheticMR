%% implement full parameter space; vectorize up front; run forloop over syntheticevalvar_np; reshape after

function [varmap,meanmap,medianmap] = syntheticevalvar_np_run(tislabel,visdim,parspace,pslabels,vardecay,plotflag)

% TDpT2=0:.1:3;
% TDpT1=0:.1:3;
% tislabel=1;
% vardecay=0;

vdparspace{1}=squeeze(eval(sprintf('parspace%s',visdim{1})));
vdparspace{2}=squeeze(eval(sprintf('parspace%s',visdim{2})));
[vdpar1,vdpar2]=ndgrid(vdparspace{1},vdparspace{2});
vdpar1=vdpar1(:);
vdpar2=vdpar2(:);
pslength=length(vdpar1);
for iii=1:pslength
    fprintf('Parameter space point: %i of %i',iii,pslength)
    [varmap(iii,:),meanmap(iii,:),medianmap(iii,:),~,~,~] = syntheticevalvar_np (tislabel,[vdpar1(iii),vdpar2(iii)],pslabels,vardecay);
end

varmap=reshape(varmap,[length(vdparspace{1}),length(vdparspace{2}),size(varmap,2)]);
meanmap=reshape(meanmap,[length(vdparspace{1}),length(vdparspace{2}),size(meanmap,2)]);
medianmap=reshape(medianmap,[length(vdparspace{1}),length(vdparspace{2}),size(medianmap,2)]);

if plotflag~=0
    figure; pcolor(TDpT2,TDpT1,varmap(:,:,1)); %shading interp;
    colorbar; xlabel('Acq. Spacing (s)'); ylabel('TD (s)'); title('M0 95% Confidence Interval Width');
    figure; pcolor(TDpT2,TDpT1,varmap(:,:,2)); %shading interp;
    colorbar; xlabel('Acq. Spacing (s)'); ylabel('TD (s)'); title('T1 95% Confidence Interval Width');
    figure; pcolor(TDpT2,TDpT1,varmap(:,:,3)); %shading interp;
    colorbar; xlabel('Acq. Spacing (s)'); ylabel('TD (s)'); title('T2 95% Confidence Interval Width');
end

end