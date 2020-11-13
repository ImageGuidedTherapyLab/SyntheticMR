function[M0,T1,T2,objval]=qalas1pfit(Mmeas1p,TR,TE_T2prep,flipAngle,nacq,dt)

% only pass 4 fitting measurements
% normalize mmeas to m5
% initialize mflex to m2
% initialize t1 to log(m5/m2)
% initialize m0 to m5?
% rescale m0 

% mpass=Mmeas1p(2:end)./Mmeas1p(end);
maxmeas=max(Mmeas1p(2:end));
mpass=Mmeas1p(2:end)/maxmeas;
[~,sortind]=sort(mpass',1,'descend');
% if isempty(xinit); xinit=[mpass(1),mpass(end),(dt(6)-dt(end-4))/log((mpass(end)-mpass(end-1))/(mpass(end)-mpass(1)))]; end;
xinit=[mpass(1),mpass(sortind(1)),(dt(6)-dt(end-4))/log((mpass(sortind(1))-mpass(sortind(2)))/(mpass(sortind(1))-mpass(1)))];
% if isempty(xinit); xinit=[0,1,1]; end;
% options=optimoptions('fmincon','MaxIterations',30000,'MaxFunctionEvaluations',30000);
[xm,objval]=fmincon(@(x) qalasobjfunflex(x,mpass,TR,flipAngle,nacq,dt),xinit,[],[],[],[],[-5.,0.,0.],[5.,5.,3.],[]);%,options);
% xm=fminsearch(@(x) qalasobjfun(x,Mmeas,TR,flipAngle,nacq,dt),xinit);
M0=xm(2)*maxmeas;
T1=xm(3);
% Mend=M0-(M0-Mmeas1p(end)./sind(flipAngle)).*exp(-dt(end)./T1);
% T2=-TE_T2prep./log((Mmeas1p(1)./sind(flipAngle))./Mend);
Mend=M0-(M0-Mmeas1p(end)).*exp(-dt(end)./T1);
T2=-TE_T2prep./log(Mmeas1p(1)./Mend);

end