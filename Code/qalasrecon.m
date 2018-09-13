
% function[M0pred,T1pred,T2pred,objval]=qalasrecon(Mmeas,xinit,imgdim,TR,TE_T2prep,flipAngle,nacq,dt)
function[M0pred,T1pred,T2pred,objval]=qalasrecon(Mmeas,imgdim,TR,TE_T2prep,flipAngle,nacq,dt)

% smeas=size(squeeze(Mmeas));
smeas=size(Mmeas);
Mmeasvec=reshape(double(Mmeas),[prod(smeas(1:imgdim)),smeas((imgdim+1):end)]);
% xinitvec=reshape(xinit,[prod(smeas(1:imgdim)),smeas((imgdim+1):(end-1)),3]);
% mmvsize=size(Mmeasvec,1);
M0predvec=zeros([prod(smeas(1:imgdim)),1]);
T1predvec=zeros([prod(smeas(1:imgdim)),1]);
T2predvec=zeros([prod(smeas(1:imgdim)),1]);
objvalvec=zeros([prod(smeas(1:imgdim)),1]);
evalind=find(~isnan(Mmeasvec(:,1)));
neval=length(evalind);
Mmeasvec=Mmeasvec(evalind,:);
parfor iii=1:neval
%     if sum(isnan(squeeze(Mmeasvec(iii,:))))>0
%         M0predvec(iii)=0;%nan;
%         T1predvec(iii)=0;%nan;
%         T2predvec(iii)=0;%nan;
%         objvalvec(iii)=0;
%     else
%         [M0predvec(iii),T1predvec(iii),T2predvec(iii),objvalvec(iii)]=qalas1pfit(Mmeasvec(iii,:),xinitvec(iii,:),TR,TE_T2prep,flipAngle,nacq,dt);
        [M0eval(iii),T1eval(iii),T2eval(iii),objeval(iii)]=qalas1pfit(Mmeasvec(iii,:),TR,TE_T2prep,flipAngle,nacq,dt);
%     end
end
M0predvec(evalind)=M0eval;
T1predvec(evalind)=T1eval;
T2predvec(evalind)=T2eval;
objvalvec(evalind)=objeval;
if imgdim>1
    M0pred=reshape(M0predvec,smeas(1:imgdim));
    T1pred=reshape(T1predvec,smeas(1:imgdim));
    T2pred=reshape(T2predvec,smeas(1:imgdim));
    objval=reshape(objvalvec,smeas(1:imgdim));
else
    M0pred=M0predvec;
    T1pred=T1predvec;
    T2pred=T2predvec;
    objval=objvalvec;
end

end


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
[xm,objval]=fmincon(@(x) qalasobjfunflex(x,mpass,TR,flipAngle,nacq,dt),xinit,[],[],[],[],[-5.,0.,0.],[5.,5.,4.],[]);%,options);
% xm=fminsearch(@(x) qalasobjfun(x,Mmeas,TR,flipAngle,nacq,dt),xinit);
M0=xm(2)*maxmeas;
T1=xm(3);
% Mend=M0-(M0-Mmeas1p(end)./sind(flipAngle)).*exp(-dt(end)./T1);
% T2=-TE_T2prep./log((Mmeas1p(1)./sind(flipAngle))./Mend);
Mend=M0-(M0-Mmeas1p(end)).*exp(-dt(end)./T1);
T2=-TE_T2prep./log(Mmeas1p(1)./Mend);

end


function[objfun]=qalasobjfunflex(input,Mmeas,TR,flipAngle,nacq,dt)

Mflex=input(1);
M0=input(2);
T1=input(3);

star=(1-exp(-TR./T1))./(1-cosd(flipAngle).*exp(-TR./T1));
%T1 sensitization
M(5)=Mflex;
M(6)=M0-(M0-M(5)).*exp(-dt(6)./T1);
% M(6)=Mmeas(2);

% Post T1 sens acquisitions
for iii=1:nacq-1
    M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

% Mopt=sind(flipAngle).*[M(2),M(6:2:end-1)];
% Mopt=[M(2),M(6:2:end-1)];
% objfun=norm(Mopt(2:end)-Mmeas(2:end));
Mopt=M(6:2:end-1);
objfun=norm(Mopt-Mmeas);

end
