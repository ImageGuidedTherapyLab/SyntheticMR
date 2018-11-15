
% function[M0pred,T1pred,T2pred,objval]=qalasrecon(Mmeas,xinit,imgdim,TR,TE_T2prep,flipAngle,nacq,dt)
function[M0pred,T1pred,T2pred,objval,Mpred]=qalasreconmap(Mmeas,imgdim,TR,TE_T2prep,flipAngle,nacq,dt)

% smeas=size(squeeze(Mmeas));
smeas=size(Mmeas);
Mmeasvec=reshape(double(Mmeas),[prod(smeas(1:imgdim)),smeas((imgdim+1):end)]);
% xinitvec=reshape(xinit,[prod(smeas(1:imgdim)),smeas((imgdim+1):(end-1)),3]);
% mmvsize=size(Mmeasvec,1);

evalind=find(~isnan(Mmeasvec(:,1)));
neval=length(evalind);
Mmeasvec=Mmeasvec(evalind,:);
parfor iii=1:neval
    try
        [M0eval(iii),T1eval(iii),T2eval(iii),objeval(iii),~]=qalas1pfit(Mmeasvec(iii,:),TR,TE_T2prep,flipAngle,nacq,dt);
    catch
        M0eval(iii)=nan;
        T1eval(iii)=nan;
        T2eval(iii)=nan;
        objeval(iii)=nan;
    end
end
M0predvec=zeros([prod(smeas(1:imgdim)),1]);
T1predvec=zeros([prod(smeas(1:imgdim)),1]);
T2predvec=zeros([prod(smeas(1:imgdim)),1]);
objvalvec=zeros([prod(smeas(1:imgdim)),1]);
% Mvec=zeros([prod(smeas(1:imgdim)),size(M,2)]);
M0predvec(evalind)=M0eval;
T1predvec(evalind)=T1eval;
T2predvec(evalind)=T2eval;
objvalvec(evalind)=objeval;
% Mvec(evalind,:)=M;
if imgdim>1
    M0pred=reshape(M0predvec,smeas(1:imgdim));
    T1pred=reshape(T1predvec,smeas(1:imgdim));
    T2pred=reshape(T2predvec,smeas(1:imgdim));
    objval=reshape(objvalvec,smeas(1:imgdim));
%     Mpred=reshape(Mvec,[smeas(1:imgdim),size(M,2)]);
else
    M0pred=M0predvec;
    T1pred=T1predvec;
    T2pred=T2predvec;
    objval=objvalvec;
%     Mpred=Mvec;
end

end


function[M0,T1,T2,objval,M]=qalas1pfit(Mmeas1p,TR,TE_T2prep,flipAngle,nacq,dt)

maxmeas=max(abs(Mmeas1p(2:end)));
mpass=Mmeas1p(2:end)/maxmeas;
% [~,sortind]=sort(mpass',1,'descend');
xinit=[mpass(1),max(abs(mpass)),.5]; %(dt(6)-dt(end-4))/log((mpass(sortind(1))-mpass(sortind(2)))/(mpass(sortind(1))-mpass(1)))];
% if isempty(xinit); xinit=[0,1,1]; end;
% options=optimoptions('fmincon','MaxIterations',30000,'MaxFunctionEvaluations',30000);
[xm,objval]=fmincon(@(x) qalasobjfunflex(x,mpass',TR,flipAngle,nacq,dt),xinit,[],[],[],[],[-5.,0.,0.],[0.,5.,4.],[]);%,options);
M0=xm(2)*maxmeas;
T1=xm(3);
% Mend=M0-(M0-Mmeas1p(end)).*exp(-dt(end)./T1);
% T2=-TE_T2prep./log(Mmeas1p(1)./Mend);

M=qalasflex(xm(1)*maxmeas,M0,T1,TR,flipAngle,nacq,dt);
T2=-TE_T2prep./log(Mmeas1p(1)./M(end));

end


function[objfun]=qalasobjfunflex(input,Mmeas,TR,flipAngle,nacq,dt)

Mflex=input(1);
M0=input(2);
T1=input(3);

star=(1-exp(-TR./T1))./(1-cosd(flipAngle).*exp(-TR./T1));
%T1 sensitization
M(5)=Mflex;
M(6)=M0-(M0-M(5)).*exp(-dt(6)./T1);

% Post T1 sens acquisitions
for iii=1:nacq-1
    M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

Mopt=M(6:2:end-1);
if size(Mopt)~=size(Mmeas); Mopt=Mopt'; end;
objfun=norm(Mopt-Mmeas);

end


function[M]=qalasflex(M5,M0,T1,TR,flipAngle,nacq,dt)

star=(1-exp(-TR./T1))./(1-cosd(flipAngle).*exp(-TR./T1));

% T1 sensitization
M(6)=M0-(M0-M5).*exp(-dt(6)./T1);

% Post T1 sens acquisitions
for iii=1:nacq-1
    M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

% T2 sensitization
M(5)=M5;
M(4)=-M5;
M(3)=M0-(M0-M(4))./exp(-dt(4)./T1);
M(2)=M0.*star-(M0.*star-M(3))./exp(-dt(3)./(T1.*star));
M(1)=M(6+2*iii);

end
