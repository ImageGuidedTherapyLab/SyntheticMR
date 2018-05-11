
function[M0,T1,T2]=qalasrecon(Mmeas,TR,TE_T2prep,flipAngle,nacq,dt)

xinit=[0,1,1];
xm=fminsearch(@(x) qalasobjfun(x,Mmeas,TR,flipAngle,nacq,dt),xinit);
M0=xm(2);
T1=xm(3);
Mend=M0-(M0-Mmeas(end)./sind(flipAngle)).*exp(-dt(end)./T1);
T2=-TE_T2prep./log((Mmeas(1)./sind(flipAngle))./Mend);

end


function[objfun]=qalasobjfun(input,Mmeas,TR,flipAngle,nacq,dt)

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

Mopt=sind(flipAngle).*[M(2),M(6:2:end-1)];
objfun=norm(Mopt(2:end)-Mmeas(2:end));

end
