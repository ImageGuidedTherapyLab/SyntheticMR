function[objfun,M]=qalasobjfunss(input,Mmeas,TR,TE_T2prep,flipAngle,nacq,dt)

M0=input(1);
T1=input(2);
T2=input(3);

[M]=qalas1(M0,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
for iii=1:100
    [M]=qalas1(M(end),M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
    merr=abs(M(1)-M(end));
    if max(merr(:))<=eps; break; end;
end

Mopt=[M(2),M(6:2:end-1)];
objfun=norm(Mopt-Mmeas);

end


function[M]=qalas1(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt)

star=(1-exp(-TR./T1))./(1-cosd(flipAngle).*exp(-TR./T1));

% T2 sensitization
M(1)=Minit;
M(2)=M(1).*exp(-TE_T2prep./T2);
M(3)=M0.*star-(M0.*star-M(2)).*exp(-dt(3)./(T1.*star));
M(4)=M0-(M0-M(3)).*exp(-dt(4)./T1);

% T1 sensitization
M(5)=-M(4);
M(6)=M0-(M0-M(5)).*exp(-dt(6)./T1);

% Post T1 sens acquisitions
for iii=1:nacq-1
    M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

end