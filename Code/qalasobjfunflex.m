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
Mopt=[M(2),M(6:2:end-1)];
objfun=double(norm(Mopt(2:end)-Mmeas(2:end)));

end