
function[objfun]=qalasobjfunND(input,Mmeas,TR,flipAngle,nacq,dt)

ndmeas=ndims(Mmeas);
M0=eval(['input(',repmat(':,',[1,ndmeas-1]),'1)']);
T1=eval(['input(',repmat(':,',[1,ndmeas-1]),'2)']);

star=(1-exp(-TR./T1))./(1-cosd(flipAngle)*exp(-TR./T1));

%% T2 sensitization
% M2
M=eval(['Mmeas(',repmat(':,',[1,ndmeas-1]),'1)']);
Mopt=M;
% M3
M=M0.*star-(M0.*star-M).*exp(-dt(3)./(T1.*star));
% M4
M=M0-(M0-M).*exp(-dt(4)./T1);

%% T1 sensitization
% M5
M=-M;             % Assume perfect inversion? 100 ms inversion time
% M6
M=M0-(M0-M).*exp(-dt(6)./T1);
Mopt=cat(ndmeas,Mopt,M);

% Post T1 sens acquisitions
for iii=1:nacq-2
    M=M0.*star-(M0.*star-M).*exp(-dt(5+2*iii)./(T1.*star));
    M=M0-(M0-M).*exp(-dt(6+2*iii)./T1);
    Mopt=cat(ndmeas,Mopt,M);
end

% Mopt=[M(2),M(6:2:6+2*(nacq-2))];

objfunerror=(Mopt-Mmeas).^2;
objfun=sum(objfunerror(~isnan(objfunerror)));

end