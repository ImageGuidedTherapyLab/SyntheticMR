
function[M,Mmeas]=qalas(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt)

if ndims(Minit)>1
    imsize=size(Minit);
    Minit=reshape(Minit,[prod(imsize),1]);
    M0=reshape(M0,[prod(imsize),1]);
    T1=reshape(T1,[prod(imsize),1]);
    T2=reshape(T2,[prod(imsize),1]);
end

[M,Mmeas]=qalas1(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
for iii=1:100
    [M,Mmeas]=qalas1(M(:,end),M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
    merr=abs(M(:,1)-M(:,end));
    if max(merr(:))<=eps; break; end;
end

if size(Minit,1)>1
    M=reshape(M,[imsize,size(M,2)]);
    Mmeas=reshape(Mmeas,[imsize,size(Mmeas,2)]);
end

end


function[M,Mmeas]=qalas1(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt)

star=(1-exp(-TR./T1))./(1-cosd(flipAngle).*exp(-TR./T1));

% T2 sensitization
M(:,1)=Minit;        % assume initialization?
M(:,2)=M(:,1).*exp(-TE_T2prep./T2);
M(:,3)=M0.*star-(M0.*star-M(:,2)).*exp(-dt(3)./(T1.*star));
M(:,4)=M0-(M0-M(:,3)).*exp(-dt(4)./T1);

% T1 sensitization
M(:,5)=-M(:,4);             % Assume perfect inversion? 100 ms inversion time
M(:,6)=M0-(M0-M(:,5)).*exp(-dt(6)./T1);

% Post T1 sens acquisitions
for iii=1:nacq-1
    M(:,5+2*iii)=M0.*star-(M0.*star-M(:,4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(:,6+2*iii)=M0-(M0-M(:,5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

% M=sind(flipAngle).*M;
% Mmeas=sind(flipAngle).*cat(4,M(:,2),M(:,6:2:end-1));
Mmeas=M(:,[2,6:2:end-1]);
% cat(4,M(:,2),M(:,6:2:end-1));

end