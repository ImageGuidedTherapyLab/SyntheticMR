% linear ODE PBPK model to fit with the experimental data

function Mopt = mnqalasmodel(data, parnames, parvals)

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length');
end
nparams = lpn;

% extract the parameters
for ii=1:nparams
  switch parnames{ii}
    case 'M0'
      M0 = parvals{ii};
    case 'T1'
      T1 = parvals{ii};
	case 'T2'
      T2 = parvals{ii};
    case 'TR'
      TR = parvals{ii};
    case 'flipAngle'
      flipAngle = parvals{ii};
  end
end

%M=Q(necho,M0,T1,T2)

%function[M]=qalas(Minit,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt)

% nrepeat=size(Mmeas,5);
nacq=5;
% TR=0.0026;                      % s
TE_T2prep=0.1;                  % s
% flipAngle=5;                    % degrees
dt=0.1*ones([1,6+2*(nacq-1)]);  % s
dt(end)=3;

star=(1-exp(-TR./T1))./(1-cosd(flipAngle)*exp(-TR./T1));

% T2 sensitization
M(1)=0;
M(2)=889;      % Fix first point?
M(3)=M0.*star-(M0.*star-M(2)).*exp(-dt(3)./(T1.*star));
M(4)=M0-(M0-M(3)).*exp(-dt(4)./T1);

% T1 sensitization
M(5)=-M(4);             % Assume perfect inversion? 100 ms inversion time
M(6)=M0-(M0-M(5))*exp(-dt(6)/T1);

% Post T1 sens acquisitions
for iii=1:nacq-2
    M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
    M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
end

Mopt=[M(2);M(6:2:6+2*(nacq-2))'];

end