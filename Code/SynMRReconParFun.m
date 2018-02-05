function[M0pred,T1pred]=SynMRReconParFun(Mmeas)
% addpath('/home/dmitchell412/Documents/MATLAB/phantom3d');
% slphan=phantom3d('Shepp-Logan',64);
% slphan=slphan(:,:,27:37);
% slphan=slphan(:,:,round(size(slphan,3)/2));
% phanval=unique(slphan);

% Denoise
Mmeas=Mmeas.*(abs(Mmeas)>0.01*max(abs(Mmeas(:))));
Mmeas(Mmeas==0)=nan;

% Phase Correct


% "True" M0, T1, & T2 values
M0val=[nan,1,0.9,0.8,1];                 %
T1val=[nan,4.000,1.400,1.000,0.100];   % s
T2val=[nan,0.600,0.100,0.075,0.100];   % s

% "True" M0, T1, & T2 maps
% M0=zeros(size(slphan));
% T1=zeros(size(slphan));
% T2=zeros(size(slphan));
% for iii=1:length(phanval)
%     M0(slphan==phanval(iii))=M0val(iii);
%     T1(slphan==phanval(iii))=T1val(iii);
%     T2(slphan==phanval(iii))=T2val(iii);
% end

% Acquisition Parameters
nacq=5;
nrepeat=size(Mmeas,5);
TR=0.0026;                      % s
TE_T2prep=0.1;                  % s
flipAngle=5;                    % degrees
dt=0.1*ones([1,6+2*(nacq-1)]);  % s
dt(end)=3;

% Mtrue=qalas(M0,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
% for lll=2:nrepeat
%     Mtrue(:,:,:,:,lll)=qalas(Mtrue(:,:,:,14,lll-1),M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
% end

% Randomize measurement
% Mrand=Mtrue+randn(size(Mtrue))/1000;
% Mmeas=cat(4,Mrand(:,:,:,2,:),Mrand(:,:,:,6:2:6+2*(nacq-2),:));

% Repeated optimizations
for lll=1:nrepeat
    
    % T2 prediction
    % M(2)=M(1).*exp(-TE_T2prep./T2);
    
%     T2pred(:,:,:,lll)=TE_T2prep./log(Mrand(:,:,:,1,lll)./Mrand(:,:,:,2,lll));
    
    % Optimization solution for M0 and T1 prediction
    xinit=[mean(M0val(2:end)),mean(T1val(2:end))];
    smeas=size(Mmeas);
    Mmeasvec=reshape(Mmeas,[prod(smeas(1:3)),smeas(4:end)]);
    parfor iii=1:size(Mmeasvec,1)
        if sum(isnan(squeeze(Mmeasvec(iii,:,lll))))>0
            M0predvec(iii)=nan;
            T1predvec(iii)=nan;
        else
            xm=fminsearch(@(x) qalasobjfun(x,squeeze(Mmeasvec(iii,:,lll)),TR,flipAngle,nacq,dt),xinit);
            M0predvec(iii)=xm(1);
            T1predvec(iii)=xm(2);
        end
    end
    M0pred(:,:,:,lll)=reshape(M0predvec,smeas(1:3));
    T1pred(:,:,:,lll)=reshape(T1predvec,smeas(1:3));
end

% Optimization averages
M0pred=mean(M0pred,4);
T1pred=mean(T1pred,4);
% T2pred=mean(T2pred,4);
