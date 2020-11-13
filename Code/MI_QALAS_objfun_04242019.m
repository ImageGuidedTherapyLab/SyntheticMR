
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun,SNR_GM,SNR_WM]=MI_QALAS_objfun_04242019(xopt,tisinput,acqparam,materialID)

%% Assign Acquisition Parameters
% Default Parameters
flipAngle=acqparam(1);
TR=acqparam(2);
TE_T2prep=acqparam(3);
Tacq=acqparam(4);
TDpT2=acqparam(5);
TDinv=acqparam(6);
nacq=acqparam(7);
TD=acqparam(8:6+nacq);
signu=acqparam(end);

TD=xopt;

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(3,1:2),tisinput(5,1:2)],[tisinput(4,1:2),tisinput(6,1:2)]);
lqp=length(xn{1}(:));

    parfor qp=1:lqp
        %     disp(sprintf('Model eval: %d of %d',qp,lqp))
        dt=[0,TE_T2prep,Tacq,TD(1),0,TDinv,Tacq,TD(2),Tacq,TD(3),Tacq,TD(4),Tacq,TD(5)];
%         dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{3}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    end

    % N-D k-space, no averaging
    nd=ndims(materialID);
    evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
    imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
    materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
    for iii=1:length(wn(:))
        imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:,iii)+(materialIDtemp==2).*Mmodel_WM(:,iii)+(materialIDtemp==3).*Mmodel_WM(:,iii);
%         if optcase~=0
            for jjj=1:size(imspace,1)
                eval(evalstr);
                %                 kspace(jjj,:,:)=fftshift(fftn(squeeze(imspace(jjj,:,:))));
            end
            ksr=real(kspace);
            ksi=imag(kspace);
        
        Ezr=Ezr+ksr*wn(iii);
        Ezi=Ezi+ksi*wn(iii);
        Sigrr=Sigrr+ksr.^2*wn(iii);
        Sigii=Sigii+ksi.^2*wn(iii);
        Sigri=Sigri+ksr.*ksi*wn(iii);
    end
% end

N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500; 
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
% signu=3.4762E-4;
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

    MIobjfun=-sum(MI(:));

            dt=[0,TE_T2prep,Tacq,TD(1),0,TDinv,Tacq,TD(2),Tacq,TD(3),Tacq,TD(4),Tacq,TD(5)];
%         dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,S_GM]=qalas1p(tisinput(1,1),tisinput(1,1),tisinput(3,1),tisinput(5,1),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,S_WM]=qalas1p(tisinput(1,2),tisinput(1,2),tisinput(3,2),tisinput(5,2),TR,TE_T2prep,flipAngle,nacq,dt);
        SNR_GM=S_GM./signu;
        SNR_WM=S_WM./signu;
    
end