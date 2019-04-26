
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun]=MI_QALAS_objfun_04232019(xopt,tisinput,synthdataT1,synthdataT2,synthdataM0,noise,acqparam,materialID,pdarg,B1inhomflag,filename)

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

TD=xopt;

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(3,1:2),tisinput(5,1:2)],[tisinput(4,1:2),tisinput(6,1:2)]);
lqp=length(xn{1}(:));
    parfor qp=1:lqp
        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
        [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{3}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
        [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    end
    
    for iii=1:5
        for jjj=1:625
    imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(iii,jjj)+(materialIDtemp==2).*Mmodel_WM(iii,jjj)+(materialIDtemp==3).*Mmodel_WM(iii,jjj);
    kspace(:,:,iii,jjj)=fftshift(fftn(imspace));
        end
    end
    
    [~,xn2,~,~,wn2]=GaussHermiteNDGauss(NumQP,[0,0],[signu,signu]);
% S2=repmat(kspace,[1,1,1,1,5,5])+repmat(real(xn2{1}),[151,181,5,625,1,1])+repmat(imag(xn2{2}),[151,181,5,625,1,1]);
zmod=xn2{1}+i*xn2{2};
    


lnterm=log(sum(repmat(wn(:)',[size(xn2{1},1),1]).*S2,2));
pterm=sum(repmat(wn(:)',[size(xn2{1},1),1]).*S2,2);
hz=-sum(wn2.*pterm.*lnterm);

    % N-D k-space, no averaging
    nd=ndims(materialID);
    evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
    imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
    materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
    for iii=1:length(wn(:))
        imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:,iii)+(materialIDtemp==2).*Mmodel_WM(:,iii)+(materialIDtemp==3).*Mmodel_CSF(:,iii);
%         if optcase~=0
            for jjj=1:size(imspace,1)
                eval(evalstr);
            end
            ksr=real(kspace);
            ksi=imag(kspace);

        Ezr=Ezr+ksr*wn(iii);
        Ezi=Ezi+ksi*wn(iii);
        Sigrr=Sigrr+ksr.^2*wn(iii);
        Sigii=Sigii+ksi.^2*wn(iii);
        Sigri=Sigri+ksr.*ksi*wn(iii);
    end

N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500; 
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
signu=3.4762E-4;

detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

    MIobjfun=-sum(MI(:));

end