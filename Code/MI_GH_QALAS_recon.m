%% Measurement model
% clear all
% close all
% format shortg

function [kspace,MI,Hz,varmap,meanmap,medianmap]=MI_GH_QALAS_recon(signal_lib,wn_t1_lib,wn_t2_lib,Msize,visdim,resized,repflag,parspace,pslabels,mkfig)

%% Input files
labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation

%% Loading tissue types
disp('loading tissue types');
switch resized
    case 1
        materialID = int32(1);
    case 1.5
        materialID = int32([1,2,3]);
    case 2
        tissuelabel = load_untouch_nii(labelfilename);
        materialID=int32(tissuelabel.img(15:165,20:200,92));
    case 2.5
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
        tissuelabel = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tissuelabel.img);
    case 0
        system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
        tissuelabel = load_untouch_nii('resampleimg.nii.gz');
        system('rm resampleimg.nii.gz');
        materialID = int32(tissuelabel.img);
    otherwise
        tissuelabel = load_untouch_nii(labelfilename);
        materialID = int32(tissuelabel.img);
end

switch resized
    case 1
        kspace(1,:,:,:,:,:,:,:,:)=reshape(signal_lib(:,1),Msize);
    case 1.5
        kspace(1,:,:,:,:,:,:,:,:,:)=reshape(signal_lib,[Msize,size(signal_lib,2)]);
        kspace=sum(kspace,ndims(kspace))./ndims(kspace);
    case 2
        siglib=reshape(signal_lib(:,1),Msize);
        siglib=squeeze(eval(sprintf('siglib(%s,:)',visdim)));
        ssiglib=size(siglib);
        siglib=siglib(:);
        lsiglib=length(siglib);
        materialID=(materialID~=0);
        parfor iii=1:lsiglib
            disp(sprintf('Fourier transform: %d of %d',iii,lsiglib))
            kspace(:,:,iii) = fftshift(fftn(siglib(iii)*materialID));
        end
    case 2.5
        
    otherwise
        siglib=reshape(signal_lib(:,1),Msize);
        siglib=squeeze(eval(sprintf('siglib(%s,:)',visdim)));
        ssiglib=size(siglib);
        siglib=siglib(:);
        lsiglib=length(siglib);
        materialID=(materialID~=0);
        parfor iii=1:lsiglib
            disp(sprintf('Fourier transform: %d of %d',iii,lsiglib))
            kspace(:,:,:,iii) = fftshift(fftn(siglib(iii)*materialID));
        end
end

disp('Performing quadrature...')
switch resized
    case 1
        ksr=real(kspace);
        ksi=imag(kspace);
        wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
        wn_k=wn_k(:);
        eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
        
        Ezr=sum(ksr.*wnmult,ndims(kspace));
        Ezi=sum(ksi.*wnmult,ndims(kspace));
        Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
        Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
        Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
    case 1.5
        ksr=real(kspace);
        ksi=imag(kspace);
        wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
        wn_k=wn_k(:);
        eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
        
        Ezr=sum(ksr.*wnmult,ndims(kspace));
        Ezi=sum(ksi.*wnmult,ndims(kspace));
        Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
        Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
        Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
    case 2
        kspace=reshape(kspace,[size(kspace,1),size(kspace,2),ssiglib]);
        ksr=real(kspace);
        ksi=imag(kspace);
        wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
        wn_k=wn_k(:);
        eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
        
        Ezr=sum(ksr.*wnmult,ndims(kspace));
        Ezi=sum(ksi.*wnmult,ndims(kspace));
        Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
        Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
        Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
    case 2.5
        wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
        wn_k=wn_k(:);
        
        siglib=reshape(signal_lib,[Msize,size(signal_lib,2)]);
        siglib=squeeze(eval(sprintf('siglib(%s,:,:)',visdim)));
        ssiglib=size(siglib);
        siglib=reshape(siglib,[prod(ssiglib(1:end-2)),ssiglib(end-1),ssiglib(end)]);

        Ezr=zeros([size(materialID),prod(ssiglib(1:end-2))]);
        Ezi=Ezr; Sigrr=Ezr; Sigii=Ezr; Sigri=Ezr;
        sizeijk=ssiglib(end-1);
        sizeh=prod(ssiglib(1:end-2));
        for kkk=1:sizeijk
            for jjj=1:sizeijk
                for iii=1:sizeijk
                    tic;
                    disp(sprintf('FFT/Quadrature: %d of %d, %d of %d, %d of %d',kkk,ssiglib(end-1),jjj,ssiglib(end-1),iii,ssiglib(end-1)))
                    for hhh=1:sizeh
                        disp(sprintf('FFT/Quadrature: %d of %d, %d of %d, %d of %d, %d/%d',kkk,ssiglib(end-1),jjj,ssiglib(end-1),iii,ssiglib(end-1),hhh,prod(ssiglib(1:end-2))))
                        kspace(:,:,:,hhh) = fftshift(fftn(siglib(hhh,iii,1).*(materialID==1) + siglib(hhh,jjj,2).*(materialID==2) + siglib(hhh,kkk,3).*(materialID==3)));
                    end
                    toc;
                    w=wn_k(iii)*wn_k(jjj)*wn_k(kkk);
                    Ezr=Ezr+w.*real(kspace);
                    Ezi=Ezi+w.*imag(kspace);
                    Sigrr=Sigrr+w.*real(kspace).^2;
                    Sigii=Sigii+w.*imag(kspace).^2;
                    Sigri=Sigri+w.*real(kspace).*imag(kspace);
                end
            end
        end
        Ezr=reshape(Ezr,[size(materialID),ssiglib(1:end-2)]);
        Ezi=reshape(Ezi,[size(materialID),ssiglib(1:end-2)]);
        Sigrr=reshape(Sigrr,[size(materialID),ssiglib(1:end-2)]);
        Sigii=reshape(Sigii,[size(materialID),ssiglib(1:end-2)]);
        Sigri=reshape(Sigri,[size(materialID),ssiglib(1:end-2)]);
    otherwise
        kspace=reshape(kspace,[size(kspace,1),size(kspace,2),size(kspace,3),ssiglib]);
        ksr=real(kspace);
        ksi=imag(kspace);
        wn_k=kron(wn_t1_lib(:,1),wn_t2_lib(:,1)');
        wn_k=wn_k(:);
        eval(sprintf('wnmult(%s:)=wn_k;',repmat('1,',[1,ndims(kspace)-1])));
        
        Ezr=sum(ksr.*wnmult,ndims(kspace));
        Ezi=sum(ksi.*wnmult,ndims(kspace));
        Sigrr=sum(ksr.^2.*wnmult,ndims(kspace));
        Sigii=sum(ksi.^2.*wnmult,ndims(kspace));
        Sigri=sum(ksr.*ksi.*wnmult,ndims(kspace));
end

N=2;
signu=1E-4;
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

% save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas5acq_mi_results_nq%d.mat',NumQP),'kspace','MI','Hz','-v7.3');

%% Reproducibility Analysis
if repflag==1
    [varmap{1},meanmap{1},medianmap{1}] = syntheticevalvar_np_run(1,vdrep,parspace,pslabels,vardecay,plotflag);
    [varmap{2},meanmap{2},medianmap{2}] = syntheticevalvar_np_run(2,vdrep,parspace,pslabels,vardecay,plotflag);
    [varmap{3},meanmap{3},medianmap{3}] = syntheticevalvar_np_run(3,vdrep,parspace,pslabels,vardecay,plotflag);
else
    varmap=0;
    meanmap=0;
    medianmap=0;
end

%% Figures
if mkfig==1
    
    if resized~=1
        % Tissue labels
        figure; imagesc(materialID); h=colorbar; ylabel(h,'Tissue Labels');
        xlabel('x'), ylabel('y'); title('Tissue Labels');
        set(gca,'xticklabel',[],'yticklabel',[]);
        saveas(gcf,'MItissuelabels.png');
        
        % MI in k-space
        iii=11; jjj=11;
        mibounds=MI(:,:,1,2,iii,jjj);
        figure; imagesc(mibounds); h=colorbar; ylabel(h,'Mutual Information');
        xlabel('k_x'); ylabel('k_y'); title(sprintf('Mutual Information for TE = %d ms, TD = %d ms',Tacq(iii)*1000,T1TD(jjj)*1000));
        set(gca,'xticklabel',[],'yticklabel',[]);
        saveas(gcf,sprintf('MIkspace_%s_%s.png',iii,jjj));
        
        % MI threshold in k-space
        threshpct=50;
        [misort,sortind]=sort(mibounds(:),'descend');
        minorm=misort-misort(end);
        minorm=minorm/sum(minorm);
        mipct=cumsum(minorm);
        mithresh=min(find(mipct>(threshpct/100)));
        mifig=zeros(size(sortind));
        mifig(sortind(1:mithresh))=1;
        mifig=reshape(mifig,size(mibounds));
        figure; imagesc(mifig); xlabel('k_x'); ylabel('k_y'); title(sprintf('%d%% of Information Content',threshpct));
        set(gca,'xticklabel',[],'yticklabel',[]);
        saveas(gcf,sprintf('MI%sinfo.png',threshpct));
        
        for jjj=1:lt1td
            for iii=1:lt1te
                mibounds=MI(:,:,1,2,iii,jjj);
                figure; imagesc(mibounds); h=colorbar; ylabel(h,'Mutual Information');
                xlabel('k_x'); ylabel('k_y'); title(sprintf('Mutual Information for TE = %d ms, TD = %d ms',Tacq(iii)*1000,T1TD(jjj)*1000));
                set(gca,'xticklabel',[],'yticklabel',[]); caxis([min(MI(:)),max(MI(:))]);
                saveas(gcf,sprintf('MIkspace_%s_%s.png',iii,jjj));
            end
        end
    end
    
    % MI in parameter space
    switch resized
        case 1
        case 2
        otherwise
            figure;
            for iii=1:size(MI,4)
                kx=ceil(size(MI,1)/2); ky=ceil(size(MI,2)/2); kz=ceil(size(MI,2)/2); %76; ky=91;
                miparam=squeeze(MI(kx,ky,1,2,:,:));
                subplot(2,3,iii); imagesc(squeeze(MI(kx,ky,kz,iii,end:-1:1,:))); %caxis([min(MI(:)),max(MI(:))]);
                xlabel('TD (ms)'); ylabel('Flip Angle (deg)'); %title('Mutual Information in Parameter Space at k-Space Center');
%                 xticklabels({'250','500','750','1000'}); xticks([6:5:21]);
%                 yticklabels({'200','180','160','140','120','100','80','60','40','20'}); yticks([1:2:19]);
            end
                h=colorbar; ylabel(h,'Mutual Information');
                saveas(gcf,'MIparam_kcenter.png');
    end
    
end