%% Measurement model
% clear all
% close all
% format shortg

function [kspace,MI,Hz]=MI_GH_QALAS(resized,mkfig)

%% Input files
labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation

%% signal model parameters
%%          GM       WM        CSF     Tumor
T1mean = [1400   , 1000    ,  4000    , 250  ]./1000; % s
T1stdd = [100    ,  100    ,   200    ,  50  ]./1000; % s
T2mean = [100    ,   75    ,   600    , 250  ]./1000; % s
T2stdd = [  5    ,    5    ,    30    ,  50  ]./1000; % s
M0mean = [  0.9  ,    0.8  ,     1.0  ,   1.2];       % relative intensity
M0stdd = [   .05 ,     .05 ,      .05 ,    .4];       % relative intensity

TR = 0.0026;              % s
TE_T2prep = 0.100;        % s
Tacq = 0.100;             % s
nacq=2;
flipAngle = 1:.5:10;      % degrees
TD1 = [0:10:200]./1000;   % s
% TD = repmat(10:50:200,[nacq,1])./1000; % s

%% Quadrature Points
NumQP=5;
[x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,[mean(T1mean(1:2))],[mean(T1stdd(1:2))]);
[x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,[mean(T2mean(1:2))],[mean(T2stdd(1:2))]);
% [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,[mean(T1mean(1:2)),T1mean(3)],[mean(T1stdd(1:2)),T1stdd(3)]);
% [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,[mean(T2mean(1:2)),T2mean(3)],[mean(T2stdd(1:2)),T2stdd(3)]);

%% Loading tissue types
disp('loading tissue types');

if resized==1
    materialID = int32(1);
elseif resized==2
    tissuelabel = load_untouch_nii(labelfilename);
    materialID=int32(tissuelabel.img(15:165,20:200,92));
else
    system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
    tissuelabel = load_untouch_nii('resampleimg.nii.gz');
    system('rm resampleimg.nii.gz');
    materialID = int32(tissuelabel.img);
end
% tissuelabel  = load_untouch_nii(labelfilename);

%% TODO - downsample materialID by factor of 2
% materialID = imresize3(int32(tissuelabel.img),0.25);
% materialID = int32(tissuelabel.img);
% if resized==2
%     materialID=materialID(15:165,20:200,92);
% end

%% loop over delay time and echo time
% lt1te=length(Tacq);
% lt1td=length(T1TD);
[x1,x2]=ndgrid(xm_t1,xm_t2);
x1=x1(:);
x2=x2(:);

[p1,p2]=ndgrid(flipAngle,TD1);
p1=p1(:);
p2=p2(:);
lpp=length(p1);

% for qp = 1:NumQP^2
%     for pp = 1:lpp
%         Mqp(pp,qp)=qalas
%     end
% end

M0map = M0mean(1) * (materialID == 1) + M0mean(2) * (materialID == 2) + M0mean(3) * (materialID == 3) ;
parfor qp = 1:NumQP^2
    % Generate parameters maps from mean
    T1map = eps * (materialID==0) + x1(qp) * or(materialID == 1 , materialID == 2) + x1(qp) * (materialID == 3) ;
    T2map = eps * (materialID==0) + x2(qp) * or(materialID == 1 , materialID == 2) + x2(qp) * (materialID == 3) ;
    
    for pp = 1:lpp
        dt = Tacq*[0,0,1,0,0,0,1,0]+p2(pp)*[0,0,0,0,0,1,0,0];
        disp(sprintf('Model eval: %d of %d, %d of %d',qp,NumQP^2,pp,lpp))
        [~,Mmodel(:,:,:,:,pp,qp)]=qalas(M0map,M0map,T1map,T2map,TR,TE_T2prep,p1(pp),nacq,dt);
    end
end

Msize=size(Mmodel);
if resized==1
    kspace=Mmodel;
else
    Mcat=prod(Msize(4:end));
    Mmodel=reshape(Mmodel,[Msize(1:3),Mcat]);
    parfor iii=1:size(Mmodel,4)
        disp(sprintf('Fourier transform: %d of %d',iii,Mcat))
        kspace(:,:,:,iii) = fftshift(fftn(Mmodel(:,:,:,iii)));
    end
end
kspace=reshape(kspace,[Msize(1:4),length(flipAngle),length(TD1),Msize(end)]);

% save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas_mi_results_nq%d.mat',NumQP),'kspace','-v7.3');

ksize=size(kspace);
ksize=ksize(1:6);
ksr=real(kspace);
ksi=imag(kspace);

wn_k=kron(wn_t1,wn_t2');
wn_k=wn_k(:);

Ezr=zeros(ksize);
Ezi=zeros(ksize);
Sigrr=zeros(ksize);
Sigii=zeros(ksize);
Sigri=zeros(ksize);
parfor iii=1:NumQP^2
    disp(sprintf('Quadrature: %d of %d',iii,NumQP^2))
    Ezr=Ezr+wn_k(iii).*ksr(:,:,:,:,:,:,iii);
    Ezi=Ezi+wn_k(iii).*ksi(:,:,:,:,:,:,iii);
    Sigrr=Sigrr+wn_k(iii).*ksr(:,:,:,:,:,:,iii).^2;
    Sigii=Sigii+wn_k(iii).*ksi(:,:,:,:,:,:,iii).^2;
    Sigri=Sigri+wn_k(iii).*ksr(:,:,:,:,:,:,iii).*ksi(:,:,:,:,:,:,iii);
end

N=2;
signu=0.1;
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MI=Hz-Hzmu;

% save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas_mi_results_nq%d.mat',NumQP),'kspace','MI','Hz','-v7.3');

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
    kx=ceil(size(MI,1)/2); ky=ceil(size(MI,2)/2); %76; ky=91;
    miparam=squeeze(MI(kx,ky,1,2,:,:));
    figure; imagesc(miparam(end:-1:1,:)); %caxis([min(MI(:)),max(MI(:))]);
    xlabel('TE (ms)'); ylabel('TD (ms)'); title('Mutual Information in Parameter Space at k-Space Center');
    xticklabels({'250','500','750','1000'}); xticks([6:5:21]);
    yticklabels({'200','180','160','140','120','100','80','60','40','20'}); yticks([1:2:19]);
    h=colorbar; ylabel(h,'Mutual Information');
    saveas(gcf,'MIparam_kcenter.png');
    
end