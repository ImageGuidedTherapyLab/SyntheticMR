function[bartrecon,kSpace,sensmag,senspha] = OrchBARTQALAS3DRecon(pfilePath)
%% CartesianRecon - Reconstruct 3D Cartesian K-Space
%
% Load Pfile

%     pfilePath='/rsrch1/ip/egates1/QALAS/20170601/P22016.7';
%     pfilePath='/rsrch1/ip/egates1/QALAS/20180427/P11264.7';

if nargin==0
%     pfilePath='/home/dmitchell412/QALASData/P20480.7';
    archive = GERecon('Archive.Load', '/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849.h5');
end
% pfile = GERecon('Pfile.Load', pfilePath);
% header = GERecon('Pfile.Header', pfile);

% acquiredSlices = header.RawHeader.nslices;% pfile.slicesPerPass;
% outputSlices = pfile.reconstructedSlicesPerPass;
% scaleFactor = outputSlices / pfile.scaleFactor3d;

xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;

nEchoes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nechoes;

pass = 1;
zRes = archive.SlicesPerPass(pass);

% kSpace = zeros(acquiredSlices, pfile.xRes, pfile.yRes, pfile.channels, pfile.echoes, pfile.passes);
kspace = complex(zeros(xRes, yRes, nChannels, zRes));

for pass = 1:pfile.passes
    for echo = pfile.echoes:-1:1
        
        
        for slice = 1:acquiredSlices
            
            sliceInfo.pass = pass;
            sliceInfo.sliceInPass = slice;
            
            for channel = 1:pfile.channels
                
                % Load K-Space
                ksptmp(slice,:,:,channel) = GERecon('Pfile.KSpace', sliceInfo, echo, channel);
                
            end
        end
        
        % espirit/pics reconstruction
        %             ksptmp=permute(ksptmp,[3,1,2,4]);
        if echo==pfile.echoes
            [calibtmp,emapstmp]=bart('ecalib -r 20',ksptmp);  % one sensitivity map for all echoes
        end
        sensmagtmp=bart('slice 4 0',calibtmp);
        sensphatmp=bart('slice 4 1',calibtmp);
        recontmp=bart('pics',ksptmp,sensmagtmp);
        
        kSpace(:,:,:,:,echo,pass)=ksptmp;
        sensmag(:,:,:,:,echo,pass)=sensmagtmp;
        senspha(:,:,:,:,echo,pass)=sensphatmp;
        bartrecon(:,:,:,:,echo,pass)=recontmp;
        
    end
    
    % QALAS Fit
    
        flipAngle = 4;           % deg
        TR = 0.005;              % s
        TE_T2prep = 0.100;       % s
        Tacq = 0.500;            % s
        TDpT2 = 0.4;             % s
        TDinv = 0.03;            % s
        nacq = 5;
        TD = [0.5,0.5,0.5,0.5];          % s
        dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];

    smeas=size(squeeze(bartrecon(:,:,:,:,:,pass)));
    Mmeasvec=reshape(double(abs(bartrecon(:,:,:,:,:,pass))),[prod(smeas(1:3)),smeas(4:end)]);
    mmvsize=size(Mmeasvec,1);
    parfor iii=1:size(Mmeasvec,1)
        if sum(isnan(squeeze(Mmeasvec(iii,:))))>0
            M0predvec(iii)=0;%nan;
            T1predvec(iii)=0;%nan;
            T2predvec(iii)=0;%nan;
        else
            [M0predvec(iii),T1predvec(iii),T2predvec(iii)]=qalasrecon(squeeze(Mmeasvec(iii,:)),TR,TE_T2prep,flipAngle,nacq,dt);
        end
    end
    M0pred(:,:,:)=reshape(M0predvec,smeas(1:3));
    T1pred(:,:,:)=reshape(T1predvec,smeas(1:3));
    T2pred(:,:,:)=reshape(T2predvec,smeas(1:3));
%     [M0,T1,T2]=qalasrecon(Mmeas,TR,TE_T2prep,flipAngle,nacq,dt)

for iii=1:20
    figure;
    subplot(2,2,1);
        imagesc(real(squeeze(kSpace(iii,:,:,1,5)))); colormap('gray'); axis off;
    subplot(2,2,2);
        imagesc(real(squeeze(bartrecon(iii,:,:,1,5)))); colormap('gray'); axis off;
    subplot(2,2,3);
        imagesc(imag(squeeze(kSpace(iii,:,:,1,5)))); colormap('gray'); axis off;
    subplot(2,2,4);
        imagesc(abs(squeeze(bartrecon(iii,:,:,1,5)))); colormap('gray'); axis off;
end

figure;
for iii=1:14
    subplot(2,7,iii);
    imagesc(abs(squeeze(sensmag(1,:,:,iii,1)))); axis off;
end
    

figure;
for iii=1:20
    subplot(4,6,iii);
    imagesc(squeeze(M0pred(iii,:,:))); axis off; colorbar;
end

figure;
for iii=1:20
    subplot(4,6,iii);
    imagesc(squeeze(T1pred(iii,:,:))); axis off; colorbar;
end

figure;
for iii=1:20
    subplot(4,6,iii);
    imagesc(squeeze(T2pred(iii,:,:)),[0,1]); axis off; colorbar;
end
