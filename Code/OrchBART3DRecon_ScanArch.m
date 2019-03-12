function[realimg,bartrecon,kSpace,sensmag,sensphs] = OrchBART3DRecon_ScanArch(scanArchivePath)

% scanArchivePath1='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
% scanArchivePath2='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';

scanArchivePath='/home/khwang/Documents/MATLAB/qalasCS/data/20181004phant/ScanArchive_713792AMR16_20181004_162834615';

if nargin==0
    scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
end
if nargin==1
    if scanArchivePath==1
        scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
    elseif scanArchivePath==2
        scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_144530716';
    elseif scanArchivePath==3
        scanArchivePath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/ScanArchive_SCRB3MR_20190124_163432701';
    elseif scanArchivePath==4
        scanArchivePath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/ScanArchive_SCRB3MR_20190124_163251485';
    elseif scanArchivePath==5
        scanArchivePath='/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/Scan Data/ScanArchive_SCRB3MR_20190124_163430054';
    else
        scanArchivePath='/home/dmitchell412/QALASData/ScanArchive_713792AMR16_20180808_141752849';
    end
end
archive = GERecon('Archive.Load', [scanArchivePath,'.h5']);

% Scan parameters
xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;

% Keep track of the current pass
pass = 1;
zRes = archive.SlicesPerPass(pass);
% nPass = archive.Passes;

% Keep track of the current echo
echo = 1;
nEcho = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nechoes;

% Allocate K-space
kspace = complex(zeros(xRes, yRes, nChannels, zRes, nEcho));

% Loop through each control, sorting frames if applicable
for i = 1:archive.ControlCount
    
    control = GERecon('Archive.Next', archive);
    
    % Sort only programmable packets in range
    if(control.opcode == 1 && ...
            control.viewNum > 0 && ...
            control.viewNum <= yRes && ...
            control.sliceNum < zRes && ...
            control.echoNum < nEcho)
        
        kspace(:,control.viewNum,:,control.sliceNum+1,control.echoNum+1) = control.Frame.Data;
        
    elseif(control.opcode == 0) % end of pass and/or scan
        
    end
end

% espirit/pics reconstruction
kSpace=permute(circshift(kspace,floor(size(kspace,4)/4),4),[4,1,2,3,10,5,6,7,8,9]);
[calib,emaps]=bart('ecalib -r 20',kSpace(:,:,:,:,nEcho)); % one sensitivity map for all echoes
sensmag=bart('slice 4 0',calib);
sensphs=bart('slice 4 1',calib);
bartrecon=bart('pics',kSpace,sensmag);
bartrecon=squeeze(permute(fftshift(bartrecon,1),[2,3,1,6,4,5])); % fftshift correction

% phase correction
magnimg=abs(bartrecon);
phaseimg=angle(bartrecon);
phasecorr=phaseimg-repmat(phaseimg(:,:,:,5),[1,1,1,5]);
realimg=magnimg.*cos(phasecorr);
imagimg=magnimg.*sin(phasecorr);

% Save results
save([scanArchivePath,'.mat'],'realimg','bartrecon','kSpace','sensmag','sensphs','-v7.3');
nii=make_nii(realimg(:,:,:,5));
save_nii(nii,[scanArchivePath,'_realimg.nii']);

% Close the archive to release it from memory
GERecon('Archive.Close', archive);

end