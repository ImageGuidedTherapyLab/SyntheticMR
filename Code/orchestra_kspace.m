function[kSpace] = orchestra_kspace(pfilePath)
   
    addpath('~/Documents/MATLAB/orchestra-sdk-1.5-2.matlab/');
%     pfilePath='/rsrch1/ip/egates1/QALAS/20170601/P22016.7';
    
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);
    
    acquiredSlices = pfile.slicesPerPass;
    outputSlices = pfile.reconstructedSlicesPerPass;
    scaleFactor = outputSlices / pfile.scaleFactor3d;
   
    kSpace = zeros(pfile.xRes, pfile.yRes, acquiredSlices, pfile.channels, pfile.echoes, pfile.passes);
    
    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
            for slice = 1:acquiredSlices
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;
                for channel = 1:pfile.channels
                    % Load K-Space
                    kSpace(:,:,slice,channel,echo,pass) = GERecon('Pfile.KSpace', sliceInfo, echo, channel);
                end
            end
            
            % Transform Across Slices
            kSpace(:,:,:,:,echo,pass) = ifft(squeeze(kSpace(:,:,:,:,echo,pass)), outputSlices, 3);
            
            % Scale
            kSpace(:,:,:,:,echo,pass) = kSpace(:,:,:,:,echo,pass) * scaleFactor;
        end
    end
end

