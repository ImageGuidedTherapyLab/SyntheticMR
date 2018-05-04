function[finalImage,kSpace] = OrchestraCartesian3DRecon()
%% CartesianRecon - Reconstruct 3D Cartesian K-Space
%
% Copyright 2016 General Electric Company. All rights reserved.
% GE Proprietary and Confidential Information. Only to be distributed with
% permission from GE. Resulting outputs are not for diagnostic purposes.
%
% Cartesian3DRecon(pfilePath)
% will reconstruct the 3D Cartesian K-Space in the given pfile. This 
% excludes pfiles with ARC enabled.
%
% Limitations: Parallel imaging, intensity correction

    % Load Pfile
    
%     pfilePath='/rsrch1/ip/egates1/QALAS/20170601/P22016.7';
    pfilePath='/rsrch1/ip/egates1/QALAS/20180427/P11264.7';

    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);
    
    acquiredSlices = pfile.slicesPerPass;
    outputSlices = pfile.reconstructedSlicesPerPass;
    scaleFactor = outputSlices / pfile.scaleFactor3d;
   
    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
    
            kSpace = zeros(pfile.xRes, pfile.yRes, acquiredSlices, pfile.channels, pfile.echoes, pfile.passes);

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
            
            % Loop for each slice/channel to create a magnitude image
            for slice = 1:outputSlices
                for channel = 1:pfile.channels

                    % Transform K-Space
                    channelImage = GERecon('Transform', kSpace(:,:,slice,channel,echo,pass));
                    
                    channelImages(:,:,channel) = channelImage;
                end
            
                % Get slice information (corners and orientation) for this slice location
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;

                info = GERecon('Pfile.Info', sliceInfo);

                % Apply Channel Combination
                combinedImage = GERecon('SumOfSquares', channelImages);

                % Create Magnitude Image
                magnitudeImage = abs(combinedImage);

                % Apply Gradwarp
                gradwarpedImage = GERecon('Gradwarp', magnitudeImage, info.Corners, 'XRMW');

                % Orient the image
                finalImage(:,:,slice,echo,pass) = GERecon('Orient', gradwarpedImage, info.Orientation);

                % Display
                imagesc(finalImage(:,:,slice,echo,pass));

                % Display
                title(['Pass: ' num2str(pass) ' Slice: ' num2str(slice) ' Echo: ' num2str(echo)]);

                % Save DICOMs
                imageNumber = (info.Number-1) * pfile.echoes + echo;
                filename = ['DICOMs/image' num2str(imageNumber) '.dcm'];
                GERecon('Dicom.Write', filename, finalImage(:,:,slice,echo,pass), imageNumber, info.Orientation, info.Corners);

                pause(0.05);
            end
        end
    end
end

