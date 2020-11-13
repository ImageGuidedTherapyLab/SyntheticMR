
for iii=1:5
    for jjj=iii:5:640
        magnimg(:,:,ceil(jjj/5),iii)=dicomread(sprintf('DICOMs/image%d.dcm',jjj));
    end
end

kSpacetest(1,:,:,:)=squeeze(kSpacefull(:,:,64,:,5));

figure;imshow3(abs(squeeze(kSpacetest)),[],[2,6])

[calib emaps]=bart('ecalib',kSpacetest);

sens=bart('slice 4 0',calib);
sens_maps=squeeze(sens);

figure;imshow3(abs(sens_maps),[],[2,6])

recon=bart('pics',kSpacetest,sens);
sense_recon=squeeze(recon);
figure;imshow(abs(sense_recon),[])