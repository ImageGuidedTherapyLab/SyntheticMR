
for iii=1:70
    t1img(:,:,iii)=dicomread(sprintf('t1mappingscan/Image#0%02d.dcm',iii));
end

t1img=4.3217.*double(t1img)-80.504;

nii=make_nii(t1img);
save_nii(nii,'t1mappingscan/t1img.nii');



t1seg=load_untouch_nii('t1mappingscan/t1imgseg.nii.gz');
t1seg=t1seg.img;

for iii=1:42
    t1mean(iii)=nanmean(t1img(t1seg==iii));
    t1std(iii)=nanstd(t1img(t1seg==iii));
    t1n(iii)=sum(t1seg(:)==iii);
end

figure; errorbar(t1mean(1:14),t1std(1:14),'o');
xlabel('PD Element'); ylabel('T1 (ms)');
figure; errorbar(t1mean(15:28),t1std(15:28),'o');
xlabel('T1 Element'); ylabel('T1 (ms)');
figure; errorbar(t1mean(29:42),t1std(29:42),'o');
xlabel('T2 Element'); ylabel('T1 (ms)');
