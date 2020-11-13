
% kspsize=[size(kSpacetest1,1),size(kSpacetest1,2)];
kspsize=[size(kSpacetest1,2),size(kSpacetest1,3)];
% spacing=[10,5,4,3:-.25:2.25,2:-.1:0.5];
spacing=[10,2.25,1.6,1.4,1.2:-.1:.5];
fullsamplesize=24;
fss1=ceil((fullsamplesize-1)/2);
fss2=floor((fullsamplesize-1)/2);

nrepeat=1;

samplepct=zeros([length(spacing),nrepeat]);
bartacc=zeros([length(spacing),nrepeat]);
M0reconacc=zeros([length(spacing),nrepeat]);
T1reconacc=zeros([length(spacing),nrepeat]);

% Fully sampled recon

kSpace4d=bart('join 10',kSpacetest1,kSpacetest2,kSpacetest3,kSpacetest4,kSpacetest5);
sens4d=bart('join 10',sens1,sens2,sens3,sens4,sens5);
discsample=ones(size(kSpace4d(:,:,:,1,:,:,:,:,:,:,:)));
writecfl('ret4d/kSpace4d',kSpace4d);
writecfl('ret4d/sens4d',sens4d);
writecfl('ret4d/discsample',discsample);
system('/opt/apps/bart-0.4.00/bart pics -d5 -i 100 -p ret4d/discsample ret4d/kSpace4d ret4d/sens4d ret4d/recontest4d');
recon4d=readcfl('ret4d/recontest4d');
% recon4d=bart('pics -i 100 -R L:7:7:0.05',kSpace4d,sens4d);

Mmeas=squeeze(recon4d);%cat(4,recon1,recon2,recon3,recon4,recon5);
Mmeasstore{1,1}=Mmeas;
for iii=1:size(Mmeas,4)
    Mfull(:,:,:,iii)=medfilt3(abs(Mmeas(:,:,:,iii)));
end
[M0pred,T1pred]=SynMRReconParFun(abs(Mmeas));

M0store{1,1}=M0pred;
T1store{1,1}=T1pred;

M0pred=medfilt3(M0pred);
T1pred=medfilt3(T1pred);

% Undersampled recons
for kkk=1:nrepeat
    
    for jjj=1:length(spacing)
        jjj
        for pdind=1:size(kSpace4d,11)
        pts=poissonDisc(kspsize,spacing(jjj));
        pts=round(pts);
        discsample=zeros(kspsize);
        for iii=1:size(pts,1)
            discsample(pts(iii,1),pts(iii,2))=1;
        end
        discsample(ceil(kspsize(1)/2)-fss1:ceil(kspsize(1)/2)+fss2,ceil(kspsize(2)/2)-fss1:ceil(kspsize(2)/2)+fss2)=1;
        % discsample=repmat(discsample,[1,1,size(kSpacetest,3),size(kSpacetest,4)]);
        discsample2(1,:,:,1,pdind)=discsample;
        end
        discsample=repmat(discsample2,[size(kSpace4d,1),1,1,size(kSpace4d,4),1]);
        samplepct(jjj,kkk)=sum(discsample(:)==1)/numel(discsample);
        
        for iii=1:size(kSpace4d,11); discsamplebart(:,:,:,1,1,1,1,1,1,1,iii)=discsample(:,:,:,1,iii); end;
        writecfl('ret4d/discsample',discsamplebart);
        system('/opt/apps/bart-0.4.00/bart pics -d5 -i 100 -p ret4d/discsample ret4d/kSpace4d ret4d/sens4d ret4d/recontest4d');
        recon4dsub=readcfl('ret4d/recontest4d');
        
%         kSpace4dsub=bart('join 10',kSpacetest1.*discsample(:,:,:,:,1),kSpacetest2.*discsample(:,:,:,:,2),kSpacetest3.*discsample(:,:,:,:,3),kSpacetest4.*discsample(:,:,:,:,4),kSpacetest5.*discsample(:,:,:,:,5));
%         recon4dsub=bart('pics -i 100 -R L:7:7:0.05',kSpace4dsub,sens4d);
%         recon1sub=bart('pics',kSpacetest1.*discsample,sens1);
%         recon2sub=bart('pics',kSpacetest2.*discsample,sens2);
%         recon3sub=bart('pics',kSpacetest3.*discsample,sens3);
%         recon4sub=bart('pics',kSpacetest4.*discsample,sens4);
%         recon5sub=bart('pics',kSpacetest5.*discsample,sens5);
        
        Mmeas=squeeze(recon4dsub);%cat(4,recon1sub,recon2sub,recon3sub,recon4sub,recon5sub);
        Mmeasstore{jjj+1,kkk}=Mmeas;
                
        [M0sub,T1sub]=SynMRReconParFun(abs(Mmeas));
        M0store{jjj+1,kkk}=M0sub;
        T1store{jjj+1,kkk}=T1sub;
        
        for iii=1:size(Mmeas,4)
            Mmeas(:,:,:,iii)=medfilt3(abs(Mmeas(:,:,:,iii)));
        end
        M0sub=medfilt3(M0sub);
        T1sub=medfilt3(T1sub);
        barterr=Mfull(:)-Mmeas(:);
        barterrmean=mean(Mfull(~isnan(Mfull)))-Mmeas(:);
        bartacc(jjj,kkk)=norm(barterr)/norm(barterrmean);
        M0reconerr=M0pred(:)-M0sub(:);
        T1reconerr=T1pred(:)-T1sub(:);
        M0reconerrmean=mean(M0pred(~isnan(M0pred)))-M0sub(:);
        T1reconerrmean=mean(T1pred(~isnan(T1pred)))-T1sub(:);
        M0reconacc(jjj,kkk)=norm(M0reconerr(~isnan(M0reconerr)))/norm(M0reconerrmean(~isnan(M0reconerrmean)));
        T1reconacc(jjj,kkk)=norm(T1reconerr(~isnan(T1reconerr)))/norm(T1reconerrmean(~isnan(T1reconerrmean)));
    end
    
    save ReconErrorTest4DResults.mat samplepct bartacc M0reconacc T1reconacc Mmeasstore M0store T1store -v7.3;
    
end