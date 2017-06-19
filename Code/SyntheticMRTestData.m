
%% Data Processing

seriesdir='/mnt/FUS4/IPVL_research/ZZZQALASTEST/20170601/1.2.840.113619.6.374.135981113113158878631862266495964625240/1.2.840.113619.2.374.5268979.3611796.17401.1496316543.636/';
[~,dicomstr]=system(sprintf('ls %s',seriesdir));
dicomls=strsplit(dicomstr);
dicomls(end)=[];
dicompath=cellfun(@(x) strcat(seriesdir,x),dicomls,'UniformOutput',false);
dicomimg=cellfun(@dicomread,dicompath,'UniformOutput',false);
dicomhdr=cellfun(@dicominfo,dicompath,'UniformOutput',false);

% Slice Location Indexing
slicelocation=arrayfun(@(x) dicomhdr{x}.SliceLocation,1:numel(dicomhdr));
[~,~,sliceindex]=unique(slicelocation);

% Echo Number Indexing
echonumber=arrayfun(@(x) dicomhdr{x}.EchoNumber,1:numel(dicomhdr));

% Image Type Indexing (0=magnitude, 1=phase, 2=real, 3=imaginary)
imtype=arrayfun(@(x) dicomhdr{x}.Private_0043_102f,1:numel(dicomhdr));

% Index Key
indexkey=[imtype',echonumber',sliceindex];
[~,sortind]=sortrows(indexkey);
imgindexsize=max(echonumber)*max(sliceindex);
sliceindexsize=max(sliceindex);

% Magnitude Images
for iii=1:max(echonumber)
    magnimg(:,:,:,iii)=double(cat(3,dicomimg{sortind(sliceindexsize*(iii-1)+1:sliceindexsize*iii)}));
end

% Real Images
for iii=1:max(echonumber)
    realimg(:,:,:,iii)=double(cat(3,dicomimg{sortind(sliceindexsize*(iii-1)+imgindexsize+1:sliceindexsize*iii+imgindexsize)}));
end

% Imaginary Images
for iii=1:max(echonumber)
    imagimg(:,:,:,iii)=double(cat(3,dicomimg{sortind(sliceindexsize*(iii-1)+2*imgindexsize+1:sliceindexsize*iii+2*imgindexsize)}));
end

% Phase Images
phaseimg=atan(imagimg./realimg);
for iii=1:max(echonumber)
    phaseimgcorrected(:,:,:,iii)=phaseimg(:,:,:,iii)-phaseimg(:,:,:,max(echonumber));
end
imagimgcorrected=realimg.*tan(phaseimgcorrected);
magnimgcorrected=sqrt(realimg.^2+imagimgcorrected.^2);

% Noise Mask
noisemask=magnimg>0.05*max(magnimg(:));
magnimg=magnimg.*noisemask;
realimg=realimg.*noisemask;
imagimg=imagimg.*noisemask;

% imagimg(:,:,1:2:end,:)=-imagimg(:,:,1:2:end,:);

% Phase Correction
for iii=1:max(echonumber)
    imagimgcorrected(:,:,:,iii)=imagimg(:,:,:,iii)-imagimg(:,:,:,max(echonumber));
    magnimgcorrected(:,:,:,iii)=sqrt(realimg(:,:,:,iii).^2+imagimgcorrected(:,:,:,iii).^2);
%     magnimgcorrected(magnimgcorrected<0.05*max(magnimgcorrected(:)))=0;
%     magnimg(:,:,:,iii)=magnimgcorrected;
end

%% Reconstruction

% Acquisition Parameters
nacq=max(echonumber);
nrepeat=1;
TR=dicomhdr{1}.RepetitionTime/1000;                 % s
TE_T2prep=0.048;                                    % s
tinv=dicomhdr{1}.InversionTime/1000;                % s
flipAngle=dicomhdr{1}.FlipAngle;                    % degrees
tacq=32*TR;                                         % s
tdelay=0;                                           % s
dt=[0,TE_T2prep,tacq,tdelay,tinv,0,repmat([tacq,tdelay],[1,nacq-1])];   % s
% dt=0.1*ones([1,6+2*(nacq-1)]);                      % s
% dt(end)=3;

Mmeas=double(magnimg(:,:,round(size(magnimg,3)/2),:));

% Repeated optimizations
for lll=1:nrepeat
    
    % T2 prediction
    % M(2)=M(1).*exp(-TE_T2prep./T2);
    
%     T2pred(:,:,:,lll)=TE_T2prep./log(Mrand(:,:,:,1,lll)./Mrand(:,:,:,2,lll));
T2pred=zeros(size(Mmeas));
    
    % Optimization solution for M0 and T1 prediction
    for iii=1:size(Mmeas,1)
        for jjj=1:size(Mmeas,2)
            for kkk=1:size(Mmeas,3)
                if sum(squeeze(Mmeas(iii,jjj,kkk,:,lll)))==0
                    M0pred(iii,jjj,kkk,lll)=nan;
                    T1pred(iii,jjj,kkk,lll)=nan;
                else
                    %                 [xm,fval,exitflag,output]=fmincon(@(x) qalasobjfun(x,squeeze(Mmeas(iii,jjj,kkk,:))',TR,TE_T2prep,flipAngle,nacq,dt),[mean(M0val(2:end)),mean(T1val(2:end))],[],[],[],[],[0,0],[1.5,5],[],optimset('TolX',1.e-10,'TolFun',1.e-10));%,'Display','iter-detailed'));
                    xm=fminsearch(@(x) qalasobjfun(x,squeeze(Mmeas(iii,jjj,kkk,:,lll))',TR,TE_T2prep,flipAngle,nacq,dt),[1,1]);
                    M0pred(iii,jjj,kkk,lll)=xm(1);
                    T1pred(iii,jjj,kkk,lll)=xm(2);
                end
            end
        end
    end
    
end

% Optimization averages
M0pred=mean(M0pred,4);
T1pred=mean(T1pred,4);
T2pred=mean(T2pred,4);

%% Output Stats
for iii=2:length(phanval)
    M0mean(iii-1)=mean(M0pred(slphan==phanval(iii)));
    M0std(iii-1)=std(M0pred(slphan==phanval(iii)));
    T1mean(iii-1)=mean(T1pred(slphan==phanval(iii)));
    T1std(iii-1)=std(T1pred(slphan==phanval(iii)));
    T2mean(iii-1)=mean(T2pred(slphan==phanval(iii)));
    T2std(iii-1)=std(T2pred(slphan==phanval(iii)));
end
M0mean
M0std
T1mean
T1std
T2mean
T2std

%% Output Images
figure
imagesc(M0(:,:,round(size(M0,3)/2)))
caxis([.7,1.1])
colormap('gray')
title('True M0')
colorbar
figure
imagesc(M0pred(:,:,round(size(M0pred,3)/2)))
caxis([.7,1.1])
colormap('gray')
title('Predicted M0')
colorbar

figure
imagesc(T1(:,:,round(size(T1,3)/2)))
caxis([0,5])
colormap('gray')
title('True T1')
colorbar
figure
imagesc(T1pred(:,:,round(size(T1pred,3)/2)))
caxis([0,5])
colormap('gray')
title('Predicted T1')
colorbar

figure
imagesc(T2(:,:,round(size(T2,3)/2)))
caxis([0,.7])
colormap('gray')
title('True T2')
colorbar
figure
imagesc(T2pred(:,:,round(size(T2pred,3)/2)))
caxis([0,.7])
colormap('gray')
title('Predicted T2')
colorbar

figure
imagesc(M0(:,:,round(size(M0,3)/2))-M0pred(:,:,round(size(M0,3)/2)))
colormap('gray')
title('M0 Error')
colorbar
figure
imagesc(T1(:,:,round(size(T1,3)/2))-T1pred(:,:,round(size(T1,3)/2)))
colormap('gray')
title('T1 Error')
colorbar
figure
imagesc(T2(:,:,round(size(T2,3)/2))-T2pred(:,:,round(size(T2,3)/2)))
colormap('gray')
title('T2 Error')
colorbar
