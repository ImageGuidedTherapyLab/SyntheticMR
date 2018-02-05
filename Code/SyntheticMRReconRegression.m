
addpath('/home/dmitchell412/Documents/MATLAB/phantom3d');
slphan=phantom3d('Shepp-Logan',64);
% slphan=slphan(:,:,round(size(slphan,3)/2));

phanval=unique(slphan);

% "True" M0, T1, & T2 values
M0val=[nan,1,0.9,0.8,1];                 %
T1val=[nan,4.000,1.400,1.000,0.100];   % s
T2val=[nan,0.600,0.100,0.075,0.100];   % s

% "True" M0, T1, & T2 maps
M0=zeros(size(slphan));
T1=zeros(size(slphan));
T2=zeros(size(slphan));
for iii=1:length(phanval)
    M0(slphan==phanval(iii))=M0val(iii);
    T1(slphan==phanval(iii))=T1val(iii);
    T2(slphan==phanval(iii))=T2val(iii);
end

% Acquisition Parameters
nacq=5;
nrepeat=3;
TR=0.0026;                      % s
TE_T2prep=0.1;                  % s
flipAngle=5;                    % degrees
dt=0.1*ones([1,6+2*(nacq-1)]);  % s
dt(end)=3;

% Signal Model
% for iii=1:size(M0,1)
%     for jjj=1:size(M0,2)
%         for kkk=1:size(M0,3)
%             Mtrue(iii,jjj,kkk,:)=qalas(M0(iii,jjj,kkk),T1(iii,jjj,kkk),T2(iii,jjj,kkk),TR,TE_T2prep,flipAngle,nacq,dt);
%         end
%     end
% end

Mtrue=qalas(M0,M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
for lll=2:nrepeat
    Mtrue(:,:,:,:,lll)=qalas(Mtrue(:,:,:,14,lll-1),M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
end

% Randomize measurement - 1% standard deviation
% normal case SNR = .3 
% try std=.01
% Mrand=Mtrue.*(randn(size(Mtrue))/100+1);
Mrand=Mtrue;%+randn(size(Mtrue))/100;
Mmeas=cat(4,Mrand(:,:,:,2,:),Mrand(:,:,:,6:2:6+2*(nacq-2),:));

% Repeated optimizations
for lll=1:nrepeat
    
% T2 prediction
% M(2)=M(1).*exp(-TE_T2prep./T2);
T2pred(:,:,:,lll)=TE_T2prep./log(Mrand(:,:,:,1,lll)./Mrand(:,:,:,2,lll));

% Optimization solution for M0 and T1 prediction
for iii=1:size(M0,1)
    for jjj=1:size(M0,2)
        for kkk=1:size(M0,3)
            if sum(isnan(squeeze(Mmeas(iii,jjj,kkk,:,lll))))>0
                M0pred(iii,jjj,kkk,lll)=nan;
                T1pred(iii,jjj,kkk,lll)=nan;
            else
%                 [xm,fval,exitflag,output]=fmincon(@(x) qalasobjfun(x,squeeze(Mmeas(iii,jjj,kkk,:))',TR,TE_T2prep,flipAngle,nacq,dt),[mean(M0val(2:end)),mean(T1val(2:end))],[],[],[],[],[0,0],[1.5,5],[],optimset('TolX',1.e-10,'TolFun',1.e-10));%,'Display','iter-detailed'));
                xm=fminsearch(@(x) qalasobjfun(x,squeeze(Mmeas(iii,jjj,kkk,:,lll))',TR,TE_T2prep,flipAngle,nacq,dt),[mean(M0val(2:end)),mean(T1val(2:end))]);
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
caxis([0,1.5])
colormap('gray')
title('True M0')
colorbar
figure
imagesc(M0pred(:,:,round(size(M0pred,3)/2)))
caxis([0,1.5])
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
caxis([0,1])
colormap('gray')
title('True T2')
colorbar
figure
imagesc(T2pred(:,:,round(size(T2pred,3)/2)))
caxis([0,1])
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

for iii=1:size(Mrand,4)
    figure
    imagesc(Mrand(:,:,round(size(Mrand,3)/2),iii,1))
    caxis([-1,1])
    colormap('gray')
    title(sprintf('M%i',iii))
    colorbar
end
