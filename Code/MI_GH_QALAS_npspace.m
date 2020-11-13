%% Measurement model
% clear all
% close all
% format shortg

function [signal_lib,wn_t1_lib,wn_t2_lib,Msize,parspace,pslabels]=MI_GH_QALAS_npspace(nparspace)

% labelfilename = '/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/ICBM_grey_white_csf.nii.gz'; % population tissue segmentation
% segsize=2;
% switch segsize
%     case 2
%         tissuelabel = load_untouch_nii(labelfilename);
%         segmap=int32(tissuelabel.img(15:165,20:200,92));
%     case 2.5
%         system(sprintf('/opt/apps/itksnap/c3d-1.1.0-Linux-x86_64/bin/c3d %s -dilate 3 8x8x8vox -interpolation NearestNeighbor -resample 36x43x36 -o resampleimg.nii.gz',labelfilename));
%         tissuelabel = load_untouch_nii('resampleimg.nii.gz');
%         system('rm resampleimg.nii.gz');
%         segmap = int32(tissuelabel.img);
%     otherwise
%         tissuelabel = load_untouch_nii(labelfilename);
%         segmap = int32(tissuelabel.img);
% end

M0varflag = 1;
T1varflag = 1;
T2varflag = 1;
%% signal model parameters
%%          GM       WM        CSF     Tumor
T1mean = [1200   , 900     ,  4000    , 250  ]./1000; % s
if T1varflag~=0
    T1stdd = [100    ,  100    ,   200    ,  50  ]./1000; % s
else
    T1stdd = [0,0,0,0];
end
T2mean = [100    ,   80    ,  1000    , 250  ]./1000; % s
if T2varflag~=0
    T2stdd = [  5    ,    4    ,    50    ,  50  ]./1000; % s
else
    T2stdd = [0,0,0,0];
end
M0mean = [  0.9  ,    0.9  ,     1.0  ,   1.2];       % relative intensity
if M0varflag~=0
    M0stdd = [   .05 ,     .05 ,      .05 ,    .4];       % relative intensity
else
    M0stdd = [0,0,0,0];
end

TR = 0.005;              % s
TE_T2prep = 0.100;        % s
Tacq = 0.500;             % s
nacq = 5;
flipAngle = 4;          % deg
TDinv=0.100;            % s

switch nparspace
    case 1
        TDpT2 = [0:20:3000]./1000;
        TDpT1 = [0:20:3000]./1000;
    case 2
        %         TD = repmat(10:300:3000,[nacq-1,1])./1000;
        %         TD1 = [800:30:1100]./1000;
        %         TD2 = [500:20:700]./1000;
        %         TD3 = [500:20:700]./1000;
        %         TD4 = [600:20:800]./1000;
        %         TD5 = [800:40:1200]./1000;
        flipAngle = 1:.5:10;      % degrees
        TD = [0.03;%[900:40:1300;
                0:10:100;
                0:10:100;
                0:10:100;
                0:10:100]./1000;
        parspace{1}=flipAngle; parspace{2}=TD; pslabels={'flipAngle','TD1','TD2','TD3','TD4','TD5'};
        % TD = repmat(10:20:200,[nacq-1,1])./1000; % s
    case 3
        % flipAngle = 5:1:15;
        flipAngle = 1:.5:10;      % degrees
        TD1 = [0:10:200]./1000;   % s
        parspace{1}=flipAngle; parspace{2}=TD1; pslabels={'flipAngle','TD1'};
        %     case 4
        %         flipAngle = 1:.5:10;
        %         TD = [0:100:1000]./1000;
end

NumQP=5;
for labelindex=1:3
    %% Quadrature Points
    switch labelindex
        case 1
            [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,T1mean(1),T1stdd(1));
            [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,T2mean(1),T2stdd(1));
        case 2
            [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,T1mean(2),T1stdd(2));
            [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,T2mean(2),T2stdd(2));
        case 3
            [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,T1mean(3),T1stdd(3));
            [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,T2mean(3),T2stdd(3));
    end
    % [x_t1,xn_t1,xm_t1,w_t1,wn_t1]=GaussHermiteNDGauss(NumQP,[mean(T1mean(1:2)),T1mean(3)],[mean(T1stdd(1:2)),T1stdd(3)]);
    % [x_t2,xn_t2,xm_t2,w_t2,wn_t2]=GaussHermiteNDGauss(NumQP,[mean(T2mean(1:2)),T2mean(3)],[mean(T2stdd(1:2)),T2stdd(3)]);
    switch nparspace
        case 1
            [p1,p2,x1,x2]=ndgrid(TDpT2,TDpT1,xm_t1,xm_t2);
            p1=p1(:);
            p2=p2(:);
            x1=x1(:);
            x2=x2(:);
            lqp=length(p1);
            
            parfor qp=1:lqp
                disp(sprintf('Model eval: %d of %d',qp,lqp))
                dt = [0,0,Tacq,p1(qp),0,TDinv,Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp)];
                [~,Mmodel(:,qp)]=qalas1p(M0mean(1),M0mean(1),x1(qp),x2(qp),TR,TE_T2prep,p1(qp),nacq,dt);
            end
            
            signal_lib(:,labelindex)=Mmodel(:);
            wn_t1_lib(:,labelindex)=wn_t1;
            wn_t2_lib(:,labelindex)=wn_t2;
            
        case 2
            [p1,p2,p3,p4,p5,p6,x1,x2]=ndgrid(flipAngle,TD(1,:),TD(2,:),TD(3,:),TD(4,:),TD(5,:),xm_t1,xm_t2);
            p1=p1(:);
            p2=p2(:);
            p3=p3(:);
            p4=p4(:);
            p5=p5(:);
            p6=p6(:);
            x1=x1(:);
            x2=x2(:);
            lqp=length(p1);
            
            parfor qp=1:lqp
                disp(sprintf('Model eval: %d of %d',qp,lqp))
                dt = [0,0,Tacq,p2(qp),0,TDinv,Tacq,p3(qp),Tacq,p4(qp),Tacq,p5(qp),Tacq,p6(qp)];
                [~,Mmodel(:,qp)]=qalas1p(M0mean(1),M0mean(1),x1(qp),x2(qp),TR,TE_T2prep,p1(qp),nacq,dt);
            end
            
            signal_lib(:,labelindex)=Mmodel(:);
            wn_t1_lib(:,labelindex)=wn_t1;
            wn_t2_lib(:,labelindex)=wn_t2;
            
        case 3
            [p1,p2,x1,x2]=ndgrid(flipAngle,TD1,xm_t1,xm_t2);
            p1=p1(:);
            p2=p2(:);
            x1=x1(:);
            x2=x2(:);
            lqp=length(p1);
            
            parfor qp=1:lqp
                disp(sprintf('Model eval: %d of %d',qp,lqp))
                dt = [0,0,Tacq,1,0,TDinv,Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp)];
                [~,Mmodel(:,qp)]=qalas1p(M0mean(1),M0mean(1),x1(qp),x2(qp),TR,TE_T2prep,p1(qp),nacq,dt);
            end
            
            signal_lib(:,labelindex)=Mmodel(:);
            wn_t1_lib(:,labelindex)=wn_t1;
            wn_t2_lib(:,labelindex)=wn_t2;
            
            %         case 4
            %             [p1,p2,x1,x2]=ndgrid(flipAngle,TD,xm_t1,xm_t2);
            %             p1=p1(:);
            %             p2=p2(:);
            %             x1=x1(:);
            %             x2=x2(:);
            %             lqp=length(p1);
            %
            %             M0model=0.*(segmap==0)+M0mean(1).*(segmap==1)+M0mean(2).*(segmap==2)+M0mean(3).*(segmap==3);
            %             parfor qp=1:lqp
            %                 T1model
            %                 disp(sprintf('Model eval: %d of %d',qp,lqp))
            %                 dt = [0,0,Tacq,1,0,TDinv,Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp),Tacq,p2(qp)];
            %                 [~,Mmodel(:,:,:,:,qp)]=qalas(M0mean(1),M0mean(1),x1(qp),x2(qp),TR,TE_T2prep,p1(qp),nacq,dt);
            %             end
            %
            %             signal_lib(:,labelindex)=Mmodel(:);
            %             wn_t1_lib(:,labelindex)=wn_t1;
            %             wn_t2_lib(:,labelindex)=wn_t2;
    end
    
end

switch nparspace
    case 1
        Msize=[nacq,length(TDpT2),length(TDpT1),size(wn_t1_lib,1)^2];%size(wn_t1_lib,1)^size(wn_t1_lib,2)*size(wn_t2_lib,1)^size(wn_t2_lib,2)];
        disp('Saving...')
        save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas5acq3tissue_np%d_siglib.mat',nparspace),'signal_lib','wn_t1_lib','wn_t2_lib','Msize','-v7.3');
    case 2
        Msize=[nacq,length(flipAngle),size(TD,2)*ones([1,size(TD,1)]),size(wn_t1_lib,1)^2];%size(wn_t1_lib,1)^size(wn_t1_lib,2)*size(wn_t2_lib,1)^size(wn_t2_lib,2)];
        disp('Saving...')
        save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas5acq3tissue_np%d_siglib.mat',nparspace),'signal_lib','wn_t1_lib','wn_t2_lib','Msize','parspace','pslabels','-v7.3');
    case 3
        Msize=[nacq,length(flipAngle),length(TD1),size(wn_t1_lib,1)^2];
        disp('Saving...')
        save(sprintf('/rsrch1/ip/dmitchell2/github/SyntheticMR/Code/qalas5acq3tissue_np%d_siglib.mat',nparspace),'signal_lib','wn_t1_lib','wn_t2_lib','Msize','parspace','pslabels','-v7.3');
end

end