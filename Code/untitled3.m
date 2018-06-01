crashplots=1;
if crashplots==1

for iii=1:size(Mmeassave,1)
    for jjj=1:size(Mmeassave,2)
            clear gmmmplot wmmmplot;
%             TD=[0.5,0.001,0.1,TD4in(iii)];
            TD=tduniform(iii,:);
            dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
            dtplot=cumsum(dt)';
            mmvec=reshape(Mmeassave{iii,jjj},[size(Mmeassave{iii,jjj},1)*size(Mmeassave{iii,jjj},2)*size(Mmeassave{iii,jjj},3),size(Mmeassave{iii,jjj},4)]);
            mmvecfull=reshape(Mmeasfullsave{iii,jjj},[size(Mmeasfullsave{iii,jjj},1)*size(Mmeasfullsave{iii,jjj},2)*size(Mmeasfullsave{iii,jjj},3),size(Mmeasfullsave{iii,jjj},4)]);
            for yyy=1:size(Mmeassave{iii,jjj},4)
                gmmmplot(:,yyy)=mmvec(materialID(:)==1,yyy);
                wmmmplot(:,yyy)=mmvec(materialID(:)==2,yyy);
            end
            
%             titlename=sprintf('Gray Matter, TD4=%2.2f, SS=%2.2f',TD4in(iii),subsampin(jjj));
            titlename=sprintf('GM %2.2f TD=[%2.2f,%2.2f,%2.2f,%2.2f] var=[%2.2f,%2.2f,%2.2f] MI=%5.4e',subsampin(jjj),TD(1),TD(2),TD(3),TD(4),varcell{iii,jjj}(1,1),varcell{iii,jjj}(2,1),varcell{iii,jjj}(3,1),-MIobjfun{iii,jjj});
            figure; hold on;
            for yyy=1:size(gmmmplot,1)
                plot(dtplot([2,6:2:end-1]),squeeze(gmmmplot(yyy,:)),'-o');
            end
            plot(dtplot,squeeze(sind(flipAngle)*Mmeasfullsave{iii,jjj}(16,58,1,:)),'r--s');
            xlabel('Time (s)'); ylabel('M'); title(titlename);
            
%             titlename=sprintf('White Matter, TD4=%2.2f, SS=%2.2f',TD4in(iii),subsampin(jjj));
            titlename=sprintf('GM %2.2f TD=[%2.2f,%2.2f,%2.2f,%2.2f] var=[%2.2f,%2.2f,%2.2f] MI=%5.4e',subsampin(jjj),TD(1),TD(2),TD(3),TD(4),varcell{iii,jjj}(1,2),varcell{iii,jjj}(2,2),varcell{iii,jjj}(3,2),-MIobjfun{iii,jjj});
            figure; hold on;
            for yyy=1:size(wmmmplot,1)
                plot(dtplot([2,6:2:end-1]),squeeze(wmmmplot(yyy,:)),'-o');
            end
            plot(dtplot,squeeze(sind(flipAngle)*Mmeasfullsave{iii,jjj}(49,36,1,:)),'r--s');
            xlabel('Time (s)'); ylabel('M'); title(titlename);
        end
    end
end