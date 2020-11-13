
function [] = plotcorr(z,axthresh,plotthresh)

% axthresh=[.85,1.05,.8,1.4,0.070,.105]

if isempty(axthresh)
    axthresh=[0,13,0,4,0,1.5,0,13,0,4,0,1.5,0,13,0,4,0,1.5];
end
if isempty(plotthresh)
    plotthresh=axthresh;
end

nparam=9;
NumTicks = 2;
fontsz=6;
labels = {'M0_{GM}','T1_{GM}','T2_{GM}','M0_{WM}','T1_{WM}','T2_{WM}','M0_{CSF}','T1_{CSF}','T2_{CSF}'};
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nparam
    for j = 1:nparam
        if i == j
            subplot(nparam,nparam,nparam*(i-1)+j);

            tmp=z(:,i);
            histogram(tmp(and(tmp>=plotthresh(2*i-1),tmp<=plotthresh(2*i))),20);
            yl=get(gca,'YLim');
            axis([axthresh(2*i-1),axthresh(2*i),yl(1),yl(2)]);
%             if (i==1)||(i==4)||(i==7)
%                 tmp=z(:,i);
%                 histogram(tmp(and(tmp>=plotthresh(1),tmp<=plotthresh(2))),20);
%                 yl=get(gca,'YLim');
%                 axis([axthresh(1),axthresh(2),yl(1),yl(2)]);
%             elseif (i==2)||(i==5)||(i==8)
%                 tmp=z(:,i);
%                 histogram(tmp(and(tmp>=plotthresh(3),tmp<=plotthresh(4))),20);
%                 yl=get(gca,'YLim');
%                 axis([axthresh(3),axthresh(4),yl(1),yl(2)]);
%             elseif (i==3)||(i==6)||(i==9)
%                 tmp=z(:,i);
%                 histogram(tmp(and(tmp>=plotthresh(5),tmp<=plotthresh(6))),20);
%                 yl=get(gca,'YLim');
%                 axis([axthresh(5),axthresh(6),yl(1),yl(2)]);
%             end
            
            xl = get(gca,'XLim');
            set(gca,'XTick',linspace(xl(1),xl(2),NumTicks),'FontSize',fontsz);
            yl=get(gca,'YLim');
            set(gca,'YTick',linspace(yl(1),yl(2),NumTicks),'FontSize',fontsz);
        end
        
        if i < j
            subplot(nparam, nparam, nparam*(i-1)+j)%count)
            
            tmpi=z(:,i); tmpj=z(:,j);
            indi=and(tmpi>=plotthresh(2*i-1),tmpi<=plotthresh(2*i));
            indj=and(tmpj>=plotthresh(2*j-1),tmpj<=plotthresh(2*j));
%             if (i==1)||(i==4)||(i==7)
%                 indi=and(tmpi>=plotthresh(1),tmpi<=plotthresh(2));
%             elseif (i==2)||(i==5)||(i==8)
%                 indi=and(tmpi>=plotthresh(3),tmpi<=plotthresh(4));
%             elseif (i==3)||(i==6)||(i==9)
%                 indi=and(tmpi>=plotthresh(5),tmpi<=plotthresh(6));
%             end
%             if (j==1)||(j==4)||(j==7)
%                 indj=and(tmpj>=plotthresh(1),tmpj<=plotthresh(2));
%             elseif (j==2)||(j==5)||(j==8)
%                 indj=and(tmpj>=plotthresh(3),tmpj<=plotthresh(4));
%             elseif (j==3)||(j==6)||(j==9)
%                 indj=and(tmpj>=plotthresh(5),tmpj<=plotthresh(6));
%             end
            
            histogram2(tmpi(and(indi,indj)), tmpj(and(indi,indj)), [20,20], 'FaceColor','flat');
            view(30,60);
            
            yl=ylim; zl=zlim;
            axis([axthresh(2*i-1),axthresh(2*i),yl(1),yl(2),zl(1),zl(2)]);
%             if (i==1)||(i==4)||(i==7)
%                 axis([axthresh(1),axthresh(2),yl(1),yl(2),zl(1),zl(2)]);
%             elseif (i==2)||(i==5)||(i==8)
%                 axis([axthresh(3),axthresh(4),yl(1),yl(2),zl(1),zl(2)]);
%             elseif (i==3)||(i==6)||(i==9)
%                 axis([axthresh(5),axthresh(6),yl(1),yl(2),zl(1),zl(2)]);
%             end
            xl=xlim;
            axis([xl(1),xl(2),axthresh(2*j-1),axthresh(2*j),zl(1),zl(2)]);
%             if (j==1)||(j==4)||(j==7)
%                 axis([xl(1),xl(2),axthresh(1),axthresh(2),zl(1),zl(2)]);
%             elseif (j==2)||(j==5)||(j==8)
%                 axis([xl(1),xl(2),axthresh(3),axthresh(4),zl(1),zl(2)]);
%             elseif (j==3)||(j==6)||(j==9)
%                 axis([xl(1),xl(2),axthresh(5),axthresh(6),zl(1),zl(2)]);
%             end
            
            xl = get(gca,'XLim');
            set(gca,'XTick',linspace(xl(1),xl(2),NumTicks),'FontSize',fontsz);
            yl=get(gca,'YLim');
            set(gca,'YTick',linspace(yl(1),yl(2),NumTicks),'FontSize',fontsz);
            zl=get(gca,'ZLim');
            set(gca,'ZTick',linspace(zl(1),zl(2),NumTicks),'FontSize',fontsz);
        end
        
        if i > j
            subplot(nparam, nparam, nparam*(i-1)+j)%count)
                        
            tmpi=z(:,i); tmpj=z(:,j);
            indi=and(tmpi>=plotthresh(2*i-1),tmpi<=plotthresh(2*i));
            indj=and(tmpj>=plotthresh(2*j-1),tmpj<=plotthresh(2*j));
%             if (i==1)||(i==4)||(i==7)
%                 indi=and(tmpi>=plotthresh(1),tmpi<=plotthresh(2));
%             elseif (i==2)||(i==5)||(i==8)
%                 indi=and(tmpi>=plotthresh(3),tmpi<=plotthresh(4));
%             elseif (i==3)||(i==6)||(i==9)
%                 indi=and(tmpi>=plotthresh(5),tmpi<=plotthresh(6));
%             end
%             if (j==1)||(j==4)||(j==7)
%                 indj=and(tmpj>=plotthresh(1),tmpj<=plotthresh(2));
%             elseif (j==2)||(j==5)||(j==8)
%                 indj=and(tmpj>=plotthresh(3),tmpj<=plotthresh(4));
%             elseif (j==3)||(j==6)||(j==9)
%                 indj=and(tmpj>=plotthresh(5),tmpj<=plotthresh(6));
%             end
            
            scatter(tmpj(and(indi,indj)), tmpi(and(indi,indj)),'r.'); 
            
            yl=ylim;
            axis([axthresh(2*j-1),axthresh(2*j),yl(1),yl(2)]);
%             if (j==1)||(j==4)||(j==7)
%                 axis([axthresh(1),axthresh(2),yl(1),yl(2)]);
%             elseif (j==2)||(j==5)||(j==8)
%                 axis([axthresh(3),axthresh(4),yl(1),yl(2)]);
%             elseif (j==3)||(j==6)||(j==9)
%                 axis([axthresh(5),axthresh(6),yl(1),yl(2)]);
%             end
            xl=xlim;
            axis([xl(1),xl(2),axthresh(2*i-1),axthresh(2*i)]);
%             if (i==1)||(i==4)||(i==7)
%                 axis([xl(1),xl(2),axthresh(1),axthresh(2)]);
%             elseif (i==2)||(i==5)||(i==8)
%                 axis([xl(1),xl(2),axthresh(3),axthresh(4)]);
%             elseif (i==3)||(i==6)||(i==9)
%                 axis([xl(1),xl(2),axthresh(5),axthresh(6)]);
%             end
            
            xl = get(gca,'XLim');
            set(gca,'XTick',linspace(xl(1),xl(2),NumTicks),'FontSize',fontsz);
            yl=get(gca,'YLim');
            set(gca,'YTick',linspace(yl(1),yl(2),NumTicks),'FontSize',fontsz);
        end
        
        if j == 1
            ylabel(labels{i},'FontSize',fontsz+6);
        end
        if i == 9
            xlabel(labels{j},'FontSize',fontsz+6);
        end
    end
end
end