load results/mi_globalsearch_10patients.mat mi;
size(mi)
[y,i]=max(-mi,[],2);
delayit=0:.02:1;
maxpoints=[delayit(i)',y];

figure;hold on;
for iii=1:10
plot(delayit,-mi(iii,:),'LineWidth',1.5);
end
plot(maxpoints(:,1),maxpoints(:,2),'rx','LineWidth',3.5);
legend('1','2','3','4','5','6','7','8','9','10');

figure; hold on;
for iii=1:10
    plot(delayit,-mi(iii,:)+mi(iii,26),'LineWidth',1.25);
end
[y,i]=max(-mi+repmat(mi(:,26),[1,51]),[],2);
delayit=0:.02:1;
maxpoints=[delayit(i)',y];
plot(maxpoints(:,1),maxpoints(:,2),'rx','LineWidth',3.5);
legend('1','2','3','4','5','6','7','8','9','10','Optimum');
xlabel('Delay Times (s)'); ylabel('Conditional Mutual Information');
saveas(gcf,'Figures/cmi_glob_10pat','png');

figure; hold on;
for iii=1:10
    plot(delayit,-mi(iii,:)+mi(iii,26),'LineWidth',1.25);
end
[y,i]=max(-mi+repmat(mi(:,26),[1,51]),[],2);
delayit=0:.02:1;
maxpoints=[delayit(i)',y];
plot(maxpoints(:,1),maxpoints(:,2),'rx','LineWidth',3.5);
legend('1','2','3','4','5','6','7','8','9','10','Optimum');
axis([0,1,-500,500]);
xlabel('Delay Times (s)'); ylabel('Conditional Mutual Information');
saveas(gcf,'Figures/cmi_glob_10pat_zoom','png');

load cmivartestresultswithftest.mat;
ps31=std(post_samples{3,1}(:,1:9));
ps32=std(post_samples{3,2}(:,1:9));
(ps32-ps31)./ps32
for iii=1:9
[h(iii),p(iii)]=vartest2(post_samples{3,1}(:,iii),post_samples{3,2}(:,iii));
end
p
p(1)
ps51=std(post_samples{5,1}(:,1:9));
ps52=std(post_samples{5,2}(:,1:9));
(ps52-ps51)./ps52
for iii=1:9
[h(iii),p(iii)]=vartest2(post_samples{5,1}(:,iii),post_samples{5,2}(:,iii));
end
p
for iii=1:9
[h(iii),p(iii)]=vartestn(post_samples{5,1}(:,iii),post_samples{5,2}(:,iii));
end


for iii=[1,3,5,7]
    for jjj=1:2
        psstd1(:,iii,jjj)=std(post_samples{iii,jjj},[],1);
    end
end
psstd1(:,2:2:end,:)=[];
relstdshift1=(psstd1(:,:,1)-psstd1(:,:,2))./psstd1(:,:,2);
savestd1=relstdshift1(1:9,:)';



load conditionalMIresults2.mat;


for iii=1:10
    for jjj=1:2
        psstd(:,iii,jjj)=std(post_samples{iii,jjj},[],1);
    end
end
relstdshift=(psstd(:,:,1)-psstd(:,:,2))./psstd(:,:,2);
savestd=[relstdshift(1:9,:)',TDin];





