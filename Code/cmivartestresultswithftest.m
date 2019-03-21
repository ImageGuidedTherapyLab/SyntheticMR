load results/mi_globalsearch_10patients.mat mi;
size(mi)
[y,i]=max(-mi,[],2);
delayit=0:.02:1;
maxpoints=[delayit(i)',y];

for iii=1:10
plot(delayit,-mi(iii,:));
end
plot(maxpoints(:,1),maxpoints(:,2),'rx');
legend('1','2','3','4','5','6','7','8','9','10');

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





