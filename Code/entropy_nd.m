function H = entropy_nd(inputdists,binwidths)

% calculate histogram counts
% nbin=256;
p=cell(size(inputdists,2),1);
for iii=1:size(inputdists,2)
%     edges=linspace(median(inputdists(:,iii))-binwidths(iii)*nbin/2,median(inputdists(:,iii))+binwidths(iii)*nbin/2,nbin+1);
%     [p{iii},~]=histcounts(inputdists(:,iii),edges,'Normalization','probability');
    [p{iii},~]=histcounts(inputdists(:,iii),'BinWidth',binwidths(iii),'Normalization','probability');
    p{iii}(p{iii}==0)=[];
end

pk=p{1};
for iii=2:size(inputdists,2)
    pk=kron(pk,p{iii});
end

% entropy
H=-sum(pk.*log(pk));

end