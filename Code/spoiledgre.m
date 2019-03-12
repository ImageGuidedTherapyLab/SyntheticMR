function S = spoiledgre(K,H,flipAngle,TR,TE,T1,T2star)
    S = (K.*H.*sind(flipAngle).*(1-exp(-TR./T1)).*exp(-TE./T2star))./(1-exp(-TR./T1).*cosd(flipAngle));
end