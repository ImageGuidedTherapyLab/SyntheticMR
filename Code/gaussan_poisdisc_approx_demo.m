
%% 
close all; clear all;

% DIMENSIONS
xdim=150;
ydim=200;
accel=1;
variance=1050;

% DISCRETE 2-D VARIABLE DENSITY POISSON DISC SAMPLING
pdisc=zeros([xdim,ydim]);
nsample=1000;
parfor iii=1:nsample
    tmpdisc=bart(sprintf('poisson -Y %i -Z %i -y %i -z %i -v -s %i',xdim,ydim,accel,accel,iii));
    pdisc=pdisc+squeeze(tmpdisc);
end
pdisc=pdisc/nsample;
pdiscn=pdisc/sum(pdisc(:));

% DISCRETE 2-D GAUSSIAN DISTRIBUTION

x=-ceil(xdim/2-1):floor(xdim/2);
y=-ceil(ydim/2-1):floor(ydim/2);

gaussx=exp(-variance).*besseli(x,variance);
gaussy=exp(-variance).*besseli(y,variance);
gauss2d=kron(gaussx',gaussy);

% SAMPLED CONTINUOUS 2-D GAUSSIAN DISTRIBUTION

[xg,yg]=ndgrid(x,y);
sig=diag([variance,variance]);
mu=[0,0];
for iii=1:numel(xg)
    gauss2dc(iii)=(det(2*pi*sig))^(-0.5)*exp(-0.5*(mu-[xg(iii),yg(iii)])*inv(sig)*(mu-[xg(iii),yg(iii)])');
end
gauss2dc=reshape(gauss2dc,size(xg));
gauss2dcn=gauss2dc/sum(gauss2dc(:));
figure;plot(squeeze(pdiscn(76,:)));hold on; plot(squeeze(gauss2dcn(76,:)));

% SAMPLED CONTINUOUS 2-D MIXTURE OF GAUSSIAN DISTRIBUTIONS




% FIGURES

cmax=max(gauss2dcn(:));
figure;
% subplot(1,3,1);
imagesc(pdiscn); colorbar; title('Mean Poisson Disc'); caxis([0,cmax]);
figure;
% subplot(1,3,2);
imagesc(gauss2d); colorbar; title('2D Discrete Gaussian'); caxis([0,cmax]);
figure;
imagesc(gauss2dcn); colorbar; title('2D Continuous Gaussian'); caxis([0,cmax]);
figure;
% subplot(1,3,3);
imagesc(pdiscn-gauss2d); colorbar; title('Error - Discrete');
figure;
imagesc(pdiscn-gauss2dcn); colorbar; title('Error - Continuous');

figure;
imagesc(pdiscn-gauss2dcn); colorbar; title('Error - Continuous'); caxis([0,cmax]);