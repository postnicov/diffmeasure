clear all;close all;clc;
%  pkg load image
%% Set of image numbers;

Nn=[50:10:320];
n2=find(Nn>140);

ImgWidth=1024;
ImgHight=1248;
%% Load images
mainfolder=cd('images');
im00=imread('Sec042.bmp');
for j=1:length(Nn);
    No{j}=num2str(Nn(j),'%03.f');
    filename=['Sec',No{j},'.bmp'];
    img{j}=imread(filename);
end
cd(mainfolder);

% Initial image
im0=im00(630:791,710:881);
mask=im0>=61;

masks=mask;
masks(50:end,:)=0;
masks=1-masks;

im0=double(im0);
mask=double(mask);

%mask(:,1:45)=0;
%mask(1:45,:)=0;

im0=double(im0.*mask);
% Co-ordinates and dispersions for the kernel
y=-81:80;
x=-86:85;
[X,Y]=meshgrid(x,y);
n=1:200;
sigma=0+n/2;
% Sequential image processing
for k=1:length(No);
  % Test of correlations
for j=n;
   ge=exp(-sqrt(X.^2+Y.^2)/sigma(j)); 
   g=exp(-(X.^2+Y.^2)/(sigma(j)^2))/(pi*sigma(j)^2); 
   Fkernel=fft2(g);
   Fkernele=fft2(ge);
   imgG=fftshift(ifft2(fft2(im0).*Fkernel));
   imgGe=fftshift(ifft2(fft2(im0).*Fkernele));
   imG=real(double(imgG));
   imGe=real(double(imgGe));
   c{k}(j)=corr2(img{k}(630:791,710:881),imG);
   ce{k}(j)=corr2(img{k}(630:791,710:881),imGe);
end
[maxc(k),mj(k)]=max(c{k});
[maxce(k),mje(k)]=max(ce{k});
   g=exp(-(X.^2+Y.^2)/(sigma(mj(k))^2))/(pi*sigma(mj(k))^2); 
   ge=exp(-sqrt(X.^2+Y.^2)/sigma(mje(k))); 
   Fkernel=fft2(g);
   Fkernele=fft2(ge);
   imout{k}=real(fftshift(ifft2(fft2(im0).*Fkernel)));
   imoute{k}=real(fftshift(ifft2(fft2(im0).*Fkernele)));
end  

% Figures

figure(1)
subplot(1,2,1)
imagesc(im00(630:791,710:881))
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(A)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(1,2,2)
imagesc(im0)
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(B)','FontSize',22)
set(gca,'FontSize',22)
axis square
print -dpdf initmasc '-S1366,768'

figure
subplot(1,2,1)
plot(Nn,maxc*100,'o','LineWidth',1)
hold on
plot(Nn,maxce*100,'*','LineWidth',1)
xlabel('t, s','FontSize',22)
ylabel('C(\sigma), %','FontSize',22)
set(gca,'FontSize',22)
subplot(1,2,2)
plot(Nn,sigma(mj).^2,'o','LineWidth',1)
hold on
plot(Nn,sigma(mje).^2,'*','LineWidth',1)

[p,S]=polyfit(Nn(n2),sigma(mj(n2)).^2,1);
plot(Nn(n2),polyval(p,Nn(n2),S),'-','color','black')

[pe,Se]=polyfit(Nn(n2),sigma(mje(n2)).^2,2);
plot([Nn(n2) 330],polyval(pe,[Nn(n2) 330],Se),'--','color','black')

D=p/4
sqrt(diag(S.C)/S.df)*S.normr

xlabel('t, s','FontSize',22)
ylabel('\sigma^2, \mu m^2','FontSize',22)
set(gca,'FontSize',22)


print -dpdf Scorr  '-S1366,768'

figure
subplot(2,3,1)
imagesc(img{1}(630:791,710:881))
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(A)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,3,2)
imagesc(img{8}(630:791,710:881))
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(B)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,3,3)
imagesc(img{16}(630:791,710:881))
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(C)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,3,4)
imagesc(imoute{1})
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(D)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,3,5)
imagesc(imoute{8})
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(E)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,3,6)
imagesc(imout{16})
xlabel('\mu m','FontSize',22);
ylabel('\mu m','FontSize',22);
title('(F)','FontSize',22)
set(gca,'FontSize',22)
axis square
print -dpdf Simcomp '-S1366,768'


%% MATLAB's determinaiton of the coefficient's confidence interval 
%fitType = fittype('a*x + b');
%f = fit(Nn(n2)',sigma(mj(n2))'.^2,fitType);
%uncertainty = confint(f,0.95)
