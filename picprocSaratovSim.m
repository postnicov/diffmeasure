clear;close all;clc;
%pkg load image

%foldername='new_data_10_30_lin';

foldername='new_data_10_300_lin';

mainfolder=cd(foldername);
% Read vascular structure
Iv=imread('vessels.png');
Iv=rgb2gray(Iv);
mask=Iv>250;
% Read leakage source
I0=imread('source.png');
I0=rgb2gray(I0);
cd(mainfolder);

% imagesc(I0)
% 260:360,250:300

%% Set of image numbers;

% secv=[6:2:40,50:10:120];
secv=[30:10:900];

n2=find(secv>400);

%% Load images
mainfolder=cd(foldername);
for j=1:length(secv);
    No{j}=num2str(secv(j),'%05.f');
    filename=['t=',No{j},'.png'];
    [img{j},map]=imread(filename);
    img{j}=ind2rgb(img{j},map);
    img{j}=rgb2gray(img{j});
    img{j}=double(img{j});
end
cd(mainfolder);
%
%% Initial image
im0=double(I0(260:360,250:300));
%% Co-ordinates and dispersions for the kernel
y=-50:50;
x=-25:25;
[X,Y]=meshgrid(x,y);
n=1:150;
sigma=n/2;
% Sequential image processing
for k=1:length(No);
  % Test of correlations
for j=n;
%   g=exp(-(X.^2+Y.^2)/(sigma(j)^2))/(pi*sigma(j)^2); 
    g=exp(-sqrt(X.^2+Y.^2)/sigma(j));   
   Fkernel=fft2(g);
   imgG=fftshift(ifft2(fft2(im0).*Fkernel));
   imG=real(imgG).*double(1-mask(260:360,250:300));
   c{k}(j)=corr2(img{k}(260:360,250:300).*double(1-mask(260:360,250:300)),imG);
end
[maxc(k),mj(k)]=max(c{k});
%   g=exp(-(X.^2+Y.^2)/(sigma(mj(k))^2))/(pi*sigma(mj(k))^2); 
      g=exp(-sqrt(X.^2+Y.^2)/(sigma(mj(k)))); 
   Fkernel=fft2(g);
   imout{k}=real(fftshift(ifft2(fft2(im0).*Fkernel)));
end  

%% Figures
figure
subplot(1,2,1)
plot(secv,maxc*100,'o','LineWidth',1)
%  secvE=secv;
%  maxcE=maxc*100;
%  sigmaE2=sigma(mj).^2;
%  save('expout.mat','secvE','maxcE','sigmaE2');
hold on
xlabel('t, a.u.','FontSize',22)
ylabel('C(\sigma), %','FontSize',22)
title('(A)','FontSize',22)
set(gca,'FontSize',22)
subplot(1,2,2)
plot(secv,sigma(mj).^2,'o')
p2=polyfit(secv(n2),sigma(mj(n2)).^2,1)
hold on
%plot(secv(3:6),polyval(p1,secv(3:6)),'--','color','red')
xlabel('t, a.u.','FontSize',22)
ylabel('\sigma^2','FontSize',22)
title('(B)','FontSize',22)
set(gca,'FontSize',22)

load expout
subplot(1,2,1)
plot(secvE,maxcE,'*','LineWidth',1)
subplot(1,2,2)
plot(secvE,sigmaE2,'*','LineWidth',1)
plot([142:900],polyval(p2,[142:900]),'-','color','black','LineWidth',1.5)

print -dpng ScorrA300 '-S1366,384'

%
figure
subplot(2,4,1)
%imagesc(img{1}(260:360,250:300).*(1-mask(260:360,250:300)))
imagesc(mask(260:360,250:300))
xlabel('px')
ylabel('px')
title('(A)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,2)
imagesc(img{13}(260:360,250:300).*(1-mask(260:360,250:300)))
xlabel('px')
ylabel('px')
title('(B)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,3)
imagesc(img{33}(260:360,250:300).*(1-mask(260:360,250:300)))
xlabel('px')
ylabel('px')
title('(C)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,4)
imagesc(img{63}(260:360,250:300).*(1-mask(260:360,250:300)))
xlabel('px')
ylabel('px')
title('(D)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,5)
%imagesc(imout{1}.*(1-mask(260:360,250:300)))
imagesc(I0(260:360,250:300))
title('(E)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,6)
imagesc(imout{13}.*(1-mask(260:360,250:300)))
xlabel('px')
ylabel('px')
title('(F)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,7)
imagesc(imout{33}.*(1-mask(260:360,250:300)))
xlabel('px')
ylabel('px')
title('(G)','FontSize',22)
set(gca,'FontSize',22)
axis square
subplot(2,4,8)
imagesc(imout{63}.*(1-mask(260:360,250:300)))
xlabel('px')
ylabel('px')
title('(H)','FontSize',22)
set(gca,'FontSize',22)
axis square
print -dpng SimcompA300E '-S1366,768'
