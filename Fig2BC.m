clear;close all;clc;FS=22;
%pkg load image

% Co-ordinates and dispersions for the kernel
y=-81:80;
x=-86:85;
[X,Y]=meshgrid(x,y);

D=2;
t=1:30;
sigma2=4*D*t;
 
for j=1:length(t);
 P{j}=exp(-(X.^2+Y.^2)/(4*D*t(j)))/(pi*4*D*t(j));
 end


%%
% Затравочные иозображения с локализацией в 1 пиксель (разрешение любого фото)
% для Гаусса и О'Шонесси-Прокаччиа
im0G=exp(-(X.^2+Y.^2))/(pi); 
sigma=linspace(1,30,200);
% Sequential image processing
for k=1:length(t);
  % Test of correlations
for j=1:length(sigma);
   g=exp(-(X.^2+Y.^2)/(sigma(j)^2))/(pi*sigma(j)^2); 
   Fkernel=fft2(fftshift(g));
% Для нормальной диффузии   
   imgGn=ifft2(fft2(im0G).*Fkernel);
   imGn=real(double(imgGn));
   cn{k}(j)=corr2(P{k},imGn);   
end
[maxcn(k),mjn(k)]=max(cn{k});
% [maxc(k),mj(k)]=max(c{k});
end  



figure
subplot(1,2,1)
for k=1:5:length(t);
   plot(sigma,cn{k}*100,'LineWidth',1.5);   
   hold on
end
text(3,15,'(A)','FontSize',FS)
ylim([0 105])
xlim([0 max(sigma)])
xlabel('\sigma, px')
ylabel('C(\sigma_j), %')
set(gca,'FontSize',FS)

subplot(1,2,2)
%plot(t,sigma(mjn).^2,'o','LineWidth',1)
plot(t,sigma(mjn).^2,'.-','LineWidth',1.5)
hold on
%plot(t,sigma2,'+','LineWidth',1)
plot(t,sigma2,'o','LineWidth',1,'MarkerSize',12)
text(25,40,'(B)','FontSize',FS)
xlim([0 30])
xlabel('t, a.u.')
ylabel('\sigma^2(t), px.^2')
set(gca,'FontSize',FS)
print -dpdf CorrMSDpoint '-S1366,384'
