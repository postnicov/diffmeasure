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
 P0{j}=exp(-(X.^2+Y.^2)/(4*D*t(j)))/(pi*4*D*t(j));
 P{j}=exp(-(X.^2+Y.^2)/(4*D*t(j)))/(pi*4*D*t(j))+exp(-((X-10).^2+(Y-10).^2)/(4*D*t(j)))/(pi*4*D*t(j));
end

subplot(2,4,1)
imagesc(P0{1});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,2)
imagesc(P0{4});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,3)
imagesc(P0{10});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,4)
imagesc(P0{30});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,5)
imagesc(P{1});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,6)
imagesc(P{4});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,7)
imagesc(P{10});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
subplot(2,4,8)
imagesc(P{30});
colorbar
axis square
xlabel('px')
ylabel('px')
set(gca,'FontSize',FS)
print -dpdf fig12points '-S1366,600'