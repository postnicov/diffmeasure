clear;close all;clc;

foldername='Lin24';

secv=[50,100,200,400,600];
sigma=linspace(1,100,200);
y2=-100:99;
x2=-100:99;
[X,Y]=meshgrid(x2,y2);

mainfolder=cd(foldername);
for j=1:length(secv);
   No{j}=num2str(secv(j),'%05.f');
   filename=['t=',No{j},'.png'];
   a=imread(filename);
   b1=double(a(:,:,1));
   b2=double(a(:,:,2));
   b3=double(a(:,:,3));
   z=b1;%+b2/255+b3/(255*255);
   z(510,491)=2*z(510,490)-z(510,489);
   z(510,492)=2*z(510,490)-z(510,489);    
   P{j}=z(510,:);
    P2=double(a(511-100:511+99,491-99:491+100,1));
    P2(100,100)=2*P2(99,100)-P2(98,100);
    P2(100,101)=2*P2(99,101)-P2(98,101);
      % Test of correlations
  for k=1:length(sigma);
     g=exp(-(X.^2+Y.^2)/(sigma(k)^2))/(pi*sigma(k)^2); 
     ge=exp(-sqrt(X.^2+Y.^2)/sigma(k));
     cn(k)=corr2(P2,g);
     cne(k)=corr2(P2,ge);
  end
[maxcn(j),mjn(j)]=max(cn);
[maxcne(j),mjne(j)]=max(cne);
end
cd(mainfolder);

subplot(1,2,1)
x=[-490:509];
for j=1:length(secv);
   semilogy(x,P{j},'-','LineWidth',1.5)
   hold on
end

xlim([-60 60])
ylim([9 260])
text(-50,130,'(A)','FontSize',22)
xlabel('x, px')
ylabel('u(x)')
set(gca,'FontSize',22)

subplot(1,2,2);
plot(secv,maxcn*100,'o-','MarkerSize',10)
hold on
plot(secv,maxcne*100,'*-','MarkerSize',10)
xlim([0 650])
ylim([96 101])
text(50,100,'(B)','FontSize',22)
xlabel('t, a.u.','FontSize',22)
ylabel('C(\sigma), %','FontSize',22)
set(gca,'FontSize',22)

%print -dpdf PointGaussExp '-S1366,384'