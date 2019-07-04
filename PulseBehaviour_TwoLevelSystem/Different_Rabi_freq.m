clc
clear all;
close all;

prompt= 'input no. of levels in your system ';
n=input(prompt);
prompt= 'input FWHM of your pulse ';
FWHM=input(prompt);
alpha=log(2)/(FWHM)^2;

prompt= 'input dipole moment of your system  ';
u=input(prompt);
P=zeros(n,n);
P(1,1)=1;
N=1000;
t=linspace(0,0.25e-12,N);
dt=t(2)-t(1);
E=linspace(0,200,1000);
disp(E(2))

q=zeros(1000);
p=zeros(1000);

for j=1:1000
    e=zeros(N);
    Omega_Rabbi=zeros(N);
   for i=2:N
       e(i-1)=E(j)*exp(-(alpha)*t(i-1)^2);
    e(i-1);
    Omega_Rabbi(i-1)= u*e(i-1)/2;
    H=[0,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
    %H1=[0,10*3.14*10^12;10*3.14*10^12,0];
    K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
    K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
    K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
    K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1)+K3))*dt;
    
   

    P(:,:,i) = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
    %%if e(i)-e(i-1)<0.001*e(i-1)
      %  break;
   % end
    
    
   
   end
    p(j)=P(2,2,N);
    p(j);
    q(j)=P(1,1,N);
end
plot(E,p,'linewidth',2);
p(10);
%plot(t,q,'k','linewidth',2)

xlabel('Amplitude','fontSize',14);
ylabel('Population','fontsize',14);
axis([0 E(1000) 0 1.1])

%axis([0 E(N) 0 1.1])
%legend({'Ground state population','excited state population'},'Location','North east')
fh = figure(1);
set(fh, 'color', 'white');
