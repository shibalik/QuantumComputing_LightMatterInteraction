clc
clear all;
close all;


prompt= 'input no. of levels in your system ';
n=input(prompt);
prompt= 'input FWHM of your pulse ';
FWHM=input(prompt);
alpha=log(4)/(FWHM)^4;

prompt= 'input dipole moment of your system  ';
u=input(prompt);

P=zeros(n,n);
P(1,1)=1;

P1=zeros(n,n);
P1(1,1)=1;

N=1000;
t=linspace(0,0.5e-12,N);
dt=t(2)-t(1);
e=zeros(N);
Omega_Rabbi=zeros(N);
for i=2:N
    e(i-1)=18*sech(2*log(2+sqrt(3))*t(i-1)/FWHM);
    
    Omega_Rabbi(i-1)= u*e(i-1)/2;
    H=[0,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
    %H1=[0,10*3.14*10^12;10*3.14*10^12,0];
    K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
    K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
    K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
    K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1)+K3))*dt;
    
   

    P(:,:,i) = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
            
  
end
m=P(1,1,N);
q=zeros(N);
p=zeros(N);
for j=1:N
    p(j)=P(2,2,j);
    q(j)=P(1,1,j);
    
end
figure(1)
P(1,1,N);
plot(t,p,t,q,'linewidth',2);

%plot(t,q,'k','linewidth',2)

xlabel('time','fontSize',14);
ylabel('Population','fontsize',14);
axis([0 t(N) 0 1.1])

axis([0 t(N) 0 1.1])
legend({'Ground state population','excited state population'},'Location','east')
fh = figure(1);
set(fh, 'color', 'white');
