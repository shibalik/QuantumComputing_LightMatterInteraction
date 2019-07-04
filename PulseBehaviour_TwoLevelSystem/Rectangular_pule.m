clc
clear all
close all

n=2;
FWHM=4.5*10^-14;
u=20*pi*10^12;

alpha=log(2)/(FWHM)^2;
P=zeros(n,n);
P(1,1)=1;
N=1000;
t=linspace(0,0.5e-12,N);
dt=t(2)-t(1);
e=zeros(N);
E=linspace(0,10,N);
p=zeros(N);
R=zeros(N);
 Omega_Rabbi=zeros(N);
for j=1:N
     P(:,:,1);
    for i=2:N
      e(i-1)=E(j);
     
    Omega_Rabbi(i-1)= u*e(i-1)/2;
    H=[0,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
    %H1=[0,10*3.14*10^12;10*3.14*10^12,0];
    K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
    K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
    K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
    K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1 )+K3))*dt;
    
    P(:,:,i) = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
       if i==1000
           R(j)=Omega_Rabbi(i-1);
           break
       end
            
   
    end
   p(j)=P(1,1,i);
 % disp(p(j));
    P=zeros(n,n);
    P(1,1,1)=1;    %values for reference
end
p;
plot(E,p,'linewidth',2);
xlabel('Rabi Frequency','fontSize',14);
ylabel('Ground State Population','fontsize',14);
axis([0 E(N) 0 1.1])
