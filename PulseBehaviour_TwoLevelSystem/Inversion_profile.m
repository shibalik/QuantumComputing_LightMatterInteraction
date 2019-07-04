clc
clear all
close all

n=2;
FWHM=3*10^-1;
u=20*pi*10;

alpha=4*log(2)/(FWHM)^2;
P=zeros(n,n);
P(1,1)=1;
N=1000;
N1=2000;
t=linspace(0,0.5,N);
dt=t(2)-t(1);
e=zeros(N);
E=linspace(-80,80,N1);
p=zeros(N1);
R=zeros(N1);
 Omega_Rabbi=zeros(N);
for j=1:N1
     P(:,:,1);
    for i=2:N
      
    e(i-1)=2.19*exp(-(alpha)*t(i-1)^2);
    e(i-1);
    Omega_Rabbi(i-1)= u*e(i-1)/2;
    H=[E(j),Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
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
   p(j)=P(2,2,i-1);
 % disp(p(j));
    P=zeros(n,n);
    P(1,1,1)=1;    %values for reference
end
p
plot(E,real(p),'linewidth',2);
xlabel('Rabi Frequency','fontSize',14);
ylabel('Ground State Population','fontsize',14);
axis([E(1) E(N1) 0 1.1]);
