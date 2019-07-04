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
E=linspace(0,20,2*N);
p=zeros(2*N);
R=zeros(2*N);
 Omega_Rabbi=zeros(N);
 B=linspace(-80*1e12,80*1e12,2*N);
for k=1:2*N
  for j=1:2*N
     P(:,:,1);
     for i=2:N
      
    e(i-1)=E(j)*exp(-(alpha)*t(i-1)^2);
    e(i-1);
    Omega_Rabbi(i-1)= u*e(i-1)/2;
    H=[B(k),Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
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
 end
[X,Y]=mesh(E,B);
figure
mesh(X,Y,p);