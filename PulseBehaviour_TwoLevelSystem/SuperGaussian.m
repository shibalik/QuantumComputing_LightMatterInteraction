clc
clear all
close all

n=2;
FWHM=4.5*10^-1;
u=2*pi;


P=zeros(n,n);
P(1,1)=1;
P1=zeros(n,n);
P1(1,1)=1;
N=1000;
t=linspace(0,0.35,N);
dt=t(2)-t(1);
e=zeros(N);
m=zeros(N);
E=linspace(0,10,N);
p=zeros(N);
q=zeros(N);
R=zeros(N);
 Omega_Rabbi=zeros(N);
for j=1:N
     P(:,:,1);
    for i=2:N
        
    m(i-1)=E(j)*exp(-(log(2))*(t(i-1)/FWHM)^2);  
    e(i-1)=E(j)*sech(2*log(2+sqrt(3))*t(i-1)/FWHM);
    e(i-1);
    Omega_Rabbi(i-1)= u*e(i-1)/2;
    Omega_Rabbi1(i-1)= u*m(i-1)/2;
    
    H=[0,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
    H1=[0,Omega_Rabbi1(i-1);Omega_Rabbi1(i-1),0];
    K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
    K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
    K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
    K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1 )+K3))*dt;
    
    L1=1i*(P1(:,:,i-1)*H1-H*P1(:,:,i-1))*dt;
    L2=1i*((P1(:,:,i-1)+0.5.*L1)*H1-H1*(P1(:,:,i-1)+0.5.*L1))*dt;
    L3=1i*((P1(:,:,i-1)+0.5.*L2)*H1-H1*(P1(:,:,i-1)+0.5.*L2))*dt;
    L4=1i*((P1(:,:,i-1)+L3)*H1-H*(P1(:,:,i-1 )+L3))*dt;
    
    P(:,:,i) = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
    P1(:,:,i) = P1(:,:,i-1)+((L1+2.*(L2+L3)+L4)/6);
    
       if i==1000
           R(j)=Omega_Rabbi(i-1);
           break
       end
            
   
    end
   p(j)=P(1,1,i);
   q(j)=P1(1,1,i);
 % disp(p(j));
    P=zeros(n,n);
    P(1,1,1)=1;%values for reference
    
     P1=zeros(n,n);
    P1(1,1,1)=1;
end
p;
plot(E,p,E,q,'linewidth',2);
xlabel('Rabi Frequency','fontSize',14);
ylabel('Ground State Population','fontsize',14);
legend({'Super Gaussian','Gaussian'},'Location','northeast')
axis([0 E(N) 0 1.1])
