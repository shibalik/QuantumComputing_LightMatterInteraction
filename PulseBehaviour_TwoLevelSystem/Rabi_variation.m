clc
clear all
close all

n=2;
FWHM=4.5*10^-14;
u=20*pi*10^12;

alpha=4*log(2)/(FWHM)^2;
P=zeros(n,n);

P1=zeros(n,n);
P2=zeros(n,n);
P2(1,1)=1;
P(1,1)=1;
P1(1,1)=1;
N=10000;
t=linspace(0,0.5e-12,N);
dt=t(2)-t(1);
e=zeros(N);
E=linspace(0,40,1000);
p=zeros(1000);
q=zeros(1000);
o=zeros(1000);
R=zeros(1000);
 Omega_Rabbi=zeros(N);
for j=1:1000
     
    for i=2:N
      
     e(i-1)=E(j)*exp(-alpha*(t(i)-FWHM)^2);
    
     Omega_Rabbi(i-1)= u*e(i-1)*0.5;
     H=[30*1e4,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
     H1=[30*1e12,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
     H2=[30*1e14,Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
     
     K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
     K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
     K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
     K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1 )+K3))*dt;
    
     l1=1i*(P1(:,:,i-1)*H1-H1*P1(:,:,i-1))*dt;
     l2=1i*((P1(:,:,i-1)+0.5.*l1)*H1-H1*(P1(:,:,i-1)+0.5.*l1))*dt;
     l3=1i*((P1(:,:,i-1)+0.5.*l2)*H1-H1*(P1(:,:,i-1)+0.5.*l2))*dt;
     l4=1i*((P1(:,:,i-1)+l3)*H1-H1*(P1(:,:,i-1 )+l3))*dt;
     
     m1=1i*(P2(:,:,i-1)*H2-H2*P2(:,:,i-1))*dt;
     m2=1i*((P2(:,:,i-1)+0.5.*m1)*H2-H2*(P2(:,:,i-1)+0.5.*m1))*dt;
     m3=1i*((P2(:,:,i-1)+0.5.*m2)*H2-H2*(P2(:,:,i-1)+0.5.*m2))*dt;
     m4=1i*((P2(:,:,i-1)+m3)*H2-H2*(P2(:,:,i-1 )+m3))*dt;
    
     P(:,:,i)  = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
     P1(:,:,i) = P1(:,:,i-1)+((l1+2.*(l2+l3)+l4)/6);
     P2(:,:,i) = P2(:,:,i-1)+((m1+2.*(m2+m3)+m4)/6);
     if i==1000
         R(j)=Omega_Rabbi(i-1);
          break
      end
            
   
    end
   p(j)=P(2,2,i);
   q(j)=P1(2,2,i);
   o(j)=P2(2,2,i);
 % disp(p(j));
    P=zeros(n,n);
    P(1,1,1)=1;    %values for reference
end
plot(E,p,E,q,E,o,'linewidth',2);
%
xlabel('Rabi Frequency','fontSize',14);
ylabel('Ground State Population','fontsize',14);
legend({'30*1e4','30*1e12','30*1e14'},'Location','southwest')
axis([0 E(1000) 0 1.1])
