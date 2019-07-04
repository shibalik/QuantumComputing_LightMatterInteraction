clc
close all

n=2;
FWHM=3*10^-2;
u=2*pi;

alpha=4*log(2)/(FWHM)^2;
P=zeros(n,n);
P(1,1)=1;
N=1000;
t=linspace(0,0.15,N);
dt=t(2)-t(1);
e=zeros(N);
E=linspace(0,120,100);
p=zeros(150,100);
R=zeros(N);
Omega_Rabbi=zeros(N);
B=linspace(-400,400,150);
[X,Y]=meshgrid(E,B);
 X(1,:)=X(2,:);
for k=1:150
  for j=1:100
    for i=1:N
    
      e(i)=E(j)*exp(-alpha*(t(i)-FWHM)^2); 
         Omega_Rabbi(i)= u*e(i)/(2);
     

      H=[B(k)-2*pi*0.22*(t(i)-FWHM)/3,Omega_Rabbi(i);Omega_Rabbi(i),0];
      %H1=[0,10*3.14*10^12;10*3.14*10^12,0];
      K1=1i*(P(:,:,i)*H-H*P(:,:,i))*dt;
      K2=1i*((P(:,:,i)+0.5.*K1)*H-H*(P(:,:,i)+0.5.*K1))*dt;
      K3=1i*((P(:,:,i)+0.5.*K2)*H-H*(P(:,:,i)+0.5.*K2))*dt;
      K4=1i*((P(:,:,i)+K3)*H-H*(P(:,:,i )+K3))*dt;
   
      P(:,:,i+1) = P(:,:,i)+((K1+2.*(K2+K3)+K4)/6);
       %  if t(i)<0.5*FWHM&&t(i)>2*FWHM
         %  Omega_Rabbi(i)=0;
   %   end
         if t(i)>4*FWHM
             break
         end

    end
    
     X(k,j)= Omega_Rabbi(i);
   p(k,j)=P(2,2,i);
   %disp(p(j));
   P=zeros(n,n);
   P(1,1,1)=1;    %values for reference
  end
  
 
   %disp(p(j));
   P=zeros(n,n); 
   P(1,1,1)=1; 
end 

figure(1)
surf(X,Y,abs(p));

xlabel('Rabi Frequency','fontSize',14);
ylabel('Detuning', 'fontSize',14);
zlabel('Population','fontsize',14)
figure(2)
contour(X,Y,abs(p));