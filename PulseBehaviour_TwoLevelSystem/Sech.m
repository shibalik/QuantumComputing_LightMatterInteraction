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
E=linspace(0,100,100);
p=zeros(100,100);
R=zeros(100);
 Omega_Rabbi=zeros(N);
 B=linspace(-180,180,100);
 [X,Y]=meshgrid(E,B);
for k=1:100
    P=zeros(n,n);
    P(1,1,1)=1;
  for j=1:100
     P(:,:,1);
     for i=2:N
      X(j,k);
       e(i-1)=X(j,k)*sech(2*log(2+sqrt(3))*(t(i-1)-FWHM)/FWHM);
    
      Omega_Rabbi(i-1)= u*e(i-1)/2;
        if t(i)<0.5*FWHM&&t(i)>2*FWHM
           Omega_Rabbi(i)=0;
      end

      H=[Y(j,k),Omega_Rabbi(i-1);Omega_Rabbi(i-1),0];
    %H1=[0,10*3.14*10^12;10*3.14*10^12,0];
      K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
      K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
      K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
      K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1 )+K3))*dt;
     
      P(:,:,i) = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
    
              if t(i)<0.5*FWHM&&t(i)>2*FWHM
           Omega_Rabbi(i)=0;
      end
     end
    
   p(j,k)=P(2,2,N);
   %disp(p(j));
    P=zeros(n,n);
    P(1,1,1)=1;    %values for reference
  end
 
   %disp(p(j));
    P=zeros(n,n);
    P(1,1,1)=1;
 end

figure(1)
mesh(X,Y,real(p));
xlabel('Rabi Frequency','fontSize',14);
ylabel('Detuning', 'fontSize',14);
zlabel('Population','fontsize',14)
figure(2)
contour(X,Y,real(p));