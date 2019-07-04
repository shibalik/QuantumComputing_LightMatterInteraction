clc
close all
FWHM=3*10^-1;
u=2*pi;

alpha=log(2)/(FWHM)^2;

N=100;
P=[1,0;0,0];
D=linspace(-80,80,N);
E=linspace(0,40,N);
t=linspace(0,0.5,N);
dt=t(2)-t(1);
p=zeros(N,N);
Ham=zeros(2,2,N,N);
for i=1:N
    P(:,:,i)=[1,0;0,0];
    Del(:,:,i)=[D(i),0;0,0];
    H(:,:,i)=[0,u/2;u/2,0];
end
for i=1:N-1
    H(:,:,:,i+1)=H(:,:,:,1).*exp(-alpha*t(i)^2);
    Del(:,:,:,i+1)=Del(:,:,:,i);
end  
for j=1:N
    Ham=Del + E(j).*H;
    
   for i=1:N
    
      
      
      
      K1=1i.*(P(:,:,:,i).*Ham(:,:,:,i)-Ham(:,:,:,i).*P(:,:,:,i)).*dt;
      K2=1i.*((P(:,:,:,i)+0.5.*K1).*Ham(:,:,:,i)-Ham(:,:,:,i).*(P(:,:,:,i)+0.5.*K1)).*dt;
      K3=1i.*((P(:,:,i)+0.5.*K2).*Ham(:,:,:,i)-Ham(:,:,:,i).*(P(:,:,:,i)+0.5.*K2)).*dt;
      K4=1i.*((P(:,:,:,i)+K3).*Ham(:,:,:,i)-Ham(:,:,:,i).*(P(:,:,:,i)+K3)).*dt;
   
      P(:,:,:,i+1) = P(:,:,:,i)+((K1+2.*(K2+K3)+K4)./6);
      H(:,:,i);
   end
   p(:,j)=P(2,2,:,N);
   p(1,j);
end
[X,Y]=meshgrid(E,D);
figure(1)
mesh(X,Y,abs(p));

xlabel('Rabi Frequency','fontSize',14);
ylabel('Detuning', 'fontSize',14);
zlabel('Population','fontsize',14)

