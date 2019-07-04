clc
clear all;
close all;

prompt= 'input no. of levels in your system ';

N=2;
H=zeros(N,N);
prompt='Rabi Frequency ';
Omega=input(prompt);
H=[0,Omega;Omega,0];
P=zeros(N,N);
P(1,1)=1;
n=1000;
t=linspace(0,0.25,n);
dt=t(2)-t(1);

for i=2:n
    K1=1i*(P(:,:,i-1)*H-H*P(:,:,i-1))*dt;
    K2=1i*((P(:,:,i-1)+0.5.*K1)*H-H*(P(:,:,i-1)+0.5.*K1))*dt;
    K3=1i*((P(:,:,i-1)+0.5.*K2)*H-H*(P(:,:,i-1)+0.5.*K2))*dt;
    K4=1i*((P(:,:,i-1)+K3)*H-H*(P(:,:,i-1)+K3))*dt;
    
    P(:,:,i) = P(:,:,i-1)+((K1+2.*(K2+K3)+K4)/6);
   
end
q=zeros(n);
p=zeros(n);
for i=1:n
    p(i)=P(2,2,i);
    q(i)=P(1,1,i);
    
end
figure(1)

plot(t,p,'linewidth',2);
hold on
plot(t,q,'k','linewidth',2)

xlabel('time','fontSize',14);
ylabel('density matrix elements','fontsize',14);
axis([0 t(n) 0 1.1])

axis([0 t(n) 0 1.1])
fh = figure(1);
set(fh, 'color', 'white');

    
