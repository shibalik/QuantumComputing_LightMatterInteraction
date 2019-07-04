clc
clear all;
close all;

N=1000;
t=linspace(0,0.25e-6,N);
dt=t(2)-t(1);
omega_Rabi=2*pi*10*1e6; 

p11(1:N)=1;
p22(1:N)=0;
p12(1:N)=0;
p21(1:N)=0;
P=zeros(1000);
C=linspace(-50*1e6,50*1e6,1000);
for j=1:N-1
 b=C(j);
  for i=2:N
    
    k1=dp22(t(i-1),omega_Rabi,p12(i-1),p21(i-1),b)*dt;
    l1=dp11(t(i-1),omega_Rabi,p12(i-1),p21(i-1),b)*dt;
    m1=dp12(t(i-1),omega_Rabi,p11(i-1),p22(i-1),p12(i-1),b)*dt;
    n1=dp21(t(i-1),omega_Rabi,p11(i-1),p22(i-1),p21(i-1),b)*dt;
    
    k2=dp22(t(i-1)+0.5*dt,omega_Rabi,p12(i-1)+0.5*m1,p21(i-1)+0.5*n1,b)*dt;
    l2=dp11(t(i-1)+0.5*dt,omega_Rabi,p12(i-1)+0.5*m1,p21(i-1)+0.5*n1,b)*dt;
    m2=dp12(t(i-1)+0.5*dt,omega_Rabi,p11(i-1)+0.5*l1,p22(i-1)+0.5*k1,p12(i-1),b)*dt;
    n2=dp21(t(i-1)+0.5*dt,omega_Rabi,p11(i-1)+0.5*l1,p22(i-1)+0.5*k1,p21(i-1),b)*dt;
    
    k3=dp22(t(i-1)+0.5*dt,omega_Rabi,p12(i-1)+0.5*m2,p21(i-1)+0.5*n2,b)*dt;
    l3=dp11(t(i-1)+0.5*dt,omega_Rabi,p12(i-1)+0.5*m2,p21(i-1)+0.5*n2,b)*dt;
    m3=dp12(t(i-1)+0.5*dt,omega_Rabi,p11(i-1)+0.5*l2,p22(i-1)+0.5*k2,p12(i-1),b)*dt;
    n3=dp21(t(i-1)+0.5*dt,omega_Rabi,p11(i-1)+0.5*l2,p22(i-1)+0.5*k2,p21(i-1),b)*dt;
    
  
    k4=dp22(t(i-1)+dt,omega_Rabi,p12(i-1)+m3,p21(i-1)+n3,b)*dt;
    l4=dp11(t(i-1)+dt,omega_Rabi,p12(i-1)+m3,p21(i-1)+n3,b)*dt;
    m4=dp12(t(i-1)+dt,omega_Rabi,p11(i-1)+l3,p22(i-1)+k3,p12(i-1),b)*dt;
    n4=dp21(t(i-1)+dt,omega_Rabi,p11(i-1)+l3,p22(i-1)+k3,p21(i-1),b)*dt;
    
    p22(i) = p22(i-1)+((k1+2*(k2+k3)+k4)/6);
    p11(i) = p11(i-1)+((l1+2*(l2+l3)+l4)/6);
    p12(i) = p12(i-1)+((m1+2*(m2+m3)+m4)/6);
    p21(i) = p21(i-1)+((n1+2*(n2+n3)+n4)/6);
  end
P(j)=max(p22);
end
figure(1)

%{plot(t,p22,'linewidth',2)
%hold on
%plot(t,p11,'k','linewidth',2)
%plot(t,abs(p12),'r','linewidth',2)
%plot(t,abs(p21),'g','linewidth',2)
%xlabel('time','fontSize',14);
%ylabel('density matrix elements','fontsize',14);
%axis([0 t(N) 0 1.1])

%axis([0 t(N) 0 1.1])
%fh = figure(1);
%set(fh, 'color', 'white');
figure(2)

plot(C,P,'linewidth',2)
xlabel('Detuning','fontsize',14);
ylabel('Population inversion of ground state','fontsize',14);
%axis([C(0) C(N)]
hold on