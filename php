clear all; close all;

dt=1e-4;
Td=0.05;
i=0;
I=0;
W=0;
V=1;
TL=0;
RR=100;

x1=0.2;%球的位置
x2=0;%球的速率
x3=60;%桿子的角度
x4=0;%桿子角速度
u=0;

mB=0.029;%球的重量
mb=0.334;%桿子重量
rB=0.0095;%球的半徑
l=0.4;%桿子長度
d=0.04;%銜接距離
JB=1.05e-6;%球的慣量距
Jb=0.0178;%桿子的慣量距
Kb=Kt=0.1491;%常數
Ra=18.91;%電阻
n=4.2;%齒輪比例
g=9.8;%加速度

G=tf([mB*g*d],[l*(JB/rB^2+mB) 0 0])%步階響應圖
step(G,10)
xlabel( 'Time(sec)');
ylabel('Amplitude');  

pMF1=[0 0.2 0.4];
pMF2=[-1 0 1];
pMFO=[-1 0 1];
K0=20;
err_temp=0;

for t=0:dt:Td
  i=i+1;
  K1=(1+mB^-1*JB.*rB^-2)^-1;
  K2=Kb*Kt*l^2.*(Ra*d^2)^-1;
  K3=n*Kb*l.*(Ra*d)^-1;
  K4(i)=(JB+Jb+mB*x1(i)*x1(i))^-1;
 
  
  err=RR-W(i)*30/pi;
  derr=(err-err_temp)/dt;
  err_temp=err;
  
  in1=err;
  in2=derr;
  
  in1=in1*K1;
if (in1<=pMF1(1))
  uMF11=1;
  uMF12=0;
  uMF13=0;
elseif ((pMF1(1)<=in1)&&(in1<pMF1(2)))
  uMF11=(pMF1(2)-in1)/(pMF1(2)-pMF1(1));
  uMF12=1-uMF11;
  uMF13=0;
elseif ((pMF1(2)<=in1)&&(in1<pMF1(3)))
  uMF11=0;
  uMF12=(pMF1(3)-in1)/(pMF1(3)-pMF1(2));
  uMF13=1-uMF12;
else
  uMF11=0;
  uMF12=0;
  uMF13=1;
endif

in2=in2*K2;
if (in2<=pMF2(1))
  uMF21=1;
  uMF22=0;
  uMF23=0;
elseif ((pMF2(1)<=in2)&&(in2<pMF2(2)))
  uMF21=(pMF2(2)-in2)/(pMF2(2)-pMF2(1));
  uMF22=1-uMF21;
  uMF23=0;
elseif ((pMF2(2)<=in2)&&(in2<pMF2(3)))
  uMF21=0;
  uMF22=(pMF2(3)-in2)/(pMF2(3)-pMF2(2));
  uMF23=1-uMF22;
else
  uMF21=0;
  uMF22=0;
  uMF23=1;
endif

  uMF11
  uMF12
  uMF13
  uMF21
  uMF22
  uMF23
  
  uMF01=max([min(uMF11,uMF21),min(uMF12,uMF21),min(uMF11,uMF22)]);
  uMF02=max([min(uMF13,uMF21),min(uMF12,uMF22),min(uMF11,uMF23)]);
  uMF03=max([min(uMF13,uMF22),min(uMF12,uMF23),min(uMF13,uMF23)]);
  
  Out=(uMF01*pMFO(1)+uMF02*pMFO(2)+uMF03*pMFO(3))/(uMF01+uMF02+uMF03)
  
  V=K0*Out;
  x1(i+1)=x2(i)*dt+x1(i);
  K4(i+1)=(JB+Jb+(mB*x1(i)*x1(i))).^-1*dt+K4(i);
  x2(i+1)=(x1(i)*x4(i)*x4(i)-g*sin(x3(i)))*K1*dt+x2(i);
  x3(i+1)=x4(i).*dt+x3(i);
  x4(i+1)=(K4(i).*cos(x3(i)).*(K2.*x4(i).*cos((l.*x3(i))/d)+K3.*cos(l*x3(i)/d).*V(i)
          -0.5*l*mb*g*x1(i))-2*mB.*x1(i).*x2(i).*x4(i).*K4(i))*dt+x4(i);
endfor

plot(T,x1(1:length(T)));
figure
%plot(T,W(1:length(T))*30/pi);
xlabel('time(sec)')
ylabel('Velocity(rad/sec)')
