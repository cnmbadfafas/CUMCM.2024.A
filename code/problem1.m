clc;clear;close all
v_longtou=1;%龙头速度
cita_chu=zeros(224,301);%cita角群
cita_chu(1,1)=32*pi;%初始极角
v=zeros(224,301);%速度群
b=0.55./(2.*pi);%螺距
a=0.55*16-32*pi*b;%
r=zeros(224,301);%极径群
r(1,1)=0.55*16;%极径
long_longtou=3.41;%龙头长度(m)
longtoubashou_d=3.41-0.275*2;%龙头前把手到龙身前把手的距离
long=2.2;%龙身和龙尾距离(m)
longshenbashou_d=2.2-0.275*2;%龙身前把手到龙身前把手的距离
longshen_quantity=221;%龙身数量
longtou_quantity=1;%龙头和龙尾数量分别为1
width=0.3;%板凳宽
dia=0.055;%孔的直径
d=0.275;%孔中心到板头的距离

%%初始位置
%龙头把手
x=zeros(224,301);%各把手初始位置的横坐标
y=zeros(224,301);%各把手初始位置的纵坐标
x(1,1)=0.55*16;%龙头初始位置据中心距离
y(1,1)=0;
syms k
% cita_chu(2,1)=acos(r(1,1).^2+(r(1,1)-b*cita_chu(1,1)+b*(cita_chu(2,1)-cita_chu(1,1))).^2-longtoubashou_d.^2)./(2.*r(1,1).*(r(1,1)-b*cita_chu(1,1)+b*(cita_chu(2,1)-cita_chu(1,1))))+cita_chu(1,1);
% eqn=k-acos((r(1,1).^2+(r(1,1)-cita_chu(1,1)+b*(k-cita_chu(1,1))).^2-longtoubashou_d.^2)./(2.*r(1,1).*(r(1,1)-cita_chu(1,1)+b*(k-cita_chu(1,1)))))-cita_chu(1,1)==0;
eqn=cos(k-cita_chu(1,1))==(r(1,1).^2+(b*k).^2-longtoubashou_d.^2)./(2.*r(1,1).*b*k);%余弦定理解角度
% options = optimoptions('fsolve','Display','iter');
[k_result] = vpasolve(eqn,k,[32*pi 34*pi]);
cita_chu(2,1)=double(k_result);
r(2,1)=b*cita_chu(2,1);
x(2,1)=r(2,1)*cos(cita_chu(2,1));
y(2,1)=r(2,1)*sin(cita_chu(2,1));
sum_low=32*pi+cita_chu(2,1)-cita_chu(1,1);
sum_high=34*pi+cita_chu(2,1)-cita_chu(1,1);
%龙身和龙尾把手
for i=2:223
eqn=cos(k-cita_chu(i,1))==(r(i,1).^2+(b*k).^2-longshenbashou_d.^2)./(2.*r(i,1).*b*k);
[k_result] = vpasolve(eqn,k,[sum_low sum_high]);
cita_chu(i+1,1)=double(k_result);
r(i+1,1)=b*cita_chu(i+1,1);
x(i+1,1)=r(i+1,1)*cos(cita_chu(i+1,1));
y(i+1,1)=r(i+1,1)*sin(cita_chu(i+1,1));
sum_low=sum_low+cita_chu(i+1,1)-cita_chu(i,1);
sum_high=sum_high+cita_chu(i+1,1)-cita_chu(i,1);
end

%%300s位置迭代计算
%龙头龙身把手位移迭代，速度
q=1;
mytheta=@(t,theta)-1./(b*sqrt(1+theta.^2));
theta0=2*pi*16; % 初始位置时候的角度
dt=0.1; % 时间步长
tspan=0:dt:300; % 求解时间点
[tt,theta]=ode45(mytheta,tspan,theta0); % 龙格库塔法求解
X1=b*theta.*cos(theta);
Y1=b*theta.*sin(theta);
for i=1:10:length(theta)
    cita_chu(1,q)=theta(i);
    r(1,q)=b*cita_chu(1,q);
    x(1,q)=r(1,q)*cos(cita_chu(1,q));
    y(1,q)=r(1,q)*sin(cita_chu(1,q));
    q=q+1;
end

for j=1:300
%龙头
% delta_cita=v_longtou./r(1,j);%每秒移动cita
% cita_chu(1,j+1)=cita_chu(1,j)-delta_cita;
% 
% r(1,j+1)=0.55./(2.*pi)*cita_chu(1,j+1);%第二个极径长度
% x(1,j+1)=r(1,j+1)*cos(cita_chu(1,j+1));
% y(1,j+1)=r(1,j+1)*sin(cita_chu(1,j+1));
%龙头速度
v(1,j+1)=v_longtou;%龙头速度

%龙身
eqn=cos(k-cita_chu(1,j+1))==(r(1,j+1).^2+(b*k).^2-longtoubashou_d.^2)./(2.*r(1,j+1).*b*k);%余弦定理解角度

[k_result] = vpasolve(eqn,k,[cita_chu(1,j+1) cita_chu(1,j+1)+2*pi]);
cita_chu(2,j+1)=double(k_result);
r(2,j+1)=b*cita_chu(2,j+1);
x(2,j+1)=r(2,j+1)*cos(cita_chu(2,j+1));
y(2,j+1)=r(2,j+1)*sin(cita_chu(2,j+1));
sum_low=cita_chu(1,j+1)+cita_chu(2,j+1)-cita_chu(1,j+1);
sum_high=cita_chu(1,j+1)+2*pi+cita_chu(2,j+1)-cita_chu(1,j+1);
 %径向速度
%   if j<300
% v_r=abs(r(2,j)-r(2,j+2))./2;
%  %切向速度
% v_q=r(2,j+1)*abs(cita_chu(2,j)-cita_chu(2,j+2))./2;
%  elseif j>=300
% v_r=abs(r(2,j+1)-r(2,j));
% v_q=r(2,j+1)*abs(cita_chu(2,j+1)-cita_chu(2,j));
%   end
% v(2,j+1)=sqrt(v_r.^2+v_q.^2);%总速度
for i=2:223
eqn=cos(k-cita_chu(i,j+1))==(r(i,j+1).^2+(b*k).^2-longshenbashou_d.^2)./(2.*r(i,j+1).*b*k);
[k_result] = vpasolve(eqn,k,[sum_low sum_high]);
cita_chu(i+1,j+1)=double(k_result);
r(i+1,j+1)=b*cita_chu(i+1,j+1);
x(i+1,j+1)=r(i+1,j+1)*cos(cita_chu(i+1,j+1));
y(i+1,j+1)=r(i+1,j+1)*sin(cita_chu(i+1,j+1));
sum_low=sum_low+cita_chu(i+1,j+1)-cita_chu(i,j+1);
sum_high=sum_high+cita_chu(i+1,j+1)-cita_chu(i,j+1);
%龙身速度
%  %径向速度
%  if j<300
% v_r=abs(r(i+1,j+2)-r(i+1,j))./2;
% v_q=r(i+1,j+1)*abs(cita_chu(i+1,j+2)-cita_chu(i+1,j))./2;
%  elseif j>=300
% v_r=abs(r(i+1,j+1)-r(i+1,j));
% v_q=r(i+1,j+1)*abs(cita_chu(i+1,j+1)-cita_chu(i+1,j));
%  end
%  %切向速度
% 
% v(i+1,j+1)=sqrt(v_r.^2+v_q.^2);%总速度
end
end

%求速度
for j=1:300
for i=1:223
 if j<300
v_r=abs(r(i+1,j+2)-r(i+1,j))./2;%径向速度
v_q=r(i+1,j+1)*abs(cita_chu(i+1,j+2)-cita_chu(i+1,j))./2;%切向速度
 elseif j>=300
v_r=abs(r(i+1,j+1)-r(i+1,j));
v_q=r(i+1,j+1)*abs(cita_chu(i+1,j+1)-cita_chu(i+1,j));
 end
 v(i+1,j+1)=sqrt(v_r.^2+v_q.^2);%总速度
end
end


%作图
figure(1)
theta=16*2*pi:-0.01:0*pi;
r1=b*theta;
x1=r1.*cos(theta);
y1=r1.*sin(theta);
figure(1)
set(gcf,'Position',[200 200 600 600]);
plot(x1,y1,'--')
axis equal
grid on
hold on
title('板凳位置变化');
xlabel('x 坐标');
ylabel('y 坐标');
plot(x(:,1),y(:,1),'-o');
hold on
plot(x(:,301),y(:,301),'-o');
legend('螺线','各板凳初始位置', '各板凳最终位置');
%%信息统计
table=zeros(448,301);
for i=1:224

    table(2*i-1,:)=x(i,:);
    table(2*i,:)=y(i,:);
end
writematrix(table,'table.xlsx');
writematrix(v,'v.xlsx');






