clc;close all;clear

b_d=1.7; % 螺距
b=b_d./2./pi; % 螺线方程的系数
longshen=3.41;
longtoubashou=longshen-0.275.*2; % 龙头把手两个孔之间的距离
longshen2=2.2;
longshenbashou=longshen2-0.275.*2; % 其他凳子把手两个孔之间的距离
v_longtou=1; % 头节点速度
N=223; 
%% 盘出螺线与盘入螺线中心对称：

cita=5.*2.*pi:-0.01:0;
r=b.*cita;
x=r.*cos(cita);
y=r.*sin(cita);
figure(1)
set(gcf,'Position',[200 200 600 600]);
c_1=rand(1,3);
plot(x,y,'--','Color',c_1,'LineWidth',1.2)
axis equal
grid on
xlabel('x')
ylabel('y')
set(gca,'FontSize',17) 
hold on

cita=cita-pi;
R2=b.*(cita+pi); % 旋转180°
x2=R2.*cos(cita);
y2=R2.*sin(cita);
c_2=rand(1,3);
plot(x2,y2,'m--','Color',c_2,'LineWidth',1.2)

R=4.5; % 掉头区域半径
x_diao=R.*cos(cita);
y_diao=R.*sin(cita);
c_3=sort(rand(1,3));
plot(x_diao,y_diao,'Color',c_3,'LineWidth',1.9)


%% 写出入螺线和出螺线与圆形边界的交点角度./交点处的切线斜率等几何信息
cita_ru=R./b;
cita_chu=R./b-pi; 
slo=(b.*sin(cita_ru)+R.*cos(cita_ru))./(b.*cos(cita_ru)-R.*sin(cita_ru));
cita_max_1=atan(-1./slo)+pi; 
cita_dengyao=atan(tan(cita_ru))+pi-cita_max_1; 
RC1_C2=R./cos(cita_dengyao); %得到r1+r2的值
RC2=RC1_C2./3;
RC2=RC2-0.0027;
RC1=RC2.*2; 
phi=2.*cita_dengyao; 
SC1=RC1.*(pi-phi); SC2=RC2.*(pi-phi); 
cita_min_1=cita_max_1-SC1./RC1;  
cita_min_2=cita_min_1-pi;
cita_max_2=cita_min_2+SC2./RC2; 
x_C1=R.*cos(cita_ru)+RC1.*cos(cita_max_1-pi);
y_C1=R.*sin(cita_ru)+RC1.*sin(cita_max_1-pi); 

x_C2=R.*cos(cita_chu)-RC2.*cos(cita_max_2);
y_C2=R.*sin(cita_chu)-RC2.*sin(cita_max_2); 
figure(1)
hold on
plot(x_C1+RC1.*cos(linspace(cita_min_1,cita_max_1,50)),y_C1+RC1.*sin(linspace(cita_min_1,cita_max_1,50)),'b','LineWidth',1.5)
plot(x_C1,y_C1,'ro')
plot(x_C2+RC2.*cos(linspace(cita_min_2,cita_max_2,50)),y_C2+RC2.*sin(linspace(cita_min_2,cita_max_2,50)),'r','LineWidth',1.5)
plot(x_C2,y_C2,'bo')
legend('盘入螺线','盘出螺线','调头边界')
axis equal

%% 盘入曲线上头节点的位置求解
mycita=@(t,cita)1./(b.*sqrt(1+cita.^2));
cita0=cita_ru; 
dt=0.1; 
dt=1./randi([5 20]);
tspan=0:dt:100; % 求解时间点
[tt,cita]=ode45(mycita,tspan,cita0); 
X1=b.*cita.*cos(cita);
Y1=b.*cita.*sin(cita);

tt_ru=tt(end:-1:1);
X_save=zeros(224,200./dt+1);
Y_save=zeros(224,200./dt+1); 
TH=zeros(224,200./dt+1);
X_save(1,1:length(X1))=X1(end:-1:1); 
Y_save(1,1:length(Y1))=Y1(end:-1:1); 
TH(1,1:length(cita))=cita(end:-1:1); 

tt_c1=dt:dt:SC1; 
cita_C1=-tt_c1./RC1+cita_max_1;
TH(1,length(cita)+(1:length(tt_c1)))=cita_C1;
X_save(1,length(X1)+(1:length(tt_c1)))=RC1.*cos(cita_C1)+x_C1; 
Y_save(1,length(Y1)+(1:length(tt_c1)))=RC1.*sin(cita_C1)+y_C1; 

tt_c2=tt_c1(end)+dt:dt:SC1+SC2; 
cita_C2=(tt_c2-SC1)./RC2+cita_min_1-pi; 
TH(1,length(cita)+length(cita_C1)+(1:length(tt_c2)))=cita_C2;
X_save(1,length(X1)+length(cita_C1)+(1:length(tt_c2)))=RC2.*cos(cita_C2)+x_C2; 
Y_save(1,length(Y1)+length(cita_C1)+(1:length(tt_c2)))=RC2.*sin(cita_C2)+y_C2; 

mycita=@(t,cita)1./(b.*sqrt(1+(cita+pi).^2));
cita0=cita_chu; 
tspan=tt_c2(end)+dt:dt:100; 
[tt,cita2]=ode45(mycita,tspan,cita0); 
X2=b.*(cita2+pi).*cos(cita2);
Y2=b.*(cita2+pi).*sin(cita2);
TH(1,length(cita)+length(cita_C1)+length(tt_c2)+1:end)=cita2;
X_save(1,length(cita)+length(cita_C1)+length(tt_c2)+1:end)=X2;
Y_save(1,length(cita)+length(cita_C1)+length(tt_c2)+1:end)=Y2; 
figure(3)
set(gcf,'Position',[300 300 600 600])
clf

for i=1:3:length(TH(1,:))
title({['t=',num2str((i-1).*dt)],'龙头把手中心的轨迹(-100s到100s)'})
plot(X_save(1,i),Y_save(1,i),'Marker','o','MarkerSize',3,'MarkerFaceColor','g')
hold on
axis equal
axis([-15 15 -15 15])
grid on
drawnow
end


t_total=-100:dt:100;
for jj=(1:length(t_total))
j=round((t_total(1)+100)./dt)+jj;
if t_total(jj)<=0 
for i=2:N+1 
d=longtoubashou.*(i<=2)+longshenbashou.*(i>2);
citaij=solve_cita1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d);
TH(i,j)=citaij;
X_save(i,j)=b.*citaij.*cos(citaij);
Y_save(i,j)=b.*citaij.*sin(citaij);

end
elseif t_total(jj)>0 && t_total(jj)<=SC1
flag0=2;
for i=2:N+1
d=longtoubashou.*(i<=2)+longshenbashou.*(i>2);
if flag0==2 
[xi,yi,citai,flag0]=solve_point_2_1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d,RC1,x_C1,y_C1,cita_max_1);
TH(i,j)=citai;
X_save(i,j)=xi;
Y_save(i,j)=yi;
else  
citaij=solve_cita1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d); 
TH(i,j)=citaij;
X_save(i,j)=b.*citaij.*cos(citaij);
Y_save(i,j)=b.*citaij.*sin(citaij);
end
end
elseif t_total(jj)>SC1 && t_total(jj)<=SC1+SC2
flag0=3;
for i=2:N+1
d=longtoubashou.*(i<=2)+longshenbashou.*(i>2); 
if flag0==3
[xi,yi,citai,flag0]=solve_point_3_2(X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d,RC1,x_C1,y_C1,RC2,x_C2,y_C2,cita_min_2);
TH(i,j)=citai;
X_save(i,j)=xi;
Y_save(i,j)=yi;
elseif flag0==2 
[xi,yi,citai,flag0]=solve_point_2_1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d,RC1,x_C1,y_C1,cita_max_1);
TH(i,j)=citai;
X_save(i,j)=xi;
Y_save(i,j)=yi;
else  
citaij=solve_cita1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d); 
TH(i,j)=citaij;
X_save(i,j)=b.*citaij.*cos(citaij);
Y_save(i,j)=b.*citaij.*sin(citaij);
end
end
else 
flag0=4;
for i=2:N+1
d=longtoubashou.*(i<=2)+longshenbashou.*(i>2); 
if flag0==4
[xi,yi,citai,flag0]=solve_point_4_3(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d,RC2,x_C2,y_C2,cita_max_2);
TH(i,j)=citai;
X_save(i,j)=xi;
Y_save(i,j)=yi;
elseif flag0==3 
[xi,yi,citai,flag0]=solve_point_3_2(X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d,RC1,x_C1,y_C1,RC2,x_C2,y_C2,cita_min_2);
TH(i,j)=citai;
X_save(i,j)=xi;
Y_save(i,j)=yi;
elseif flag0==2 
[xi,yi,citai,flag0]=solve_point_2_1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d,RC1,x_C1,y_C1,cita_max_1);
TH(i,j)=citai;
X_save(i,j)=xi;
Y_save(i,j)=yi;
else  
citaij=solve_cita1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d); 
TH(i,j)=citaij;
X_save(i,j)=b.*citaij.*cos(citaij);
Y_save(i,j)=b.*citaij.*sin(citaij);
end
end
end
end

%% 可视化调头过程！
figure(100)
clf;
set(gcf,'Position',[200 200 600 600])
for j=(1:2:length(t_total))+round((t_total(1)+100)./dt)
plot(X_save(:,j),Y_save(:,j),'k-o')
title({['t=',num2str(dt.*(j-1)-100)],'轨迹图'})
axis equal
grid on
xlabel('x')
ylabel('y')
axis([-17 17 -17 17])
drawnow

end

Vx=zeros(size(X_save));
Vy=Vx;
Vx(:,1)=(X_save(:,2)-X_save(:,1))./dt; 
Vx(:,end)=(X_save(:,end)-X_save(:,end-1))./dt; 
Vx(:,2:end-1)=(X_save(:,3:end)-X_save(:,1:end-2))./2./dt; 

Vy(:,1)=(Y_save(:,2)-Y_save(:,1))./dt; 
Vy(:,end)=(Y_save(:,end)-Y_save(:,end-1))./dt; 
Vy(:,2:end-1)=(Y_save(:,3:end)-Y_save(:,1:end-2))./2./dt; 

V=sqrt(Vx.^2+Vy.^2);


figure
plot(t_total,V(1,:),'b.','LineWidth',1.4)
ylim([0 1.1])
xlabel('时间')
ylabel('头把手的速度')
title('速度')

format long
nn=1./dt; 
index=1:nn:length(t_total); 
Dataxy=zeros(2.*(N+1),length(index)); 
Dataxy(1:2:end,:)=round(X_save(:,index),6); 
Dataxy(2:2:end,:)=round(Y_save(:,index),6);
Datav=round(V(:,index),6); 
% 写入文件
filename = 'result44_test.xlsx';
sheet = 1;
xlRange = 'B2';
xlswrite(filename,Dataxy,sheet,xlRange)  
sheet = 2;
xlRange = 'B2';
xlswrite(filename,Datav,sheet,xlRange)  


nn2=50./dt; 
index2=1:nn2:length(t_total); 
index_row=[1 2:50:224 224];
Dataxy2=zeros(2.*length(index_row),length(index2));
Dataxy2(1:2:end,:)=round(X_save(index_row,index2),6);
Dataxy2(2:2:end,:)=round(Y_save(index_row,index2),6) 
Datav2=round(V(index_row,index2),6) 

filename = 'result44_test_2.xlsx';
sheet = 1;
xlRange = 'B2';
xlswrite(filename,Dataxy2,sheet,xlRange)  % 写入位置坐标
sheet = 2;
xlRange = 'B2';
xlswrite(filename,Datav2,sheet,xlRange)  % 写入速度值

