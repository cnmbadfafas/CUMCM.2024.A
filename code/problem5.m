clc;close all;clear

b_d=1.7; % 螺距
b=b_d./2./pi; % 螺线方程的系数 
longtouchangdu=3.41;
longtoubashou=longtouchangdu-0.275.*2; % 龙头把手两个孔之间的距离
longshenchangdu=2.2;
longshenbashou=longshenchangdu-0.275.*2; % 其他凳子把手两个孔之间的距离
V0=1:0.1:2;        
Vmax=zeros(1,length(V0)); 
N=223;
for vv=1:length(V0) 
v0=V0(vv); % 当前的速度
R=4.5;
cita_ru=R./b;
cita_chu=R./b-pi;   
slo=(b.*sin(cita_ru)+R.*cos(cita_ru))./(b.*cos(cita_ru)-R.*sin(cita_ru));
cita_max_1=atan(-1./slo)+pi; 
cita_dengyao=atan(tan(cita_ru))+pi-cita_max_1; 
RC1_C2=R./cos(cita_dengyao); 
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

dt=0.1; 
Ttotal=(SC1+SC2)./v0; 
T=0:dt:Ttotal;
X_save=nan.*zeros(224,length(T)); 
Y_save=nan.*zeros(224,length(T)); 
TH=nan.*zeros(224,length(T));

tt_c1=0:dt:SC1./v0; 
cita_C1=-v0.*tt_c1./RC1+cita_max_1;
TH(1,(1:length(tt_c1)))=cita_C1;
X_save(1,(1:length(tt_c1)))=RC1.*cos(cita_C1)+x_C1;
Y_save(1,(1:length(tt_c1)))=RC1.*sin(cita_C1)+y_C1; 
tt_c2=tt_c1(end)+dt:dt:(SC1+SC2)./v0; 
cita_C2=v0.*(tt_c2-SC1./v0)./RC2+cita_min_1-pi; 
TH(1,length(cita_C1)+(1:length(tt_c2)))=cita_C2;
X_save(1,length(cita_C1)+(1:length(tt_c2)))=RC2.*cos(cita_C2)+x_C2; 
Y_save(1,length(cita_C1)+(1:length(tt_c2)))=RC2.*sin(cita_C2)+y_C2; 


t_total=0:dt:(SC1+SC2)./v0; 
for jj=(1:length(t_total))
j=jj;
if t_total(jj)<0
for i=2:N+1
d=longtoubashou.*(i<=2)+longshenbashou.*(i>2);
citaij=solve_cita1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d);
TH(i,j)=citaij;
X_save(i,j)=b.*citaij.*cos(citaij);
Y_save(i,j)=b.*citaij.*sin(citaij);

end
elseif t_total(jj)>=0 && t_total(jj)<=SC1./v0
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
elseif t_total(jj)>SC1./v0 && t_total(jj)<=SC1./v0+SC2./v0
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
citaij=solve_cita1(b_d,X_save(i-1,j),Y_save(i-1,j),TH(i-1,j),d); % 子函数求解下一个孔的角度值
TH(i,j)=citaij;
X_save(i,j)=b.*citaij.*cos(citaij);
Y_save(i,j)=b.*citaij.*sin(citaij);
end
end
end
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
Vmax(vv)=max(max(V)); 

end

%%可视化
p=polyfit(V0,Vmax,1);
figure
p1=plot(V0,Vmax,'ro','LineWidth',1.1);
xlabel('龙头速度')
ylabel('队伍把手最大速度')
hold on
color=rand(1,3);
p2=plot(linspace(V0(1),V0(end),100),polyval(p,linspace(V0(1),V0(end),100)),'Color',color);
plot([V0(1) V0(end)],[2 2],'m--')
legend('队伍最大速度与龙头速度数据点','线性拟合关系')

v0max=(2-p(2))./p(1);
disp(['龙头最大行进速度为:',num2str(v0max)])
