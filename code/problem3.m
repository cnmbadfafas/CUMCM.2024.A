clc;clear;close all

v_longtou=1;%龙头速度
cita_chu=zeros(224,1);%cita角群
b=0.55./(2.*pi);%初始螺距
a=0.55*16-32*pi*b;%螺线中心为圆心，故a=0;
r=zeros(224,1);%极径群
r(1,1)=4.5;%极径，龙头
long_longtou=3.41;%龙头长度(m)
longtoubashou_d=3.41-0.275*2;%龙头前把手到龙身前把手的距离
long=2.2;%龙身和龙尾距离(m)
longshenbashou_d=2.2-0.275*2;%龙身前把手到龙身前把手的距离
longshen_quantity=221;%龙身数量
longtou_quantity=1;%龙头和龙尾数量分别为1
width=0.3;%板凳宽
dia=0.055;%孔的直径
d=0.275;%孔中心到板头的距离
x=zeros(224,1);%各把手横坐标群
y=zeros(224,1);%各把手纵坐标群
cita_chu(1,1)=r(1,1)./b;%初始极角
x(1,1)=r(1,1).*cos(cita_chu(1,1));%龙头初始位置据中心距离
y(1,1)=r(1,1).*sin(cita_chu(1,1));
time=0;
%%初始位置
db=0.002./(2.*pi);%螺距步长
flag=0;%触碰标志位
syms k
eqn=cos(k-cita_chu(1,1))==(r(1,1).^2+(b*k).^2-longtoubashou_d.^2)./(2.*r(1,1).*b*k);%余弦定理解角度
[k_result] = vpasolve(eqn,k,[cita_chu(1,1) cita_chu(1,1)+2*pi]);
cita_chu(2,1)=double(k_result);
r(2,1)=b*cita_chu(2,1);
x(2,1)=r(2,1)*cos(cita_chu(2,1));
y(2,1)=r(2,1)*sin(cita_chu(2,1));
sum_low=cita_chu(1,1)+cita_chu(2,1)-cita_chu(1,1);
sum_high=cita_chu(1,1)+2*pi+cita_chu(2,1)-cita_chu(1,1);

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

%%盘入迭代
while flag==0
    time=time+1;
    b=b-db;
    cita_chu(1,1)=r(1,1)./b;%变化后龙头极角
    x(1,1)=r(1,1).*cos(cita_chu(1,1));%龙头变化后位置坐标
    y(1,1)=r(1,1).*sin(cita_chu(1,1));
    eqn=cos(k-cita_chu(1,1))==(r(1,1).^2+(b*k).^2-longtoubashou_d.^2)./(2.*r(1,1).*b*k);%余弦定理解角度
    % options = optimoptions('fsolve','Display','iter');
    [k_result] = vpasolve(eqn,k,[cita_chu(1,1) cita_chu(1,1)+2*pi]);
    cita_chu(2,1)=double(k_result);
    r(2,1)=b*cita_chu(2,1);
    x(2,1)=r(2,1)*cos(cita_chu(2,1));
    y(2,1)=r(2,1)*sin(cita_chu(2,1));
    sum_low=cita_chu(1,1)+cita_chu(2,1)-cita_chu(1,1);
    sum_high=cita_chu(1,1)+2*pi+cita_chu(2,1)-cita_chu(1,1);
    
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

    %线性规划过程
     %对龙头:斜率以及四个点的坐标
     k1=(y(1)-y(2))./(x(1)-x(2));
     Alpha=atan(k1);
     if (y(2)-y(1)>=0 && x(2)-x(1)>=0) || (y(2)-y(1)<=0 && x(2)-x(1)>=0)
     X_A1=x(1)+(longtoubashou_d+0.275)*cos(Alpha)+0.15*sin(Alpha);
     X_A4=x(1)+(longtoubashou_d+0.275)*cos(Alpha)-0.15*sin(Alpha);
     X_A2=x(1)-(0.275)*cos(Alpha)+0.15*sin(Alpha);
     X_A3=x(1)-(0.275)*cos(Alpha)-0.15*sin(Alpha);

     Y_A1=y(1)+(longtoubashou_d+0.275)*sin(Alpha)-0.15*cos(Alpha);
     Y_A4=y(1)+(longtoubashou_d+0.275)*sin(Alpha)+0.15*cos(Alpha);
     Y_A2=y(1)-(0.275)*sin(Alpha)-0.15*cos(Alpha);
     Y_A3=y(1)-(0.275)*sin(Alpha)+0.15*cos(Alpha);
     
     elseif (y(2)-y(1)<=0 && x(2)-x(1)<=0) || (y(2)-y(1)>=0 && x(2)-x(1)<=0)

     X_A1=x(1)-(longtoubashou_d+0.275)*cos(Alpha)-0.15*sin(Alpha);
     X_A4=x(1)-(longtoubashou_d+0.275)*cos(Alpha)+0.15*sin(Alpha);
     X_A2=x(1)+(0.275)*cos(Alpha)-0.15*sin(Alpha);
     X_A3=x(1)+(0.275)*cos(Alpha)+0.15*sin(Alpha);

     Y_A1=y(1)-(longtoubashou_d+0.275)*sin(Alpha)+0.15*cos(Alpha);
     Y_A4=y(1)-(longtoubashou_d+0.275)*sin(Alpha)-0.15*cos(Alpha);
     Y_A2=y(1)+(0.275)*sin(Alpha)+0.15*cos(Alpha);
     Y_A3=y(1)+(0.275)*sin(Alpha)-0.15*cos(Alpha);
      
     end
     f=[];
     i=cita_chu(:);
     for m=1:224
     if i(m)-cita_chu(1)>=2*pi && i(m)-cita_chu(2)<=2*pi
        f =[f; find(cita_chu==i(m))];
     end
     end
     X_lisan=[];
     Y_lisan=[];
     if numel(f)==2 %找到两个节点
     %中间
     k3=(y(f(1))-y(f(2)))./(x(f(1))-x(f(2)));
     k7=-1./k1;
     Alpha1=atan(k3);

     if (y(f(2))-y(f(1))>=0 && x(f(2))-x(f(1))>=0) || (y(f(2))-y(f(1))<=0 && x(f(2))-x(f(1))>=0)
     X3_A4=x(f(1))+(longshenbashou_d+0.275)*cos(Alpha1)-0.15*sin(Alpha1);
     X3_A3=x(f(1))-(0.275)*cos(Alpha1)-0.15*sin(Alpha1);

     Y3_A4=y(f(1))+(longshenbashou_d+0.275)*sin(Alpha1)+0.15*cos(Alpha1);
     Y3_A3=y(f(1))-(0.275)*sin(Alpha1)+0.15*cos(Alpha1);

     elseif (y(f(2))-y(f(1))<=0 && x(f(2))-x(f(1))<=0) || (y(f(2))-y(f(1))>=0 && x(f(2))-x(f(1))<=0)

     X3_A4=x(f(1))-(longshenbashou_d+0.275)*cos(Alpha1)+0.15*sin(Alpha1);
     X3_A3=x(f(1))+(0.275)*cos(Alpha1)+0.15*sin(Alpha1);

     Y3_A4=y(f(1))-(longshenbashou_d+0.275)*sin(Alpha1)-0.15*cos(Alpha1);
     Y3_A3=y(f(1))+(0.275)*sin(Alpha1)-0.15*cos(Alpha1);
     end
     X_lisan=[X_lisan;(linspace(X3_A3,X3_A4,20))'];
     Y_lisan=[Y_lisan;(linspace(Y3_A3,Y3_A4,20))'];
     %下面
     k2=(y(f(1)-1)-y(f(1)))./(x(f(1)-1)-x(f(1)));
     Alpha2=atan(k2);

     if (y(f(1))-y(f(1)-1)>=0 && x(f(1))-x(f(1)-1)>=0) || (y(f(1))-y(f(1)-1)<=0 && x(f(1))-x(f(1)-1)>=0)
     X2_A4=x(f(1)-1)+(longshenbashou_d+0.275)*cos(Alpha2)-0.15*sin(Alpha2);
     X2_A3=x(f(1)-1)-(0.275)*cos(Alpha2)-0.15*sin(Alpha2);

     Y2_A4=y(f(1)-1)+(longshenbashou_d+0.275)*sin(Alpha2)+0.15*cos(Alpha2);
     Y2_A3=y(f(1)-1)-(0.275)*sin(Alpha2)+0.15*cos(Alpha2);

     elseif (y(f(1))-y(f(1)-1)<=0 && x(f(1))-x(f(1)-1)<=0) || (y(f(1))-y(f(1)-1)>=0 && x(f(1))-x(f(1)-1)<=0)

     X2_A4=x(f(1)-1)-(longshenbashou_d+0.275)*cos(Alpha2)+0.15*sin(Alpha2);
     X2_A3=x(f(1)-1)+(0.275)*cos(Alpha2)+0.15*sin(Alpha2);

     Y2_A4=y(f(1)-1)-(longshenbashou_d+0.275)*sin(Alpha2)-0.15*cos(Alpha2);
     Y2_A3=y(f(1)-1)+(0.275)*sin(Alpha2)-0.15*cos(Alpha2);
     end
     
     X_lisan=[X_lisan;(linspace(X2_A3,X2_A4,20))'];
     Y_lisan=[Y_lisan;(linspace(Y2_A3,Y2_A4,20))'];
     %上面
     k4=(y(f(2))-y(f(2)+1))./(x(f(2))-x(f(2)+1));
     Alpha3=atan(k4);

     if (y(f(2)+1)-y(f(2))>=0 && x(f(2)+1)-x(f(2))>=0) || (y(f(2)+1)-y(f(2))<=0 && x(f(2)+1)-x(f(2))>=0)
     X4_A4=x(f(2))+(longshenbashou_d+0.275)*cos(Alpha3)-0.15*sin(Alpha3);
     X4_A3=x(f(2))-(0.275)*cos(Alpha3)-0.15*sin(Alpha3);

     Y4_A4=y(f(2))+(longshenbashou_d+0.275)*sin(Alpha3)+0.15*cos(Alpha3);
     Y4_A3=y(f(2))-(0.275)*sin(Alpha3)+0.15*cos(Alpha3);

     elseif (y(f(2)+1)-y(f(2))<=0 && x(f(2)+1)-x(f(2))<=0) || (y(f(2)+1)-y(f(2))>=0 && x(f(2)+1)-x(f(2))<=0)

     X4_A4=x(f(2))-(longshenbashou_d+0.275)*cos(Alpha3)+0.15*sin(Alpha3);
     X4_A3=x(f(2))+(0.275)*cos(Alpha3)+0.15*sin(Alpha3);

     Y4_A4=y(f(2))-(longshenbashou_d+0.275)*sin(Alpha3)-0.15*cos(Alpha3);
     Y4_A3=y(f(2))+(0.275)*sin(Alpha3)-0.15*cos(Alpha3);
     end

     X_lisan=[X_lisan;(linspace(X4_A3,X4_A4,20))'];
     Y_lisan=[Y_lisan;(linspace(Y4_A3,Y4_A4,20))'];

     elseif numel(f)==1 %找到一个节点

     %上面
     k5=(y(f)-y(f+1))./(x(f)-x(f+1));
     Alpha4=atan(k5);

     if (y(f+1)-y(f)>=0 && x(f+1)-x(f)>=0) || (y(f+1)-y(f)<=0 && x(f+1)-x(f)>=0)
     X5_A4=x(f)+(longshenbashou_d+0.275)*cos(Alpha4)-0.15*sin(Alpha4);
     X5_A3=x(f)-(0.275)*cos(Alpha4)-0.15*sin(Alpha4);

     Y5_A4=y(f)+(longshenbashou_d+0.275)*sin(Alpha4)+0.15*cos(Alpha4);
     Y5_A3=y(f)-(0.275)*sin(Alpha4)+0.15*cos(Alpha4);

     elseif  (y(f+1)-y(f)<=0 && x(f+1)-x(f)<=0) || (y(f+1)-y(f)>=0 && x(f+1)-x(f)<=0)

     X5_A4=x(f)-(longshenbashou_d+0.275)*cos(Alpha4)+0.15*sin(Alpha4);
     X5_A3=x(f)+(0.275)*cos(Alpha4)+0.15*sin(Alpha4);

     Y5_A4=y(f)-(longshenbashou_d+0.275)*sin(Alpha4)-0.15*cos(Alpha4);
     Y5_A3=y(f)+(0.275)*sin(Alpha4)-0.15*cos(Alpha4);


     end

     X_lisan=[X_lisan;(linspace(X5_A3,X5_A4,20))'];
     Y_lisan=[Y_lisan;(linspace(Y5_A3,Y5_A4,20))'];

     %下面
     k6=(y(f-1)-y(f))./(x(f-1)-x(f));
     Alpha5=atan(k6);

     if (y(f)-y(f-1)>=0 && x(f)-x(f-1)>=0) || (y(f)-y(f-1)<=0 && x(f)-x(f-1)>=0)
     X6_A4=x(f-1)+(longshenbashou_d+0.275)*cos(Alpha5)-0.15*sin(Alpha5);
     X6_A3=x(f-1)-(0.275)*cos(Alpha5)-0.15*sin(Alpha5);

     Y6_A4=y(f-1)+(longshenbashou_d+0.275)*sin(Alpha5)+0.15*cos(Alpha5);
     Y6_A3=y(f-1)-(0.275)*sin(Alpha5)+0.15*cos(Alpha5);

     elseif (y(f)-y(f-1)<=0 && x(f)-x(f-1)<=0) || (y(f)-y(f-1)>=0 && x(f)-x(f-1)<=0)

     X6_A4=x(f-1)-(longshenbashou_d+0.275)*cos(Alpha5)+0.15*sin(Alpha5);
     X6_A3=x(f-1)+(0.275)*cos(Alpha5)+0.15*sin(Alpha5);

     Y6_A4=y(f-1)-(longshenbashou_d+0.275)*sin(Alpha5)-0.15*cos(Alpha5);
     Y6_A3=y(f-1)+(0.275)*sin(Alpha5)-0.15*cos(Alpha5);
     end

     X_lisan=[X_lisan;(linspace(X6_A3,X6_A4,20))'];
     Y_lisan=[Y_lisan;(linspace(Y6_A3,Y6_A4,20))'];

     end

     for i=1:length(X_lisan)
     if Y_lisan(i)-Y_A1>=k1*(X_lisan(i)-X_A1) && Y_lisan(i)-Y_A2>=k7*(X_lisan(i)-X_A2) && Y_lisan(i)-Y_A3<=k1*(X_lisan(i)-X_A3) && Y_lisan(i)-Y_A4<=k7*(X_lisan(i)-X_A4)
     
         flag=1;
         break;
     end
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
legend('螺线','各板凳最终位置');

%螺距
b_real=b*2*pi+db;
