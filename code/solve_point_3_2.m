

function [x,y,cita,flag0]=solve_point_3_2(x1,y1,cita1,d,RC1,x_c1,y_c1,RC2,x_c2,y_c2,cita_min)

delta_cita=2.*asin(d./2./RC2);
if delta_cita<=cita1-cita_min 
flag0=3; 
cita=cita1-delta_cita;
x=x_c2+RC2.*cos(cita);
y=y_c2+RC2.*sin(cita); 
else
di=sqrt((x1-x_c1)^2+(y1-y_c1)^2); 
delta_cita=acos((di^2+RC1^2-d^2)./2./di./RC1); 
cita_C1_di=atan((y1-y_c1)./(x1-x_c1)); 
cita=cita_C1_di+delta_cita; 
flag0=2; 
x=x_c1+RC1.*cos(cita);
y=y_c1+RC1.*sin(cita);
end
end