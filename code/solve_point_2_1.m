


function [x,y,cita,flag0]=solve_point_2_1(b_d,x1,y1,cita1,d,RC1,x_c,y_c,cita_max)

b=b_d./2./pi;
delta_cita=2.*asin(d./2./RC1); 
if delta_cita<=cita_max-cita1 
flag0=2; 
cita=cita1+delta_cita;
x=x_c+RC1.*cos(cita);
y=y_c+RC1.*sin(cita); 
else
cita=solve_cita1(b_d,x1,y1,4.5./b,d); 
flag0=1; 
x=b.*cita.*cos(cita);
y=b.*cita.*sin(cita); 
end
end