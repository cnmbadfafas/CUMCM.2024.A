
function [x,y,cita,flag0]=solve_point_4_3(b_d,x1,y1,cita1,d,RC2,x_c,y_c,cita_max)

b=b_d./2./pi;
cita=solve_cita2(b_d,x1,y1,cita1,d); 
if cita>=4.5./b-pi 
flag0=4;
x=b.*(cita+pi).*cos(cita);
y=b.*(cita+pi).*sin(cita); 
else
fun=@(t)(x_c+RC2.*cos(cita_max-t)-x1).^2+(y_c+RC2.*sin(cita_max-t)-y1).^2-d^2;
q=-0.1;
options = optimoptions('fsolve','Display','off');
delta_cita=fsolve(fun,cita_max+q,options); 
cita=cita_max-delta_cita; 
flag0=3; 
x=x_c+RC2.*cos(cita);
y=y_c+RC2.*sin(cita); 
    
end
end