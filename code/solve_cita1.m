

function cita=solve_cita1(b_d,x1,y1,cita1,d) 
b=b_d./2./pi;
fun=@(cita)(b.*cita.*cos(cita)-x1).^2+(b.*cita.*sin(cita)-y1).^2-d^2;
q=0.01;
options = optimoptions('fsolve','Display','off'); 
cita=fsolve(fun,cita1+q,options); 
while cita<=cita1 || abs(b.*cita-b.*cita1)>b_d./2 
q=q+0.1;
cita=fsolve(fun,cita+q,options); 
end 
end