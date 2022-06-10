function dydt = BackStance(t,y,u,k_vals)

g = u(1);
m = u(2);
I = u(3);
T_stance = u(4);
T_air = u(5);
amp_y = u(6);
amp_x = u(7);
v = u(8);
x_foot = u(9);
hip_des = u(10);

kP_z = k_vals(1);
kD_z = k_vals(2);
kD_x = k_vals(3);
kP_th = k_vals(4);
kD_th = k_vals(5);

%Open Loop or Closed-loop linear system
B_y = amp_y*sin((pi()/T_stance)*(t-(T_air+T_stance)));
B_x = amp_x*sin((pi()/T_stance)*(t-(T_air+T_stance)));

hip_height = y(1) - 0.5*0.7*sin(y(5));
hip_height_dot = y(2) - y(6)*0.7*0.5;

if ((t-(T_air+T_stance))<(0.2*T_stance))
    g_feedback = (1/(0.2*T_stance))*(t-(T_air+T_stance));
elseif ((t-(T_air+T_stance))>0.8*T_stance)
    g_feedback = -(1/(0.2*T_stance))*(t-(T_air+T_stance)) + 5;
else
    g_feedback = 1;
end

F_hip = g_feedback*(-kP_z*(hip_height-hip_des)-kD_z*(hip_height_dot));
F_v = -kD_x*(y(4)-v);
F_th = (1/(y(3)-x_foot))*(kP_th*y(5) + kD_th*y(6));

B_y = B_y + F_hip + F_th;
B_x = B_x + F_v;

dydt = zeros(6,1);
dydt(1) = y(2);
dydt(2) = (B_y/m)-g;
dydt(3) = y(4);
dydt(4) = (B_x/m);
dydt(5) = y(6);
dydt(6) = (B_y*(x_foot - y(3)) + B_x*(y(1)))/I;
end

