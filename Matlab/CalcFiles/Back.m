function x = Back(t,eta,u)
y_o = eta(1);
y_dot_o = eta(2);
x_o = eta(3);
x_dot_o = eta(4);
th_o = eta(5);
th_dot_o = eta(6);
g = u(1);
m = u(2);
I = u(3);
T_stance = u(4);
T_mod = pi()/T_stance;
amp_y = u(5);
amp_x = u(6);
t_o = u(7);
body_length = u(8);
T_air = u(9);

x_foot = 0.2 - (body_length/2)*cos(th_o) + x_o;

x = [y_o + y_dot_o*(t - t_o) - (g*t^2)/2 + (g*t_o*(2*t - t_o))/2 + (amp_y*sin(T_mod*(T_air + T_stance - t)))/(T_mod^2*m) - (amp_y*(sin(T_mod*(T_air + T_stance - t_o)) - T_mod*t*cos(T_mod*(T_air + T_stance - t_o)) + T_mod*t_o*cos(T_mod*(T_air + T_stance - t_o))))/(T_mod^2*m);
    (amp_y*cos(T_mod*(T_air + T_stance - t_o)) - amp_y*cos(T_mod*(T_air + T_stance - t)) + T_mod*m*y_dot_o - T_mod*g*m*t + T_mod*g*m*t_o)/(T_mod*m); 
    x_o + x_dot_o*(t - t_o) - (amp_x*sin(T_mod*(T_air + T_stance - t_o)) - T_mod*amp_x*(t*cos(T_mod*(T_air + T_stance - t_o)) - t_o*cos(T_mod*(T_air + T_stance - t_o))))/(T_mod^2*m) + (amp_x*sin(T_mod*(T_air + T_stance - t)))/(T_mod^2*m);
    (amp_x*cos(T_mod*(T_air + T_stance - t_o)) - amp_x*cos(T_mod*(T_air + T_stance - t)) + T_mod*m*x_dot_o)/(T_mod*m);
    th_o - (cos(T_mod*(T_air + T_stance - t))*(4*amp_x*y_dot_o - 4*amp_y*x_dot_o + 4*amp_x*g*t_o))/(2*I*T_mod^3) + (cos(T_mod*(T_air + T_stance - t_o))*(4*amp_x*y_dot_o - 4*amp_y*x_dot_o + 4*amp_x*g*t_o))/(2*I*T_mod^3) + (t*(2*amp_x*g*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*I*T_mod^3*th_dot_o - 2*T_mod*amp_y*x_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_x*y_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o)))/(2*I*T_mod^3) - (t_o*(2*amp_x*g*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*I*T_mod^3*th_dot_o - 2*T_mod*amp_y*x_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_x*y_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o)))/(2*I*T_mod^3) + (sin(T_mod*(T_air + T_stance - t))*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o + 2*T_mod*amp_x*y_o + (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o - 2*T_mod*amp_x*t_o*y_dot_o - T_mod*amp_x*g*t_o^2))/(2*I*T_mod^3) - (sin(T_mod*(T_air + T_stance - t_o))*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o + 2*T_mod*amp_x*y_o + (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o - 2*T_mod*amp_x*t_o*y_dot_o - T_mod*amp_x*g*t_o^2))/(2*I*T_mod^3) + (t*sin(T_mod*(T_air + T_stance - t))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) - (t_o*sin(T_mod*(T_air + T_stance - t_o))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) - (amp_x*g*t^2*sin(T_mod*(T_air + T_stance - t)))/(2*I*T_mod^2) + (amp_x*g*t_o^2*sin(T_mod*(T_air + T_stance - t_o)))/(2*I*T_mod^2) + (2*amp_x*g*t*cos(T_mod*(T_air + T_stance - t)))/(I*T_mod^3) - (2*amp_x*g*t_o*cos(T_mod*(T_air + T_stance - t_o)))/(I*T_mod^3);
    th_dot_o - (cos(T_mod*(T_air + T_stance - t))*((2*amp_y*x_foot - 2*amp_y*x_o + 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o - 2*amp_x*t_o*y_dot_o - amp_x*g*t_o^2)/T_mod + (2*amp_x*g)/T_mod^3))/(sym(2)*I) + (cos(T_mod*(T_air + T_stance - t_o))*((2*amp_y*x_foot - 2*amp_y*x_o + 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o - 2*amp_x*t_o*y_dot_o - amp_x*g*t_o^2)/T_mod + (2*amp_x*g)/T_mod^3))/(sym(2)*I) - (sin(T_mod*(T_air + T_stance - t))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) + (sin(T_mod*(T_air + T_stance - t_o))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) - (t*cos(T_mod*(T_air + T_stance - t))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod) + (t_o*cos(T_mod*(T_air + T_stance - t_o))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod) + (amp_x*g*t^2*cos(T_mod*(T_air + T_stance - t)))/(2*I*T_mod) - (amp_x*g*t_o^2*cos(T_mod*(T_air + T_stance - t_o)))/(2*I*T_mod) + (amp_x*g*t*sin(T_mod*(T_air + T_stance - t)))/(I*T_mod^2) - (amp_x*g*t_o*sin(T_mod*(T_air + T_stance - t_o)))/(I*T_mod^2)];
end

