function x = Front(t,eta,u)
y_o = eta(1);
y_dot_o = eta(2);
x_o = eta(3);
x_dot_o = eta(4);
th_o = eta(5);
th_dot_o = eta(6);
g = u(1);
m = u(2);
I = u(3);
T_mod = pi()/u(4);
amp_y = u(5);
amp_x = u(6);
t_o = u(7);
% body_length = u(8);
x_foot = u(10);

% x = [y_o + y_dot_o*(t - t_o) - (2*amp_y*sin(T_mod*t) - 2*T_mod*amp_y*t*cos(T_mod*t) + 2*T_mod*amp_y*t_o*cos(T_mod*t) - T_mod^2*g*m*t^2 + 2*T_mod^2*g*m*t*t_o)/(2*T_mod^2*m) + (2*amp_y*sin(T_mod*t_o) + T_mod^2*g*m*t_o^2)/(2*T_mod^2*m);
%     (amp_y*cos(T_mod*t_o) - amp_y*cos(T_mod*t) + T_mod*m*y_dot_o - T_mod*g*m*t + T_mod*g*m*t_o)/(T_mod*m);
%     (amp_x*sin(T_mod*t) - amp_x*sin(T_mod*t_o) + T_mod^2*m*x_o - T_mod*amp_x*t*cos(T_mod*t) + T_mod*amp_x*t_o*cos(T_mod*t) + T_mod^2*m*t*x_dot_o - T_mod^2*m*t_o*x_dot_o)/(T_mod^2*m);
%     (amp_x*cos(T_mod*t) - amp_x*cos(T_mod*t_o) + T_mod*m*x_dot_o)/(T_mod*m);
%     th_o + th_dot_o*(t - t_o) - ((sin(T_mod*t)*(- amp_x*g*T_mod^2*t_o^2 + 6*amp_y*x_foot*T_mod^2 + 2*amp_x*g))/T_mod^3 - (2*t*cos(T_mod*t)*(3*amp_y*x_foot*T_mod^2 + amp_x*g))/T_mod^2 + (amp_x*g*t^3*cos(T_mod*t_o))/3 + (2*t_o*cos(T_mod*t)*(3*amp_y*x_foot*T_mod^2 + amp_x*g))/T_mod^2 - (amp_x*g*t^2*sin(T_mod*t))/T_mod + amp_x*g*t*t_o^2*cos(T_mod*t_o) - amp_x*g*t^2*t_o*cos(T_mod*t_o) + (2*amp_x*g*t*t_o*sin(T_mod*t))/T_mod)/(6*I*T_mod) + (6*amp_x*g*sin(T_mod*t_o) + 18*T_mod^2*amp_y*x_foot*sin(T_mod*t_o) + T_mod^3*amp_x*g*t_o^3*cos(T_mod*t_o))/(18*I*T_mod^4) + (amp_y*x_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod) + (amp_x*y_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod) + (amp_y*x_dot_o*(t - t_o)^2*(cos(T_mod*t) - cos(T_mod*t_o)))/(6*I*T_mod) + (amp_x*y_dot_o*(t - t_o)^2*(cos(T_mod*t) - cos(T_mod*t_o)))/(6*I*T_mod);
%     th_dot_o + (cos(T_mod*t_o)*(2*amp_y*x_foot + (amp_x*g)/T_mod^2) + (amp_x*g*t_o^2*cos(T_mod*t_o))/2)/(2*I*T_mod) - (cos(T_mod*t)*(2*amp_y*x_foot + (amp_x*g)/T_mod^2) - (amp_x*g*t^2*cos(T_mod*t_o))/2 + amp_x*g*t*t_o*cos(T_mod*t_o) + (amp_x*g*t*sin(T_mod*t))/T_mod - (amp_x*g*t_o*sin(T_mod*t))/T_mod)/(2*I*T_mod) + (amp_y*x_o*(cos(T_mod*t) - cos(T_mod*t_o)))/(I*T_mod) + (amp_x*y_o*(cos(T_mod*t) - cos(T_mod*t_o)))/(I*T_mod) + (amp_y*x_dot_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod) + (amp_x*y_dot_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod)];

% x = [y_o + y_dot_o*(t - t_o) - (2*amp_y*sin(T_mod*t) - 2*T_mod*amp_y*t*cos(T_mod*t) + 2*T_mod*amp_y*t_o*cos(T_mod*t) - T_mod^2*g*m*t^2 + 2*T_mod^2*g*m*t*t_o)/(2*T_mod^2*m) + (2*amp_y*sin(T_mod*t_o) + T_mod^2*g*m*t_o^2)/(2*T_mod^2*m);
%     (amp_y*cos(T_mod*t_o) - amp_y*cos(T_mod*t) + T_mod*m*y_dot_o - T_mod*g*m*t + T_mod*g*m*t_o)/(T_mod*m);
%     (amp_x*sin(T_mod*t) - amp_x*sin(T_mod*t_o) + T_mod^2*m*x_o - T_mod*amp_x*t*cos(T_mod*t) + T_mod*amp_x*t_o*cos(T_mod*t) + T_mod^2*m*t*x_dot_o - T_mod^2*m*t_o*x_dot_o)/(T_mod^2*m);
%     (amp_x*cos(T_mod*t) - amp_x*cos(T_mod*t_o) + T_mod*m*x_dot_o)/(T_mod*m);
%     th_o + th_dot_o*(t - t_o) - ((sin(T_mod*t)*(10*amp_x*g + 6*T_mod^2*amp_y + 30*T_mod^2*amp_y*x_o - 5*T_mod^2*amp_x*g*t_o^2))/T_mod^3 + (5*amp_x*g*t^3*cos(T_mod*t_o))/3 - (2*t*cos(T_mod*t)*(5*amp_x*g + 3*T_mod^2*amp_y + 15*T_mod^2*amp_y*x_o))/T_mod^2 + (2*t_o*cos(T_mod*t)*(5*amp_x*g + 3*T_mod^2*amp_y + 15*T_mod^2*amp_y*x_o))/T_mod^2 - (5*amp_x*g*t^2*sin(T_mod*t))/T_mod + 5*amp_x*g*t*t_o^2*cos(T_mod*t_o) - 5*amp_x*g*t^2*t_o*cos(T_mod*t_o) + (10*amp_x*g*t*t_o*sin(T_mod*t))/T_mod)/(30*I*T_mod) + (18*T_mod^2*amp_y*sin(T_mod*t_o) + 30*amp_x*g*sin(T_mod*t_o) + 90*T_mod^2*amp_y*x_o*sin(T_mod*t_o) + 5*T_mod^3*amp_x*g*t_o^3*cos(T_mod*t_o))/(90*I*T_mod^4) + (amp_y*x_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod) + (amp_x*y_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod) + (amp_y*x_dot_o*(t - t_o)^2*(cos(T_mod*t) - cos(T_mod*t_o)))/(6*I*T_mod) + (amp_x*y_dot_o*(t - t_o)^2*(cos(T_mod*t) - cos(T_mod*t_o)))/(6*I*T_mod); 
%     th_dot_o - ((cos(T_mod*t)*(5*amp_x*g + 2*T_mod^2*amp_y + 10*T_mod^2*amp_y*x_o))/T_mod^2 - (5*amp_x*g*t^2*cos(T_mod*t_o))/2 + 5*amp_x*g*t*t_o*cos(T_mod*t_o) + (5*amp_x*g*t*sin(T_mod*t))/T_mod - (5*amp_x*g*t_o*sin(T_mod*t))/T_mod)/(10*I*T_mod) + (cos(T_mod*t_o)*(10*amp_x*g + 4*T_mod^2*amp_y + 20*T_mod^2*amp_y*x_o + 5*T_mod^2*amp_x*g*t_o^2))/(20*I*T_mod^3) + (amp_y*x_o*(cos(T_mod*t) - cos(T_mod*t_o)))/(I*T_mod) + (amp_x*y_o*(cos(T_mod*t) - cos(T_mod*t_o)))/(I*T_mod) + (amp_y*x_dot_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod) + (amp_x*y_dot_o*(t - t_o)*(cos(T_mod*t) - cos(T_mod*t_o)))/(2*I*T_mod)];

x = [y_o + y_dot_o*(t - t_o) + (2*amp_y*sin(T_mod*t_o) + 2*T_mod*amp_y*t*cos(T_mod*t_o) - 2*T_mod*amp_y*t_o*cos(T_mod*t_o) - T_mod^2*g*m*t_o^2 + 2*T_mod^2*g*m*t*t_o)/(2*T_mod^2*m) - (2*amp_y*sin(T_mod*t) + T_mod^2*g*m*t^2)/(2*T_mod^2*m);
    (amp_y*cos(T_mod*t_o) - amp_y*cos(T_mod*t) + T_mod*m*y_dot_o - T_mod*g*m*t + T_mod*g*m*t_o)/(T_mod*m);
    (amp_x*sin(T_mod*t) - amp_x*sin(T_mod*t_o) + T_mod^2*m*x_o - T_mod*amp_x*t*cos(T_mod*t_o) + T_mod*amp_x*t_o*cos(T_mod*t_o) + T_mod^2*m*t*x_dot_o - T_mod^2*m*t_o*x_dot_o)/(T_mod^2*m);
    (amp_x*cos(T_mod*t) - amp_x*cos(T_mod*t_o) + T_mod*m*x_dot_o)/(T_mod*m); 
    th_o - (sin(T_mod*t)*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o - 2*T_mod*amp_x*y_o - (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o + 2*T_mod*amp_x*t_o*y_dot_o + T_mod*amp_x*g*t_o^2) - cos(T_mod*t)*(4*amp_y*x_dot_o + 4*amp_x*y_dot_o + 4*amp_x*g*t_o) - t*(2*I*T_mod^3*th_dot_o - 2*amp_x*g*cos(T_mod*t_o) + 2*T_mod*amp_y*x_dot_o*sin(T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_mod*t_o) - 2*T_mod^2*amp_x*y_o*cos(T_mod*t_o)) + 4*amp_x*g*t*cos(T_mod*t) - 2*T_mod*t*sin(T_mod*t)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o) + T_mod*amp_x*g*t^2*sin(T_mod*t))/(2*I*T_mod^3) + (sin(T_mod*t_o)*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o - 2*T_mod*amp_x*y_o - (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o + 2*T_mod*amp_x*t_o*y_dot_o + T_mod*amp_x*g*t_o^2) - cos(T_mod*t_o)*(4*amp_y*x_dot_o + 4*amp_x*y_dot_o + 4*amp_x*g*t_o) - t_o*(2*I*T_mod^3*th_dot_o - 2*amp_x*g*cos(T_mod*t_o) + 2*T_mod*amp_y*x_dot_o*sin(T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_mod*t_o) - 2*T_mod^2*amp_x*y_o*cos(T_mod*t_o)) + 4*amp_x*g*t_o*cos(T_mod*t_o) - 2*T_mod*t_o*sin(T_mod*t_o)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o) + T_mod*amp_x*g*t_o^2*sin(T_mod*t_o))/(2*I*T_mod^3);
    th_dot_o - (cos(T_mod*t)*((2*amp_y*x_foot - 2*amp_y*x_o - 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o + 2*amp_x*t_o*y_dot_o + amp_x*g*t_o^2)/T_mod - (2*amp_x*g)/T_mod^3))/(sym(2)*I) + (cos(T_mod*t_o)*((2*amp_y*x_foot - 2*amp_y*x_o - 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o + 2*amp_x*t_o*y_dot_o + amp_x*g*t_o^2)/T_mod - (2*amp_x*g)/T_mod^3))/(sym(2)*I) - (sin(T_mod*t)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod^2) + (sin(T_mod*t_o)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod^2) + (t*cos(T_mod*t)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod) - (t_o*cos(T_mod*t_o)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod) + (amp_x*g*t*sin(T_mod*t))/(I*T_mod^2) - (amp_x*g*t_o*sin(T_mod*t_o))/(I*T_mod^2) - (amp_x*g*t^2*cos(T_mod*t))/(2*I*T_mod) + (amp_x*g*t_o^2*cos(T_mod*t_o))/(2*I*T_mod)];
end

