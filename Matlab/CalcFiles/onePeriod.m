function x = onePeriod(eta,u)
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
amp_y = u(5);
amp_x = u(6);
t_o = u(7);
body_length = u(8);
T_air = u(9);

x_foot = 0.2 + (body_length/2)*cos(th_o) + x_o;
T_mod = pi()/T_stance;
t = T_stance;

%% Front Stance
x_front = [y_o + y_dot_o.*(t - t_o) + (2*amp_y*sin(T_mod*t_o) + 2*T_mod*amp_y*t*cos(T_mod*t_o) - 2*T_mod*amp_y*t_o*cos(T_mod*t_o) - T_mod^2*g*m*t_o^2 + 2*T_mod^2*g*m*t*t_o)/(2*T_mod^2*m) - (2*amp_y*sin(T_mod*t) + T_mod^2*g*m*t^2)/(2*T_mod^2*m);
    (amp_y*cos(T_mod.*t_o) - amp_y*cos(T_mod*t) + T_mod.*m*y_dot_o - T_mod*g*m*t + T_mod*g*m*t_o)/(T_mod*m);
    (amp_x*sin(T_mod*t) - amp_x.*sin(T_mod*t_o) + T_mod^2.*m*x_o - T_mod*amp_x*t*cos(T_mod*t_o) + T_mod*amp_x*t_o*cos(T_mod*t_o) + T_mod^2*m*t*x_dot_o - T_mod^2*m*t_o*x_dot_o)/(T_mod^2*m);
    (amp_x*cos(T_mod*t) - amp_x.*cos(T_mod*t_o) + T_mod*m*x_dot_o)/(T_mod*m); 
    th_o - (sin(T_mod*t).*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o - 2*T_mod*amp_x*y_o - (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o + 2*T_mod*amp_x*t_o*y_dot_o + T_mod*amp_x*g*t_o^2) - cos(T_mod*t)*(4*amp_y*x_dot_o + 4*amp_x*y_dot_o + 4*amp_x*g*t_o) - t*(2*I*T_mod^3*th_dot_o - 2*amp_x*g*cos(T_mod*t_o) + 2*T_mod*amp_y*x_dot_o*sin(T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_mod*t_o) - 2*T_mod^2*amp_x*y_o*cos(T_mod*t_o)) + 4*amp_x*g*t*cos(T_mod*t) - 2*T_mod*t*sin(T_mod*t)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o) + T_mod*amp_x*g*t^2*sin(T_mod*t))/(2*I*T_mod^3) + (sin(T_mod*t_o)*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o - 2*T_mod*amp_x*y_o - (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o + 2*T_mod*amp_x*t_o*y_dot_o + T_mod*amp_x*g*t_o^2) - cos(T_mod*t_o)*(4*amp_y*x_dot_o + 4*amp_x*y_dot_o + 4*amp_x*g*t_o) - t_o*(2*I*T_mod^3*th_dot_o - 2*amp_x*g*cos(T_mod*t_o) + 2*T_mod*amp_y*x_dot_o*sin(T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_mod*t_o) - 2*T_mod^2*amp_x*y_o*cos(T_mod*t_o)) + 4*amp_x*g*t_o*cos(T_mod*t_o) - 2*T_mod*t_o*sin(T_mod*t_o)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o) + T_mod*amp_x*g*t_o^2*sin(T_mod*t_o))/(2*I*T_mod^3);
    th_dot_o - (cos(T_mod*t)*((2*amp_y*x_foot - 2*amp_y*x_o - 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o + 2*amp_x*t_o*y_dot_o + amp_x*g*t_o^2)/T_mod - (2*amp_x*g)/T_mod^3))/(2*I) + (cos(T_mod*t_o)*((2*amp_y*x_foot - 2*amp_y*x_o - 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o + 2*amp_x*t_o*y_dot_o + amp_x*g*t_o^2)/T_mod - (2*amp_x*g)/T_mod^3))/(2*I) - (sin(T_mod*t)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod^2) + (sin(T_mod*t_o)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod^2) + (t*cos(T_mod*t)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod) - (t_o*cos(T_mod*t_o)*(amp_y*x_dot_o + amp_x*y_dot_o + amp_x*g*t_o))/(I*T_mod) + (amp_x*g*t*sin(T_mod*t))/(I*T_mod^2) - (amp_x*g*t_o*sin(T_mod*t_o))/(I*T_mod^2) - (amp_x*g*t^2*cos(T_mod*t))/(2*I*T_mod) + (amp_x*g*t_o^2*cos(T_mod*t_o))/(2*I*T_mod)];

%% Flight 1
y_o = x_front(1);
y_dot_o = x_front(2);
x_o = x_front(3);
x_dot_o = x_front(4);
th_o = x_front(5);
th_dot_o = x_front(6);

t_o = T_stance;
t = T_stance + T_air;

x_flight_one = [y_o + y_dot_o.*(t - t_o) - (g*(t - t_o)^2)/2;
    y_dot_o - g*t + g.*t_o;
    x_o + x_dot_o.*(t - t_o);
    x_dot_o;
    th_o + th_dot_o.*(t - t_o);
    th_dot_o];

%% Back Stance
y_o = x_flight_one(1);
y_dot_o = x_flight_one(2);
x_o = x_flight_one(3);
x_dot_o = x_flight_one(4);
th_o = x_flight_one(5);
th_dot_o = x_flight_one(6);

t_o = T_stance + T_air;
t = (2.*T_stance) + T_air;
x_foot = 0.2 - (body_length/2)*cos(th_o) + x_o;

x_back = [y_o + y_dot_o.*(t - t_o) - (g.*t^2)/2 + (g*t_o*(2*t - t_o))/2 + (amp_y*sin(T_mod*(T_air + T_stance - t)))/(T_mod^2*m) - (amp_y*(sin(T_mod*(T_air + T_stance - t_o)) - T_mod*t*cos(T_mod*(T_air + T_stance - t_o)) + T_mod*t_o*cos(T_mod*(T_air + T_stance - t_o))))/(T_mod^2*m);
    (amp_y*cos(T_mod.*(T_air + T_stance - t_o)) - amp_y*cos(T_mod*(T_air + T_stance - t)) + T_mod*m*y_dot_o - T_mod*g*m*t + T_mod*g*m*t_o)/(T_mod*m); 
    x_o + x_dot_o.*(t - t_o) - (amp_x*sin(T_mod*(T_air + T_stance - t_o)) - T_mod*amp_x*(t*cos(T_mod*(T_air + T_stance - t_o)) - t_o*cos(T_mod*(T_air + T_stance - t_o))))/(T_mod^2*m) + (amp_x*sin(T_mod*(T_air + T_stance - t)))/(T_mod^2*m);
    (amp_x*cos(T_mod.*(T_air + T_stance - t_o)) - amp_x*cos(T_mod*(T_air + T_stance - t)) + T_mod*m*x_dot_o)/(T_mod*m);
    th_o - (cos(T_mod.*(T_air + T_stance - t))*(4*amp_x*y_dot_o - 4*amp_y*x_dot_o + 4*amp_x*g*t_o))/(2*I*T_mod^3) + (cos(T_mod*(T_air + T_stance - t_o))*(4*amp_x*y_dot_o - 4*amp_y*x_dot_o + 4*amp_x*g*t_o))/(2*I*T_mod^3) + (t*(2*amp_x*g*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*I*T_mod^3*th_dot_o - 2*T_mod*amp_y*x_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_x*y_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o)))/(2*I*T_mod^3) - (t_o*(2*amp_x*g*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*I*T_mod^3*th_dot_o - 2*T_mod*amp_y*x_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod*amp_x*y_dot_o*sin(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_y*x_foot*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) - 2*T_mod^2*amp_y*x_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o) + 2*T_mod^2*amp_x*y_o*cos(T_air*T_mod + T_mod*T_stance - T_mod*t_o)))/(2*I*T_mod^3) + (sin(T_mod*(T_air + T_stance - t))*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o + 2*T_mod*amp_x*y_o + (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o - 2*T_mod*amp_x*t_o*y_dot_o - T_mod*amp_x*g*t_o^2))/(2*I*T_mod^3) - (sin(T_mod*(T_air + T_stance - t_o))*(2*T_mod*amp_y*x_foot - 2*T_mod*amp_y*x_o + 2*T_mod*amp_x*y_o + (6*amp_x*g)/T_mod + 2*T_mod*amp_y*t_o*x_dot_o - 2*T_mod*amp_x*t_o*y_dot_o - T_mod*amp_x*g*t_o^2))/(2*I*T_mod^3) + (t*sin(T_mod*(T_air + T_stance - t))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) - (t_o*sin(T_mod*(T_air + T_stance - t_o))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) - (amp_x*g*t^2*sin(T_mod*(T_air + T_stance - t)))/(2*I*T_mod^2) + (amp_x*g*t_o^2*sin(T_mod*(T_air + T_stance - t_o)))/(2*I*T_mod^2) + (2*amp_x*g*t*cos(T_mod*(T_air + T_stance - t)))/(I*T_mod^3) - (2*amp_x*g*t_o*cos(T_mod*(T_air + T_stance - t_o)))/(I*T_mod^3);
    th_dot_o - (cos(T_mod.*(T_air + T_stance - t))*((2*amp_y*x_foot - 2*amp_y*x_o + 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o - 2*amp_x*t_o*y_dot_o - amp_x*g*t_o^2)/T_mod + (2*amp_x*g)/T_mod^3))/(2*I) + (cos(T_mod*(T_air + T_stance - t_o))*((2*amp_y*x_foot - 2*amp_y*x_o + 2*amp_x*y_o + 2*amp_y*t_o*x_dot_o - 2*amp_x*t_o*y_dot_o - amp_x*g*t_o^2)/T_mod + (2*amp_x*g)/T_mod^3))/(2*I) - (sin(T_mod*(T_air + T_stance - t))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) + (sin(T_mod*(T_air + T_stance - t_o))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod^2) - (t*cos(T_mod*(T_air + T_stance - t))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod) + (t_o*cos(T_mod*(T_air + T_stance - t_o))*(amp_x*y_dot_o - amp_y*x_dot_o + amp_x*g*t_o))/(I*T_mod) + (amp_x*g*t^2*cos(T_mod*(T_air + T_stance - t)))/(2*I*T_mod) - (amp_x*g*t_o^2*cos(T_mod*(T_air + T_stance - t_o)))/(2*I*T_mod) + (amp_x*g*t*sin(T_mod*(T_air + T_stance - t)))/(I*T_mod^2) - (amp_x*g*t_o*sin(T_mod*(T_air + T_stance - t_o)))/(I*T_mod^2)];

%% Flight 2
y_o = x_back(1);
y_dot_o = x_back(2);
x_o = x_back(3);
x_dot_o = x_back(4);
th_o = x_back(5);
th_dot_o = x_back(6);

t_o = (2*T_stance) + T_air;
t = (2*T_stance) + (2*T_air);

x = [y_o + y_dot_o.*(t - t_o) - (g*(t - t_o)^2)/2, y_dot_o - g.*t + g*t_o, x_o + x_dot_o.*(t - t_o),x_dot_o,th_o + th_dot_o*(t - t_o),th_dot_o];
% x = double(x_flight);

end

