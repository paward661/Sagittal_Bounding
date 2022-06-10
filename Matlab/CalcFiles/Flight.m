function x = Flight(t,eta,u)
y_o = eta(1);
y_dot_o = eta(2);
x_o = eta(3);
x_dot_o = eta(4);
th_o = eta(5);
th_dot_o = eta(6);
g = u(1);
t_o = u(7);

x = [y_o + y_dot_o*(t - t_o) - (g*(t - t_o)^2)/2;
    y_dot_o - g*t + g*t_o;
    x_o + x_dot_o*(t - t_o);
    x_dot_o;
    th_o + th_dot_o*(t - t_o);
    th_dot_o];
end

