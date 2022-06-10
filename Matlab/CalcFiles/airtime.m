function [T_air_mod] = airtime(theta_init,theta_dot,y_init,y_dot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
l = 0.7;
g = 9.81;
p = [-0.5*g y_dot y_init-0.5];
t_guess = roots(p)
if t_guess(1) > 0
    t = t_guess(1);
else
    t = t_guess(2);
end
% T_air = 0.26;
% t = T_air;
epsilon = 1;
n = 0;
figure;
hold on;
plot(n,t,"s","MarkerSize",12);
while (epsilon > 0.000001)
    n = n+1;
    f = -0.5*g*(t^2) + y_dot*t + y_init - 0.5*l*sin(theta_dot*t + theta_init) - 0.5;
    f_prime = -g*t + y_dot - 0.5*l*theta_dot*cos(theta_dot*t+theta_init);
    t_n = t - f/f_prime;
    epsilon = abs(t_n-t);
    t = t_n;
    plot(n,t,"s","MarkerSize",12);
end
T_air_mod = t;
end

