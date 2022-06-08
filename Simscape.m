%% Sagittal Bounding Model
clear; close all; clc;

% Notes
% Tau_H for first back swing with opposite friction forces is -0.00763

%% Old lines of code
% Tried this on run 4:
%     if tau_H_mod < 0.25
%         tau_H_mod = 0.25;
%     end
% Tried this on run 13:
%      if tau_H_mod < (5*mod_time)
%          tau_H_mod = 5*mod_time;
%      end
% Tried this on run 15:
%     if tau_H_mod < (7*mod_time)
%         tau_H_mod = 7*mod_time;
%     end

% After Airborne phase
% T_swing_new = T_stance_new + T_air_est_H + T_air_est_F;


%% Constant Parameters
% System Parameters
m = 33;             % [kg] mass of hip
I = 2.9;            % [kg m^2]
r_1 = 0.4;          % [m] thigh length
r_2 = 0.3;          % [m] shank length
body_length = 0.7;  % [m]

% Environment Parameters
g = 9.80665;                               % [m/s^2] gravitational constant
contact_stiffness = 400/0.001;              % Approximated at weight (N) / desired displacement (m)
contact_damping = contact_stiffness/10;     % Tuned based on contact stiffness value
mu_s = 0.7;  %0.7                                 % Static friction coefficient: Around that of rubber-asphalt (actually made it so slipping was impossible)
mu_k = 0.7; %0.6                                 % Kinetic friction coefficient: Lower than the static coefficient
mu_vth = 0.1;                               % Friction velocity threshold (m/s)
TransRegion = 0.8e-4;

% Animation Parameters
dis = 3;            % [m] distance of travel for animation
start = 1;          % [] animation data start point
%% Foot Trajectory
Beta_x = [-0.2,-0.259,-0.275,-0.384, 0.261,-0.017, 0.248, 0.267, 0.259, 0.2];
Beta_y = [-0.5,-0.45,-0.406,-0.065,-1.031,0.095,-0.545,-0.374,-0.45,-0.5];
beta = [0, (1/3)*-4/(-0.531), 1-((1/3)*-4/-0.531), 1];

%% Initial Conditions
theta_1i = deg2rad(-101.5172); %-34.88                                                 % [rad] intial theta 1 angle
theta_2i = deg2rad(80.4059);   %-80.4059                                               % [rad] intial theta 2 angle
th_dot_i = 0;  %-1.15                                                             % [rad/s] initial angular velocity of the mass
th_i = -0.08;  %-0.07                                                             % [rad] initial mass angle
x_hip = -r_1*cos(theta_1i) - r_2*cos(theta_1i + theta_2i);                  % [m] initial x hip position
y_hip = -r_1*sin(theta_1i) - r_2*sin(theta_1i + theta_2i);                  % [m] initial y hip position

x_i = -0.2 - (body_length/2)*cos(th_i);
y_i = 0.5 -(body_length/2)*sin(th_i); 
x_dot_i = 4;                                                                % [m/s] initial center of mass horizontal velocity
y_dot_i = 0;                                                                % [m/s] initial center of mass vertical velocity
v = 4;

%potential back stance ICs
% theta_1i_B = deg2rad(-34.88);
% theta_2i_B = deg2rad(-80.4059);
theta_1i_B = theta_1i;
theta_2i_B = theta_2i;
%theta_1i_B = deg2rad(-78.5);
%theta_2i_B = deg2rad(-80.4);

%% Gait timing
l = 0.4;                            % [m] stride length
T_swing = 0.22;                     % [s] swing time
T_stance = l/v;                     % [s] stance time
T = T_swing + T_stance;             % [s] Total gait time
amp_y = m*g*T*pi()/(4*T_stance);   % [N] Vertical force profile amplitude
% amp_y_mod = 1.1; %1.1
% amp_y = amp_y*amp_y_mod;
amp_x = 50;                         % [N] Suggested from MIT paper (50)
time = 1.15;                         % [s] Simulation time 2.05 for x_dot plot 1.15 for theta phase plot
T_air = (T_swing - T_stance)/2;     % [s] Time in the air

%% Back IC
offset = 0.0295;
tau_F = (T_air+T_stance+offset)/T_swing;
tau_F_mod = bezier(beta,tau_F);
x_foot_des_F = bezier(Beta_x,tau_F_mod);
y_foot_des_F = bezier(Beta_y,tau_F);
[theta_1i, theta_2i] = inverse_kine(theta_1i, theta_2i, x_foot_des_F, y_foot_des_F, r_1, r_2);
tau_H = T_air/T_swing;
tau_H_mod = bezier(beta,tau_H);
x_foot_des_H = bezier(Beta_x,tau_H_mod);
y_foot_des_H = bezier(Beta_y,tau_H);
[theta_1i_B, theta_2i_B] = inverse_kine(theta_1i_B, theta_2i_B, x_foot_des_H, y_foot_des_H, r_1, r_2);
%% Impedance Control Gains
K = 500;    % Stiffness
B = 150;    % Damping
C = 30;
D = 15;

% MIT Controller
kP_z = 800; % 1000 without k_time %800 is other number
kD_z = 120; % 120
kD_x = 60; %60 %80 %100
kP_th = 25; %30 %25
kD_th = 15; %15 %11
% kP_th = 0;
% kD_th = 0;
k_time = 0.3; % 0.3

hip_des = 0.48; %0.48
% th_des = -0.03; %-0.07
% th_dot_des = -0.15; %-0.5
th_des = 0; %-0.07
th_dot_des = 0.15; %-0.5

% Flight Controller
kD_flt = 224.2605; %121.15
kP_flt = 1731.4; %457.563
N_flt = 4.3619;%36.2
kI_flt = -0.012; % -0.012

kI_knee = 303.1846;
kD_knee = 0;%1.4519e-14
kP_knee = 1.1527e+03;
N_knee = 0; %3.1840e-14

%System Damping
d_joint1 = 0.5;
d_joint2 = 0.2;
d_body = 0;

%% Simulation

u = [g,m,I,r_1,r_2,T_stance,T_air,amp_y,amp_x];

SSout = sim('Trying.slx');
% SSout = sim('T_stanceChangeRun18.slx');

%% Plots
% [A,~,idx]=unique(SSout.touch_time,'stable');
% n=accumarray(idx(:),1);
% vals_touch=A(find(n~=1));
% 
% t = 0:0.0005:time;
% t_first = vals_touch(2):0.0005:(T_stance+vals_touch(2));
% B_Y_1 = amp_y*sin((pi()/T_stance)*(t_first-vals_touch(2)));
% B_x_1 = amp_x*sin((pi()/T_stance)*(t_first-vals_touch(2)));
% 
% % t_second = vals_touch(3):0.0005:(T_stance+vals_touch(3));
% % B_Y_2 = amp_y*sin((pi()/T_stance)*(t_second-vals_touch(3)));
% % B_x_2 = amp_x*sin((pi()/T_stance)*(t_second-vals_touch(3)));
% 
% figure;
% plot(SSout.tout,SSout.normal);
% hold on;
% plot(t_first,B_Y_1);
% plot(SSout.tout, SSout.friction);
% plot(t_first,B_x_1);
% legend("Normal Force", "Ideal Normal Force", "Friction Force","Ideal Friction Force");
% hold off;
% 
figure;
hold on;
plot(SSout.th(1),SSout.th_dot(1),"b","Marker",".",'MarkerSize',20,'LineStyle','none');
plot(SSout.th(end),SSout.th_dot(end),"g","Marker",".",'MarkerSize',20,'LineStyle','none');
plot(SSout.th,SSout.th_dot,'k');
set(gca,'FontSize',14)
xlabel("Angular Position, [rad]","FontSize",18);
ylabel("Angular Velocity, [rad/s]","FontSize",18);
legend("Start","End","","Location","best");
hold off;

figure;
plot(SSout.tout,SSout.x_y_dot(:,1),'k');
set(gca,'FontSize',14)
% set(gca,'XTick',0:0.35:1.4)
% set(gca,'YTick',-0.5:1:4.5)
xlabel("Time, [s]","FontSize",18);
ylabel("Horizontal Velocity, [m/s]","FontSize",18);
% 
% figure;
% hold on;
% plot(SSout.tout,SSout.H_hip_torque);
% plot(SSout.tout, SSout.th_1_H_err);
% xlabel("Time [s]");
% ylabel("Airborne hip angle error and torque");
% legend("Torque", "Angle Error", 'Location','best');
% hold off;
% 
% figure;
% hold on;
% plot(SSout.tout,SSout.F_hip_torque);
% plot(SSout.tout, SSout.th_1_F_err);
% xlabel("Time [s]");
% ylabel("Airborne hip angle error and torque");
% legend("Torque", "Angle Error", 'Location','best');
% hold off;

% 
% t = 0:0.0005:(T_stance+SSout.touch_time(end));
% T_mod = (pi()/T_stance);
% 
% [A,~,idx]=unique(SSout.x_o,'stable');
% n=accumarray(idx(:),1);
% vals_x=A(find(n~=1));
% x_o = vals_x(2);
% 
% [A,~,idx]=unique(SSout.y_o,'stable');
% n=accumarray(idx(:),1);
% vals_y=A(find(n~=1));
% y_o = vals_y(1);
% 
% y_dot_o =0;
% 
% x_des = x_o + (t-SSout.touch_time(end))*x_dot_i - (amp_x*(t-SSout.touch_time(end)))/(T_mod*m) + (amp_x*sin(T_mod*(t-SSout.touch_time(end))))/(T_mod^2*m);
% y_des = y_o + ((T_mod*(-g)*(t-SSout.touch_time(end)).^2)/2 + T_mod*y_dot_o*(t-SSout.touch_time(end)))/T_mod + (amp_y*(t-SSout.touch_time(end)) - (amp_y.*sin(T_mod.*(t-SSout.touch_time(end))))/T_mod)/(T_mod*m);
% 
% figure;
% min1 = min(x_des);
% min2 = min(SSout.x_y(:,1));
% y_axis(1,1) = min([min1,min2]);
% max1 = max(x_des);
% max2 = max(SSout.x_y(:,1));
% y_axis(1,2) = max([max1,max2]);
% hold on;
% plot(t,x_des);
% plot(SSout.tout,SSout.x_y(:,1));
% plot([SSout.touch_time(end),SSout.touch_time(end)],[-5,5]);
% axis([0,t(end),y_axis]);
% legend("x des", "x act",'Location','best');
% hold off;
% 
% x_dot_des = (amp_x*cos(T_mod*(t-SSout.touch_time(end))) - amp_x + T_mod*m*x_dot_i)/(T_mod*m);
% y_dot_des = (amp_y - amp_y*cos(T_mod*(t-SSout.touch_time(end))) + T_mod*m*y_dot_o + T_mod*(-g)*m*(t-SSout.touch_time(end)))/(T_mod*m);
% figure;
% min1 = min(x_dot_des);
% min2 = min(SSout.x_y_dot(:,1));
% y_axis(1,1) = min([min1,min2]);
% max1 = max(x_dot_des);
% max2 = max(SSout.x_y_dot(:,1));
% y_axis(1,2) = max([max1,max2]);
% hold on;
% plot(t,x_dot_des);
% plot(SSout.tout,SSout.x_y_dot(:,1));
% plot([SSout.touch_time(end),SSout.touch_time(end)],[-5,5]);
% axis([0,t(end),y_axis]);
% legend("x dot des","x dot act",'Location','best');



