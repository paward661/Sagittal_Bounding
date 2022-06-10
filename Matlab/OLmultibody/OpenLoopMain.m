%% Sagittal Bounding Model
clear; close all; clc;

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

theta_1i_B = theta_1i;
theta_2i_B = theta_2i;

%% Gait timing
l = 0.4;                            % [m] stride length
T_swing = 0.22;                     % [s] swing time
T_stance = l/v;                     % [s] stance time
T = T_swing + T_stance;             % [s] Total gait time
amp_y = m*g*T*pi()/(4*T_stance);    % [N] Vertical force profile amplitude
amp_x = 50;                         % [N] Suggested from MIT paper (50)
time = 4*T;                         % [s] Simulation time (hits a double back stance at sixth cycle)
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

SSout = sim('OpenLoop.slx');

%% Plots
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
set(gca,'XTick',0:0.35:1.4)
set(gca,'YTick',-0.5:1:4.5)
xlabel("Time, [s]","FontSize",18);
ylabel("Horizontal Velocity, [m/s]","FontSize",18);

figure;
plot(SSout.tout,SSout.x_y(:,2),'k');
set(gca,'FontSize',14)
set(gca,'XTick',0:0.35:1.4)
set(gca,'YTick',0.46:0.05:0.66)
xlabel("Time, [s]","FontSize",18);
ylabel("Vertical Position, [m]","FontSize",18);
axis([0,1.4,0.46,0.66])


