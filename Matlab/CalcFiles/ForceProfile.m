%% Sagittal Bounding Model
clear;

%% Constant Parameters
% System Parameters
m = 33;             % [kg] mass of hip
I = 2.9;            % [kg m^2]
r_1 = 0.4;          % [m] thigh length
r_2 = 0.3;          % [m] shank length
body_length = 0.7;  % [m]

% Environment Parameters
g = 9.81;           % [m/s^2] gravitational constant

% Animation Parameters
dis = 3;            % [m] distance of travel for animation
start = 1;          % [] animation data start point

%% Initial Conditions
theta_1i = deg2rad(-34.88);                                                                % [rad] intial theta 1 angle
theta_2i = deg2rad(-80.4);                                                                % [rad] intial theta 2 angle
th_dot_i = deg2rad(-3);                                                                  % [rad/s] initial angular velocity of the mass
th_i = deg2rad(-7);                                                                      % [rad] initial mass angle
x_hip = -r_1*cos(theta_1i) - r_2*cos(theta_1i + theta_2i);
y_hip = -r_1*sin(theta_1i) - r_2*sin(theta_1i + theta_2i);
x_i = -r_1*cos(theta_1i) - r_2*cos(theta_1i + theta_2i) - (body_length/2)*cos(th_i);    % [m] initial x hip position
y_i = -r_1*sin(theta_1i) - r_2*sin(theta_1i + theta_2i) - (body_length/2)*sin(th_i);    % [m] initial y hip position
x_dot_i = 4;                                                                            % [m/s] initial center of mass horizontal velocity
y_dot_i = -0.2;                                                                            % [m/s] initial center of mass vertical velocity
v = 4;
theta_1i_B = -(pi()+theta_1i);
theta_2i_B = -theta_2i;
%theta_1i_B = deg2rad(-78.5);
%theta_2i_B = deg2rad(-80.4);

%% Gait timing
l = 0.4;               % [m] stride length
T_swing = 0.22;        % [s]
T_stance = l/v;  % [s]
T = T_swing + T_stance;               % [s]
amp_y = m*g*T*pi()/(4*T_stance);
amp_x = 50; % Suggested from MIT paper (50)
time = 0.24;
T_air = (T_swing - T_stance)/2;

%T_stance = l/x_dot     % [s] stance time
%T_swing = 0.22;     % [s] swing time

%% Force Profile
stance_array = 0:0.001:T_stance;
t= 0:0.001:T;
B_y_1(1:length(stance_array)) = amp_y*sin((pi()/T_stance)*t(1:length(stance_array)));
B_y_1(length(stance_array)+1:length(t)) = 0;

swing_array = 0:0.001:T_swing;
B_y_2(1:length(swing_array)) = 0;
B_y_2(length(swing_array)+1:length(t)) = amp_y*sin((pi()/T_stance)*(t(length(swing_array)+1:length(t))-T_swing));

figure;
hold on;
plot(t,B_y_1);
plot(t,B_y_2);
set(gca,'FontSize',14)
xlabel("Time [s]","FontSize",18);
ylabel("Normal Reaction Force [N]","FontSize",18);
