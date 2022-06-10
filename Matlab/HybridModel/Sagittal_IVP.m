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
% th_dot_i = -3.0207;%-3.0161    % [rad/s] initial angular velocity of the mass
th_dot_i = -3.02;
% th_i = -0.1662;%-0.1014        % [rad] initial mass angle
th_i = -0.16;
% x_i = -0.5734;%-0.5444              % [m] initial horizontal mass position
x_i = 0;
y_i = 0.5227;%0.5184               % [m] initial vertical mass position
x_dot_i = 4;                % [m/s] initial center of mass horizontal velocity
% y_dot_i = -0.2942;%-0.2943        % [m/s] initial center of mass vertical velocity
y_dot_i = 0;
v = 4;                      % [m/s] nominal horizontal velocity
x_foot = x_i + 0.2 + cos(th_i)*(body_length/2);

%% Gait timing
l = 0.4;                % [m] stride length
T_swing = 0.22;         % [s]
T_stance = l/v;         % [s]
T = T_swing + T_stance; % [s]
amp_y = m*g*T*pi()/(4*T_stance);
amp_x = 50; % Suggested from MIT paper (50)
% time = 0.24; 
T_air = (T_swing - T_stance)/2;

%% Controller
kP_z = 1000; % 1000 without k_time %800 is other number
kD_z = 120; % 120
kD_x = 60; %60 %80 %100
kP_th = 30; %30 %25
kD_th = 15; %15 %11
k_vals = [kP_z, kD_z, kD_x, kP_th, kD_th];
% k_vals = zeros(1,5);
hip_des = 0.48;
%% Simulation
% u = [g,m,I,T_stance,T_air,amp_y,amp_x,v,x_foot,hip_des];
% ic = [y_i,y_dot_i,x_i,x_dot_i,th_i,th_dot_i];
ic = [0.491966776748697,-0.297349049165308,25.559537067905340,4.048125270901569,-0.091333588311559,-3.027146704483003];
u = [9.810000000000000,33,2.900000000000000,0.100000000000000,0.060000000000000,8.136222317972990e+02,50,4,26.108078263157600,0.480000000000000]; 
n = 5;
options = odeset('MaxStep',1e-4);
figure;
%% FrontStance
for i = 1:n

% y = ode45(@(t,y) FrontStance(t,y,u),[0,T_stance],ic,options);
y = ode15s(@(t,y) FrontStance(t,y,u,k_vals),[0,T_stance],ic,options);

%% Airborne

ic_1 = [y.y(1,end),y.y(2,end),y.y(3,end),y.y(4,end),y.y(5,end),y.y(6,end)];

% x = ode45(@(t,y) airborne(t,y,u),[T_stance,T_stance+T_air],ic_1,options);
x = ode15s(@(t,y) airborne(t,y,u),[T_stance,T_stance+T_air],ic_1,options);

%% BackStance

%x_foot = x.y(3,end)+ r_1*cos(theta_1i) + r_2*cos(theta_1i + theta_2i) - (body_length/2)*cos(x.y(5,end));
% x_foot = x.y(3,end)+ r_1*cos(theta_1i_B) + r_2*cos(theta_1i_B + theta_2i_B) - (body_length/2)*cos(x.y(5,end));

x_foot = x.y(3,end) + 0.2 - (body_length/2)*cos(x.y(5,end));

u(9) = x_foot;

ic_2 = [x.y(1,end),x.y(2,end),x.y(3,end),x.y(4,end),x.y(5,end),x.y(6,end)];

% z = ode45(@(t,y) BackStance(t,y,u),[T_stance+T_air,T - T_air], ic_2, options);
z = ode15s(@(t,y) BackStance(t,y,u,k_vals),[T_stance+T_air,T - T_air], ic_2, options);

%% Airborne

ic_3 = [z.y(1,end),z.y(2,end),z.y(3,end),z.y(4,end),z.y(5,end),z.y(6,end)];

% w = ode45(@(t,y) airborne(t,y,u),[T-T_air,T],ic_3, options);
w = ode15s(@(t,y) airborne(t,y,u),[T-T_air,T],ic_3, options);

%% Plot of whole cycle
ic = [w.y(1,end),w.y(2,end),w.y(3,end),w.y(4,end),w.y(5,end),w.y(6,end)];
u(9) = w.y(3,end) + 0.2 + cos(w.y(5,end))*(body_length/2);
hold on;
% plot(y.y(5,:),y.y(6,:));
% plot(x.y(5,:),x.y(6,:));
% plot(z.y(5,:),z.y(6,:));
% plot(w.y(5,:),w.y(6,:));

% th_plt = [y.y(5,:),x.y(5,:),z.y(5,:),w.y(5,:)];
% th_dot_plt = [y.y(6,:),x.y(6,:),z.y(6,:),w.y(6,:)];
% plot(th_plt,th_dot_plt,'k');
y_plt = [y.y(1,:),x.y(1,:),z.y(1,:),w.y(1,:)];
y_dot_plt = [y.y(2,:),x.y(2,:),z.y(2,:),w.y(2,:)];
plot(y_plt,y_dot_plt,'k');

end
%% Orbit Plot
% plot(th_i,th_dot_i,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(y.y(5,end),y.y(6,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(x.y(5,end),x.y(6,end),"r", "Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(z.y(5,end),z.y(6,end),"y", "Marker",".",'MarkerSize',20,'LineStyle','none');
% set(gca,'FontSize',14)
% xlabel("Angular Position [rad]","FontSize",18);
% ylabel("Angular Velocity [rad/s]","FontSize",18);
% legend("","","","","q_0","q_1","q_2","q_3","Location","best");
% hold off;

%% Open Loop theta phase plot n = 4
% plot(th_i,th_dot_i,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(w.y(5,end),w.y(6,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
% L((n)+1) = "start";
% L((n)+2) = "end";
% set(gca,'FontSize',14)
% xlabel("Angular Position [rad]","FontSize",18);
% ylabel("Angular Velocity [rad/s]","FontSize",18);
% set(gca,'YTick',-3:1.5:3)
% legend(L,"Location","best");
% axis([-0.3,0.2,-4,3])
% hold off;

%% Open Loop y phase plot n = 4
% plot(y_i,y_dot_i,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(w.y(1,end),w.y(2,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
% L((n)+1) = "start";
% L((n)+2) = "end";
% set(gca,'FontSize',14)
% ylabel("Mass Vertical Velocity [m/s]","FontSize",18);
% xlabel("Mass Height [m]","FontSize",18);
% set(gca,'XTick',0.5:0.01:0.53)
% axis([0.505,0.53,-0.4,0.4])
% legend(L,"Location","best");
% hold off;

%% Open Loop Unstable  n = 15
% plot(th_i,th_dot_i,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(w.y(5,end),w.y(6,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
% L((n)+1) = "start";
% L((n)+2) = "end";
% set(gca,'FontSize',14)
% xlabel("Angular Position [rad]","FontSize",18);
% ylabel("Angular Velocity [rad/s]","FontSize",18);
% set(gca,'YTick',-3:1.5:3)
% legend(L,"Location","best");
% hold off;

%% Closed Loop theta phase n = 25
% plot(th_i,th_dot_i,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(w.y(5,end),w.y(6,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
% L((n)+1) = "start";
% L((n)+2) = "end";
% set(gca,'FontSize',14)
% xlabel("Angular Position [rad]","FontSize",18);
% ylabel("Angular Velocity [rad/s]","FontSize",18);
% set(gca,'YTick',-3:1.5:3)
% legend(L,"Location","best");
% hold off;

%% Closed Loop y phase n = 25
% plot(y_plt,y_dot_plt,'k', "LineWidth",4);
% plot(y_i,y_dot_i,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
% plot(w.y(1,end),w.y(2,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
% L((n+1)) = "Stable Orbit";
% L((n)+2) = "start";
% L((n)+3) = "end";
% set(gca,'FontSize',14)
% ylabel("Mass Vertical Velocity [m/s]","FontSize",18);
% xlabel("Mass Height [m]","FontSize",18);
% set(gca,'YTick',-0.8:0.35:0.6)
% legend(L,"Location","best");
% hold off;

%% Closed Loop y phase n = 21-25
plot(0.491966776748697,-0.297349049165308,"b","Marker",".",'MarkerSize',20,'LineStyle','none');
plot(w.y(1,end),w.y(2,end),"g", "Marker",".",'MarkerSize',20,'LineStyle','none');
L((n)+1) = "start";
L((n)+2) = "end";
set(gca,'FontSize',14)
ylabel("Mass Vertical Velocity [m/s]","FontSize",18);
xlabel("Mass Height [m]","FontSize",18);
set(gca,'XTick',0.47:0.01:0.5)
axis([0.47,0.5,-0.5,0.5])
legend(L,"Location","best");
hold off;

