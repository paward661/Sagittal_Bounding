function [theta_1 , theta_2] = inverse_kine(theta_1g , theta_2g , foot_x_inv , foot_y_inv,r_1,r_2)

% inverse_kine Inverse Kinematic Function
%   Takes two input guesses and returns the actual values values based off
%   of reducing error. Theta_2g and Theta_3g are the last calculated values
%   of theta 2 and theta 3 respectively that we will use as a point to
%   start iterating from. r_1 is the known magnitude vector 1 and theta_1
%   is the angle of that vector.

%% Inverse Kinematics for In Contact

% Initial Guesses

th_1 = theta_1g; % Theta 1
th_2 = theta_2g; % Theta 2

% Our reasonable tolerance for position
epsilon = 5e-5;

% Setting a base epsilon to insure that code is read
epsilon1 = 1;
epsilon2 = 1;

% Setting a run-out timer
runs = 0;

%% Looping to converge on answer

while (epsilon1>= epsilon || epsilon2 >= epsilon && runs <= 10)
    
    runs = runs +1;
    
    % Defining Sin and Cos for our angles
    s12 = sin(th_2 + th_1);
    s1 = sin(th_1);
    c12 = cos(th_2 + th_1);
    c1 = cos(th_1);
    
    %Equations of Motion based off our angles
    delt_x = r_1*c1 + r_2*c12 - foot_x_inv;
    delt_y = r_1*s1 + r_2*s12 - foot_y_inv;
    
    %Partial derivatives
    dx_dth1 =-r_1*s1-r_2*s12;
    dx_dth2 =-s12*r_2;
    dy_dth1 =r_1*c1+r_2*c12;
    dy_dth2 =c12*r_2;
    
    %Defining Jacobian
    Jac = [dx_dth1 , dx_dth2 ; dy_dth1 , dy_dth2];
    
    %Setting old thetas and storing as "o"
    th_1o = th_1;
    th_2o = th_2;
    
    %Reguessing 
    xn_yn = [th_1;th_2] - ((Jac)\[delt_x;delt_y]);
    th_1 = xn_yn(1);
    th_2 = xn_yn(2);
    
    %Subtracting to refine our guesses
    epsilon1 = abs(th_1 - th_1o);
    epsilon2 = abs(th_2 - th_2o);
    


end
if (runs >= 10)
    
    error("inverse kinematics failed");
end

theta_1 = th_1;
theta_2 = th_2;

end