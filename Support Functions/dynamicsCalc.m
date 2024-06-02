function [y_s,theta_i,phi_i] = dynamicsCalc(x_vals,stat,tvec)

%Initialize
theta0 = (stat-1)*pi/6; %determine theta
rE = 6378;
mu = 4*10^5;
omega_e = 2*pi/86400;

%Allocate xvalues for H calculation
x1 = x_vals(1);
x2 = x_vals(2);
x3 = x_vals(3);
x4 = x_vals(4);

%Equations to determine a stations distance from the sat
Xi = rE * cos(omega_e * tvec + theta0);
Yi = rE * sin(omega_e * tvec + theta0);
Xi_dot = -rE * omega_e * sin(omega_e * tvec + theta0);
Yi_dot = rE * omega_e * cos(omega_e * tvec + theta0);


%More things to calculate
rho = sqrt((x1 - Xi)^2 + (x3 - Yi)^2);
rho_dot = ((x1 - Xi)*(x2 - Xi_dot) + (x3 - Yi)*(x4 - Yi_dot)) / rho;
phi = atan2((x3 - Yi),(x1 - Xi));

y_s = [rho rho_dot phi]';

theta_i = atan2(Xi,Yi);
phi_i = atan2((x_vals(3)- Yi),(x_vals(1)- Xi));

thetai = atan2(Yi,Xi);
    if abs(angdiff(thetai,phi)) > pi/2
       y_s(1) = nan;
       y_s(2) = nan;
       y_s(3) = nan;
   end
end

