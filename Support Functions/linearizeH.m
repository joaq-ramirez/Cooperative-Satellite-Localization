function [H] = linearizeH(i,x_vals,tvec,stat)
%Initialize
theta0 = (stat-1)*pi/6;% determine theta
rE = 6378;
mu = 4*10^5;
omega_e = 2*pi/86400;

%Allocate xvalues for H calculation
x1 = x_vals(1);
x2 = x_vals(2);
x3 = x_vals(3);
x4 = x_vals(4);

%Equations to determine a stations distance from the sat
Xi = rE * cos(omega_e * tvec(i) + theta0);
Yi = rE * sin(omega_e * tvec(i) + theta0);
Xi_dot = -rE * omega_e * sin(omega_e * tvec(i) + theta0);
Yi_dot = rE * omega_e * cos(omega_e * tvec(i) + theta0);

%More things to calculate
rho = sqrt((x1 - Xi)^2 + (x3 - Yi)^2);
rho_dot = ((x1 - Xi)*(x2 - Xi_dot) + (x3 - Yi)*(x4 - Yi_dot)) / rho;
phi = atan2((x3 - Yi),(x1 - Xi));

%Finally the output
H = [(x1-Xi)/rho 0 (x3-Yi)/rho 0;
            (x2-Xi_dot)/rho+(Xi-x1)*rho_dot/rho^2 (x1-Xi)/rho (x4-Yi_dot)/rho+(Yi-x3)*rho_dot/rho^2 (x3-Yi)/rho;
            (Yi-x3)/((x1-Xi)^2 + (x3-Yi)^2) 0 (x1-Xi)/((x1-Xi)^2 + (x3-Yi)^2) 0];

end

