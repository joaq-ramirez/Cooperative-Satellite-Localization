function [F,G] = linearizeMat(x_vals,dt)

x1 = x_vals(1);
x2 = x_vals(2);
x3 = x_vals(3);
x4 = x_vals(4);

rE = 6378;
mu = 4*10^5;
omega_e = 2*pi/86400;

A = [0 1 0 0;
            mu*(2*x1^2 -x3^2)/(x1^2 + x3^2)^(5/2) 0 3*mu*x1*x3/(x1^2 + x3^2)^(5/2) 0
            0 0 0 1;
            3*mu*x1*x3/(x1^2 + x3^2)^(5/2) 0 mu*(2*x3^2 -x1^2)/(x1^2 + x3^2)^(5/2) 0];
F = eye(4) + A*dt;

B = [0 0;
    1 0;
    0 0;
    0 1];

G = dt*B; 


end

