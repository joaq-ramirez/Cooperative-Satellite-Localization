function [fout] = odeLKF(tk,xk,wk)
mu = 4*10^5;
r = sqrt(xk(1)^2 + xk(3)^2);

dxdt = zeros(2,1);
dxdt(1) = xk(2);
dxdt(2) = -mu*xk(1)/(r^3) + wk(1);

dydt = zeros(2,1);
dydt(1) = xk(4);
dydt(2) = -mu*xk(3)/(r^3) + wk(2);

fout = [dxdt(1) ; dxdt(2) ; dydt(1) ; dydt(2)];

end

