function [x_new,y_new] = LKF(F,H,x_nom0,dt)
%LKF Linearized Kalman Filter

%Calculate period for nominal state
mu = 4*10^5; 
P = 2*pi*sqrt(6678^3/mu);

for i = 1:length(tvec)
%Recalc state nominal state solution for time step
    x1 = 6678*cos(2*pi/P * tvec(i));
    x2 = -6678*sin(2*pi/P * tvec(i));
    x3 = 6678*sin(2*pi/P * tvec(i));
    x4 = 6678*cos(2*pi/P * tvec(i));
    
%Find updated F matrix and H matrix

    
    
perturb_x0p1 = F*
x_new = x_nom0 + perturb_x0p1


end


end

