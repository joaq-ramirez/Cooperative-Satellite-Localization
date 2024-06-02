% House keeping
clc
clear all
close all

%% Extended Kalman Filter
load('orbitdeterm_finalproj_KFdata.mat')
rE = 6378;
mu = 4*10^5;
omega_e = 2*pi/86400;

perturb_x0 = [0,0.075,0,-0.021]';
x0 = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';

B = [0 0;
    1 0;
    0 0;
    0 1];
Gamma = [0 0;
        1 0;
        0 0;
        0 1];
delT = 10;

x = [];
    
P = 2*pi*sqrt(6678^3/mu);
y_tf = [];
y_tot = nan(3,length(tvec));
station = zeros(1,length(tvec));
x_tot = [];
for j = 1:12
    theta0 = (j-1)*pi/6;
    y = [];
    x_tot = [x_tot; x];
    x = x0+perturb_x0(:,1);
    for i = 1:length(tvec)
        x_state = x(:,i);
        x1 = 6678*cos(2*pi/P * tvec(i));
        x2 = -6678*sin(2*pi/P * tvec(i));
        x3 = 6678*sin(2*pi/P * tvec(i));
        x4 = 6678*cos(2*pi/P * tvec(i));
        [F,G] = linearizeMat([x1 x2 x3 x4]', delT);
        
        perturb_x0_new = F*perturb_x0(:,i);
    %     x_state(1) = x1;
    %     x_state(3) = x3;
        x_state_new = [x1 x2 x3 x4]' + perturb_x0_new;
        x = [x x_state_new];

%         x1 = x_state(1);
%         x2 = x_state(2);
%         x3 = x_state(3);
%         x4 = x_state(4);
        
        Xi = rE * cos(omega_e * tvec(i) + theta0);
        Yi = rE * sin(omega_e * tvec(i) + theta0);
        Xi_dot = -rE * omega_e * sin(omega_e * tvec(i) + theta0);
        Yi_dot = rE * omega_e * cos(omega_e * tvec(i) + theta0);
        
        rho = sqrt((x1 - Xi)^2 + (x3 - Yi)^2);
        rho_dot = ((x1 - Xi)*(x2 - Xi_dot) + (x3 - Yi)*(x4 - Yi_dot)) / rho;
        phi = atan2((x3 - Yi),(x1 - Xi));
        
        H = [(x1-Xi)/rho 0 (x3-Yi)/rho 0;
            (x2-Xi_dot)/rho+(Xi-x1)*rho_dot/rho^2 (x1-Xi)/rho (x4-Yi_dot)/rho+(Yi-x3)*rho_dot/rho^2 (x3-Yi)/rho;
            (Yi-x3)/((x1-Xi)^2 + (x3-Yi)^2) 0 (x1-Xi)/((x1-Xi)^2 + (x3-Yi)^2) 0];
        perturb_y = H*perturb_x0(:,i);
        ystate = [rho rho_dot phi]' + perturb_y;
        y = [y ystate];

        perturb_x0 = [perturb_x0 perturb_x0_new];
        
        thetai = atan2(Yi,Xi);
        if abs(angdiff(thetai,phi)) > pi/2
           y(1,i) = nan;
           y(2,i) = nan;
           y(3,i) = nan;
        else
            station(i) = j;
        end
    end
    repy = isnan(y_tot);
    y_tot(repy) = y(repy);
    y_tf(:,:,j) = y;
end

for i = 1:length(tvec)
    [~,xode] = ode45(@(t,xODE) odeLKF(tvec(i),xODE,[0 0]),[0 10],x(:,i));
end

P0 = diag([0.01, 0.001, 0.01, 0.001]);

Qk = 10^-7*Qtrue;
x_k = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';
Omk = 10*[0 0 ; 1 0; 0 0 ; 0 1];
Sw = chol(Qk,'lower');
Sv = chol(Rtrue,'lower');
% Monte Carlo Sim

N = 10;
n = 2;

alpha = 0.05;

for m = 1:N

    x_k = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';

for i = 1:length(tvec)

    %Simulate truth
    w_k = mvnrnd(zeros(2,1),Qk); % of Qtrue
    t_k = tvec(i);

    [~,x_nom] = ode45(@(t,x) odeLKF(t_k,x,[0;0]),[0 10],x_k);
    
    x_k = x_nom(end,:)';

    xnomp(:,i) = x_nom(end,:)';

    w = Sw*randn(2,1);
    x_t(:,i) = xnomp(:,i) + Omk*w; % Truth solution
    %Simulate measurement truth
    for j = 1:12
        v = Sv*randn(3,1);
        y_t(:,i,j) = dynamicsCalc(x_t(:,i),j,t_k)+v;
        if any(~isnan(y_t(:,i,j)))
            station(:,i,1) = j;            
        end
    end
end
    [x_state,Pp,NEES(m,:),NIS(m,:)] = EKF(x0,P0,tvec,station,Qk,Rtrue,y_t,xnomp);
    xdiff = x_t - x_state(:,2:end);

end
%% NEES/NIS Stat Calc



figure
for i = 1:4
    subplot(4,1,i)
    hold on
    plot(tvec, x_state(i,2:end))
    hold off
end
sgtitle('EKF Simulated States')
 
NEESmean = mean(NEES);
NISmean = mean(NIS);

%Plotting noiseless and noisey signal
figure
sgtitle('Noisy Simulated Truth States')
for i = 1:4   
    subplot(4,1,i)
    hold on
    plot(tvec, x_t(i,:))
    hold off
end

figure()
sgtitle('Noisy Simulated Data')
subplot(3,1,1)
hold on
for i = 1:12
    scatter(tvec,y_t(1,:,i),'x')
end
ylabel('$\rho [km]$','Interpreter','latex')
hold off
subplot(3,1,2)
hold on
for i = 1:12
    scatter(tvec,y_t(2,:,i))
end
ylabel('$\dot{\rho} [\frac{km}{s}]$','Interpreter','latex')
hold off
subplot(3,1,3)
hold on
for i = 1:12
    scatter(tvec,y_t(3,:,i))
end
ylabel('$\phi [rad]$','Interpreter','latex')
xlabel('Time [s]')
hold off

figure
sgtitle('EKF State Estimation Errors')
for i = 1:4
    subplot(4,1,i)
    hold on
    plot(tvec, xdiff(i,:))
    sigma = 2*sqrt(Pp(i,i,:));
    sigma = reshape(sigma,1,1401);
    plot(tvec, sigma,'r--')
    plot(tvec, -sigma,'r--')
    hold off
end

r1 = chi2inv(alpha/2,N*4)/N;
r2 = chi2inv(1-alpha/2,N*4)/N;
figure
hold on
title('EKF NEES')
yline(r1,'r--')
yline(r2,'r--')
scatter(tvec,NEESmean,'b')
hold off

r1 = chi2inv(alpha/2,N*3)/N;
r2 = chi2inv(1-alpha/2,N*3)/N;
figure
hold on
title('EKF NIS')
scatter(tvec,NISmean,'b')
yline(r1,'r--')
yline(r2,'r--')
hold off

%% 6

[x_statedat,Pp,~,~] = EKFdata(x0,P0,tvec,station,Qtrue,Rtrue,y_t,x_t);
xdiff = x_t - x_statedat;

figure
sgtitle('EKF State Estimation Errors with data')
for i = 1:4
    subplot(4,1,i)
    hold on
    plot(tvec, x_statedat(i,:))
    plot(tvec,x_t(i,:),'-')
%     sigma = 2*sqrt(Pp(i,i,:));
%     sigma = reshape(sigma,1,1401);
%     plot(tvec, sigma,'r--')
%     plot(tvec, -sigma,'r--')
    hold off
end
sgtitle('EKF State Simulation with Data')

figure
for i = 1:4
    subplot(4,1,i)
    hold on
    plot(tvec, xdiff(i,:))
    sigma = 2*sqrt(Pp(i,i,:));
    sigma = reshape(sigma,1,1400);
    plot(tvec(1:1400), sigma,'r--')
    plot(tvec(1:1400), -sigma,'r--')
    hold off
end
sgtitle('EKF State Estimation Errors with Data')