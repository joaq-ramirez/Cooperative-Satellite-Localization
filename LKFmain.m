%Housekeeping
clear
clc
close all
%% Part 1
rng(100)
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

x = x0+perturb_x0;

P = 2*pi*sqrt(6678^3/mu);
figure
for j = 1:12
    theta0 = (j-1)*pi/6;
    y = [];
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
       end
    end
    subplot(3,1,1)
    plot(tvec,y(1,:),'o')
    hold on
    xlabel('Time (s)')
    ylabel('$\rho (km)$','Interpreter','latex')
    subplot(3,1,2)
    plot(tvec,y(2,:),'o')
    hold on
    xlabel('Time (s)')
    ylabel('$\dot{\rho} (km)$','Interpreter','latex')
    subplot(3,1,3)
    plot(tvec,y(3,:),'o')
    hold on
    xlabel('Time (s)')
    ylabel('$\phi (rad)$','Interpreter','latex')
    sgtitle('Linearized Output')
end


figure
subplot(4,1,1)
plot(tvec,x(1,1:length(tvec)))
xlabel('Time (s)')
ylabel('X (km)')
subplot(4,1,2)
plot(tvec,x(2,1:length(tvec)))
xlabel('Time (s)')
ylabel('$\dot{X}$ (km)','Interpreter','latex')
subplot(4,1,3)
plot(tvec,x(3,1:length(tvec)))
xlabel('Time (s)')
ylabel('Y (km)')
subplot(4,1,4)
plot(tvec,x(4,1:length(tvec)))
xlabel('Time (s)')
ylabel('$\dot{Y}$ (km)','Interpreter','latex')
sgtitle('Linearized States')

figure
subplot(4,1,1)
plot(tvec,perturb_x0(1,1:length(tvec)))
xlabel('Time (s)')
ylabel('\delta X (km)')
subplot(4,1,2)
plot(tvec,perturb_x0(2,1:length(tvec)))
xlabel('Time (s)')
ylabel('$\delta \dot{X}$ (km)','Interpreter','latex')
subplot(4,1,3)
plot(tvec,perturb_x0(3,1:length(tvec)))
xlabel('Time (s)')
ylabel('\delta Y (km)')
subplot(4,1,4)
plot(tvec,perturb_x0(4,1:length(tvec)))
xlabel('Time (s)')
ylabel('$\delta \dot{Y}$ (km)','Interpreter','latex')
sgtitle('Linearized Perturbations')

%% 3

dT = 10;

tspan = [0 1400*dT];

mu = 398600;

r0 = 6678;

x0 = [r0; 0; 0; r0 * sqrt(mu/r0^3)];

dx0 = [0;0.075;0;-0.021];

opts = odeset('reltol',1e-12,'abstol',1e-12);

[t,x] = ode45(@(t,x) NLODE(t,x,mu),tspan,x0+dx0,opts);

figure()
sgtitle('Nonlinear Simulation States')
subplot(4,1,1)
hold on
plot(t,x(:,1))
ylabel('$X [km]$','Interpreter','latex')
hold off
subplot(4,1,2)
hold on
plot(t,x(:,2))
ylabel('$\dot{X} [\frac{km}{s}]$ ','Interpreter','latex')
hold off
subplot(4,1,3)
hold on
plot(t,x(:,3))
ylabel('$Y [km]$','Interpreter','latex')
hold off
subplot(4,1,4)
hold on
plot(t,x(:,4))
xlabel('Time [s]')
ylabel('$\dot{Y} [\frac{km}{s}]$','Interpreter','latex')
hold off

tmax = length(t);

RE = 6378;

wE = 2*pi/86400;

for i = 1:12

    thetas0(i) = (i-1)*pi/6;

    for k = 1:tmax
        tk = t(k);
        Xs = RE*cos(wE*tk+thetas0(i));
        Ys = RE*sin(wE*tk+thetas0(i));
        Xsd = -RE*wE*sin(wE*tk+thetas0(i));
        Ysd = RE*wE*cos(wE*tk+thetas0(i));
        thetas(i,k) = atan2(Ys,Xs);
        rho(i,k) = sqrt((x(k,1)-Xs)^2 + (x(k,3)-Ys)^2);
        rhod(i,k) = ((x(k,1)-Xs)*(x(k,2)-Xsd)+(x(k,3)-Ys)*(x(k,4)-Ysd))/rho(i,k);
        phi(i,k) = atan2(x(k,3)-Ys,x(k,1)-Xs);

        if abs(angdiff(thetas(i,k),phi(i,k))) > pi/2
            rho(i,k) = nan;
            rhod(i,k) = nan;
            phi(i,k) = nan;
        end


    end
end

figure()
sgtitle('Nonlinear Simulation Outputs')
subplot(3,1,1)
hold on
for i = 1:12
    scatter(t,rho(i,:),'x')
end
ylabel('$\rho [km]$','Interpreter','latex')
hold off
subplot(3,1,2)
hold on
for i = 1:12
    scatter(t,rhod(i,:))
end
ylabel('$\dot{\rho} [\frac{km}{s}]$','Interpreter','latex')
hold off
subplot(3,1,3)
hold on
for i = 1:12
    scatter(t,phi(i,:))
end
ylabel('$\phi [rad]$','Interpreter','latex')
xlabel('Time [s]')
hold off


%% Part II STOCHASTIC NONLINEAR FILTER
%data vars: Qtrue, Rtrue, tvec, ydata

%Important Notes:
%No input u
dt = delT;

% Qk = [0.02 0; 0 0.02];
Qk = Qtrue;


Omk = dt*[0 0 ; 1 0; 0 0 ; 0 1];
Sw = chol(Qk,'lower');
Sv = chol(Rtrue,'lower');

% Monte Carlo Sim

N = 10;
n = 2;

alpha = 0.05;

for m = 1:N

    x_k = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';

%Noisy Data sim of x
for i = 1:length(tvec)

    %Simulate truth
    w_k = mvnrnd(zeros(2,1),Qk); % of Qtrue
    t_k = tvec(i);

    [~,x_nom] = ode45(@(t,x) odeLKF(t_k,x,[0;0]),[0 dt],x_k);

    x_nom = x_nom'; 
    
    x_k = x_nom(:,end);

    xnomp(:,i) = x_nom(:,end);

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

[x_lkf,P_lkf,NEES(m,:),NIS(m,:)] = LKF(x_t,y_t,station,tvec);

xdiff = x_t - x_lkf;

end

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
sgtitle('LKF State Estimation Errors')
for i = 1:4
    subplot(4,1,i)
    hold on
    plot(tvec, xdiff(i,:))
    sigma = 2*sqrt(P_lkf(i,i,:));
    sigma = reshape(sigma,1,1401);
    plot(tvec, sigma,'r--')
    plot(tvec, -sigma,'r--')
    hold off
end

n = 4;

r1 = chi2inv(alpha/2,N*n)/N;
r2 = chi2inv(1-alpha/2,N*n)/N;

figure
hold on
title('NEES')
yline(r1,'r--')
yline(r2,'r--')
scatter(tvec,NEESmean,'b')
hold off

n = 3;

r1 = chi2inv(alpha/2,N*n)/N;
r2 = chi2inv(1-alpha/2,N*n)/N;

figure
hold on
title('NIS')
yline(r1,'r--')
yline(r2,'r--')
scatter(tvec,NISmean,'b')
hold off

%% 6

[x_lkf,P_lkf,NEES,NIS] = LKFdata(x_t,station,tvec);

xdiff = x_t(:,2:1401) - x_lkf;


figure
sgtitle('LKF output')
for i = 1:4   
    subplot(4,1,i)
    hold on
    plot(tvec(2:1401), x_lkf(i,:))
    hold off
end


figure
sgtitle('LKF State Estimation Errors')
for i = 1:4
    subplot(4,1,i)
    hold on
    plot(tvec(2:1401), xdiff(i,:))
    sigma = 2*sqrt(P_lkf(i,i,:));
    sigma = reshape(sigma,1,1400);
    plot(tvec(2:1401), sigma,'r--')
    plot(tvec(2:1401), -sigma,'r--')
    hold off
end