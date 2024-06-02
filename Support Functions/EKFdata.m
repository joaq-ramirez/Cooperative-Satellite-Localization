function [x,Pptot,NEES,NIS] = EKFdata(x0,Pp,tvec,station,Qk,Rp,ystat,xNom)
    
    load('orbitdeterm_finalproj_KFdata.mat')
    Omegak = 10 * [0 0;
        1 0;
        0 0;
        0 1];
    x = nan(length(x0),length(tvec));
    x(:,1) = x0;
    Pptot = [];
    for i = 1:length(tvec)-1
        y_d = ydata{i+1};
        
        if any(any(isnan(y_d))) || any(any(isempty(y_d)))
%             yk = ystat(:,:,station(i));
%             yk = yk(:,i);
            stationCur = station(i);
        else
            stationCur = y_d(4);
            yk = y_d(1:3,1);
        end
        [~,xode] = ode45(@(t,xODE) odeLKF(tvec(i),xODE,[0 0]),[0 10],x(:,i));
        x_m = xode(end,:);
        [Fk,~] = linearizeMat(x(:,i),10);
        Pm = Fk*Pp*Fk' + Omegak*Qk*Omegak';
        Hp = linearizeH(i,x_m,tvec,stationCur);
        Kp = Pm*Hp' * inv(Hp*Pm*Hp'+Rp);
        [ycalc, ~, ~] = dynamicsCalc2(x_m,stationCur,tvec(i));
    %         if isnan(ycalc)
    %             y_m = yk(:,i);
    %         else
    %             y_m = ycalc;
    %         end
        y_m = ycalc;
        ey = yk - y_m;
        x_p = x_m' + Kp*ey;
        x(:,i+1) = x_p;
        Pp = (eye(length(x0)) - Kp *Hp)*Pm;
        Pptot(:,:,i) = Pp;
    
        xdiff = xNom(:,i) - x_p;
    
        NEES(i) = xdiff' * inv(Pp) * xdiff;
        NIS(i) = ey' * inv(Hp*Pm*Hp'+Rp) * ey;
    i
    end
end