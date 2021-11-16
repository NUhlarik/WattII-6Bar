close all
clear all
clc
%%
%Constants & Initial Guess

%Angles
Theta1 = 0; %deg
Theta2 = 0; %deg
Theta3star = 60; %deg
Theta4star = 110; %deg
Theta5star = 320; %deg
Theta6star = 40; %deg
Theta11 = 0; %deg

%Radii
R1 = 203.20/1000; %mm
R2 = 57.15/1000; %mm
R3 = 184.15/1000; %mm
R4 = 177.8/1000; %mm
R5 = 50.80/1000; %mm
R6 = 127.00/1000; %mm
R11 = 101.60/1000; %mm
R33 = (R3 / 2)/1000; %mm
R44 = 127.00/1000; %mm
R55 = 25.40/1000; %mm

% Masses
m2 = 2; %kg
m3 = 5.5; %kg
m4 = 7.5; %kg
m5 = 1.5; %kg
m6 = 6; %kg

% Second Moments of Mass
Ig2 = 0.0067; %kgm^2
Ig3 = 0.0433; %kgm^2
Ig4 = 0.2426; %kgm^2
Ig5 = 0.0009; %kgm^2
Ig6 = 0.0634; %kgm^2

% Balance Link Acceleration
A2x = 0;
A2y = 0;
A4x = 0;
A4y = 0;
A6x = 0;
A6y = 0;


g = 9.810; %m/s^2

resolution = 0.01; %deg

Omega2a = 25; %rad/s
Omega2b = 50; %rad/s
Alpha2 = 0; %rad/s^2

Theta3 = double.empty;
Theta4 = double.empty;
Theta5 = double.empty;
Theta6 = double.empty;

Det_1 = double.empty;
Det_2 = double.empty;

Theta3prime = double.empty;
Theta4prime = double.empty;
Theta5prime = double.empty;
Theta6prime = double.empty;

Theta3Doubleprime = double.empty;
Theta4Doubleprime = double.empty;
Theta5Doubleprime = double.empty;
Theta6Doubleprime = double.empty;

Xp = double.empty;
Yp = double.empty;

Xprime = double.empty;
Yprime = double.empty;
Xdprime = double.empty;
Ydprime = double.empty;
Rpprime = double.empty;

Xg = double.empty;
Yg = double.empty;
Xgprime = double.empty;
Ygprime = double.empty;
XgDoubleprime = double.empty;
YgDoubleprime = double.empty;

Utx = double.empty;
Uty = double.empty;
Unx = double.empty;
Uny = double.empty;
Rhop = double.empty;
Xcc = double.empty;
Ycc = double.empty;

F65x = double.empty;
F65y = double.empty;
F34x = double.empty;
F34y = double.empty;
F16x = double.empty;
F16y = double.empty;
F12x = double.empty;
F12y = double.empty;
F45x = double.empty;
F45y = double.empty;
F14x = double.empty;
F14y = double.empty;
F23x = double.empty;
F23y = double.empty;
T12 = double.empty;

DeltaTheta3star = 1;
DeltaTheta4star = 1;
DeltaTheta5star = 1;
DeltaTheta6star = 1;

%% Postures
for i = 1:360
    % Section 1 Newton-Raphson
    iteration = 0;
    while ((abs(DeltaTheta3star) > resolution) || (abs(DeltaTheta4star) > resolution))
        
        iteration = iteration + 1;
        
        epsilon_x = R2 .* cosd(Theta2) + R3 .* cosd(Theta3star) - R4 .* cosd(Theta4star) - R1 .* cosd(Theta1);
        epsilon_y = R2 .* sind(Theta2) + R3 .* sind(Theta3star) - R4 .* sind(Theta4star) - R1 .* sind(Theta1);
        
        DeltaTheta3star = (((epsilon_x .* R4 .* cosd(Theta4star)) + (epsilon_y .* R4 .* sind(Theta4star))) ./ (R3 .* R4 .* sind(Theta3star - Theta4star))) .* (180/pi); %deg
        DeltaTheta4star = (((epsilon_y .* R3 .* sind(Theta3star)) + (epsilon_x .* R3 .* cosd(Theta3star))) ./ (R3 .* R4 .* sind(Theta3star - Theta4star))) .* (180/pi); %deg
        
        Theta3star = Theta3star + DeltaTheta3star; %deg
        Theta4star = Theta4star + DeltaTheta4star; %deg
        
        if (iteration > 10)
            disp('Error Section 1')
            break
        end
    end
    
    Theta3(i) = Theta3star;
    Theta4(i) = Theta4star;
    
    DeltaTheta3star = 1;
    DeltaTheta4star = 1;
    
    Det_1(i) = R3 .* R4 .* sind(Theta3star - Theta4star) .* (180 / pi); %deg
    
    % Section 2 Newton-Raphson
    iteration = 0;
    Theta44 = Theta4star; %deg
    
    while (abs(DeltaTheta5star) > resolution || abs(DeltaTheta6star) > resolution)
        
        iteration = iteration + 1;
        
        epsilon_x2 = R6 .* cosd(Theta6star) - R5 .* cosd(Theta5star) - R44 .* cosd(Theta44) - R11 .* cosd(Theta11);
        epsilon_y2 = R6 .* sind(Theta6star) - R5 .* sind(Theta5star) - R44 .* sind(Theta44) - R11 .* sind(Theta11);
        
        DeltaTheta5star = ((((-epsilon_x2) .* R6 .* cosd(Theta6star)) + (-epsilon_y2 .* R6 .* sind(Theta6star))) ./ (R5 .* R6 .* sind(Theta5star - Theta6star))) .* (180 / pi); %deg
        DeltaTheta6star = ((((-epsilon_y2) .* R5 .* sind(Theta5star)) + (-epsilon_x2 .* R5 .* cosd(Theta6star))) ./ (R5 .* R6 .* sind(Theta5star - Theta6star))) .* (180 / pi);
        
        Theta5star = Theta5star + DeltaTheta5star; %deg
        Theta6star = Theta6star + DeltaTheta6star; %deg
        
        if (iteration > 10)
            disp('Error Section 2')
            break
        end
    end
    
    if (Theta5star > 360)
        Theta5star = Theta5star - 360;
        
    elseif (Theta5star < 0)
        Theta5star = Theta5star + 360;
    end
    
    Det_2(i) = R5 .* R6 .* sind(Theta5star - Theta6star); %deg
    
    Theta5(i) = Theta5star;
    Theta6(i) = Theta6star;
    
    DeltaTheta5star = 1;
    DeltaTheta6star = 1;
    
    Theta2 = Theta2 + 1;
end
%% Kinematic Coefficients of Links
Theta2 = 0;

for i = 1:360
    %Link 3 and 4
    A = [-R3 * sind(Theta3(i)), R4 * sind(Theta4(i)); R3 * cosd(Theta3(i)), -R4 * cosd(Theta4(i))];
    B = [R2 * sind(Theta2); -R2 * cosd(Theta2)];
    
    X = linsolve(A,B);
    Theta3prime(i) = X(1);
    Theta4prime(i) = X(2);
    
    B = [R2 * cosd(Theta2) + R3 * cosd(Theta3(i)) * Theta3prime(i)^2 - R4 * cosd(Theta4(i)) * Theta4prime(i)^2; R2 * sind(Theta2) + R3 * sind(Theta3(i)) * Theta3prime(i)^2 - R4 * sind(Theta4(i)) * Theta4prime(i)^2];
    X = linsolve(A,B);
    Theta3Doubleprime(i) = X(1);
    Theta4Doubleprime(i) = X(2);
    
    %Link 5 and 6
    A = [R5 * sind(Theta5(i)), -R6 * sind(Theta6(i)); -R5 * cosd(Theta5(i)), R6 * cosd(Theta6(i))];
    B = [-R44 * sind(Theta4(i)) * Theta4prime(i); R44 * cosd(Theta4(i)) * Theta4prime(i)];
    
    X = linsolve(A,B);
    Theta5prime(i) = X(1);
    Theta6prime(i) = X(2);
    
    B = [-R44 * (cosd(Theta4(i)) * Theta4prime(i)^2 + sind(Theta4(i)) * Theta4Doubleprime(i)) - R5 * cosd(Theta5(i)) * Theta5prime(i)^2 + R6 * cosd(Theta6(i)) * Theta6prime(i)^2; ...
        R44 * (-sind(Theta4(i)) * Theta4prime(i)^2 + cosd(Theta4(i)) * Theta4Doubleprime(i)) - R5 * sind(Theta5(i)) * Theta5prime(i)^2 + R6 * sind(Theta6(i)) * Theta6prime(i)^2];
    
    X = linsolve(A,B);
    Theta5Doubleprime(i) = X(1);
    Theta6Doubleprime(i) = X(2);
    
    Theta2 = Theta2 + 1;
end

%Velocity
Omega3a = Theta3prime * Omega2a;
Omega4a = Theta4prime * Omega2a;
Omega5a = Theta5prime * Omega2a;
Omega6a = Theta6prime * Omega2a;

Omega3b = Theta3prime * Omega2b;
Omega4b = Theta4prime * Omega2b;
Omega5b = Theta5prime * Omega2b;
Omega6b = Theta6prime * Omega2b;

%Acceleration
Alpha3a = Theta3Doubleprime * Omega2a^2 + Theta3prime * Alpha2;
Alpha4a = Theta4Doubleprime * Omega2a^2 + Theta4prime * Alpha2;
Alpha5a = Theta5Doubleprime * Omega2a^2 + Theta5prime * Alpha2;
Alpha6a = Theta6Doubleprime * Omega2a^2 + Theta6prime * Alpha2;

Alpha3b = Theta3Doubleprime * Omega2b^2 + Theta3prime * Alpha2;
Alpha4b = Theta4Doubleprime * Omega2b^2 + Theta4prime * Alpha2;
Alpha5b = Theta5Doubleprime * Omega2b^2 + Theta5prime * Alpha2;
Alpha6b = Theta6Doubleprime * Omega2b^2 + Theta6prime * Alpha2;

%% Position Analysis of Point P
Theta2 = 0;

for i = 1:360
    Xp(i) = R1 * cosd(Theta1) + R44 * cosd(Theta4(i)) + R55 * cosd(Theta5(i));
    Yp(i) = R1 * sind(Theta1) + R44 * sind(Theta4(i)) + R55 * sind(Theta5(i));
    
    % Kinematic Coefficients of Point P
    Xprime(i) = -R44 * sind(Theta4(i)) * Theta4prime(i) - R55 * sind(Theta5(i)) * Theta5prime(i);
    Yprime(i) = R44 * cosd(Theta4(i)) * Theta4prime(i) + R55 * cosd(Theta5(i)) * Theta5prime(i);
    
    Xdprime(i) = -R44 * (cosd(Theta4(i)) * Theta4prime(i)^2 + sind(Theta4(i)) * Theta4Doubleprime(i)) ...
        -R55 * (cosd(Theta5(i)) * Theta5prime(i)^2 + sind(Theta5(i)) * Theta5Doubleprime(i));
    Ydprime(i) = -R44 * (sind(Theta4(i)) * Theta4prime(i)^2 - cosd(Theta4(i)) * Theta4Doubleprime(i)) ...
        -R55 * (sind(Theta5(i)) * Theta5prime(i)^2 - cosd(Theta5(i)) * Theta5Doubleprime(i));
    
    % Kinematic Coefficients of Point G
    Xg(i) = R1 * cosd(Theta1) + R4 * cosd(Theta4(i)) - R33 * cosd(Theta3(i));
    Yg(i) = R1 * sind(Theta1) + R4 * sind(Theta4(i)) - R33 * sind(Theta3(i));
    
    Xgprime(i) = -R4 * sind(Theta4(i)) * Theta4prime(i) + R33 * sind(Theta3(i)) * Theta3prime(i);
    Ygprime(i) = R4 * cosd(Theta4(i)) * Theta4prime(i) - R33 * cosd(Theta3(i)) * Theta3prime(i);
    
    XgDoubleprime(i) = -R4 * cosd(Theta4(i)) * Theta4prime(i)^2 - R4 * sind(Theta4(i)) * Theta4Doubleprime(i) ...
        + R33 * cosd(Theta3(i)) * Theta3prime(i)^2 + R33 * sind(Theta3(i)) * Theta3Doubleprime(i);
    YgDoubleprime(i) = -R4 * sind(Theta4(i)) * Theta4prime(i)^2 + R4 * cosd(Theta4(i)) * Theta4Doubleprime(i) ...
        + R33 * sind(Theta3(i)) * Theta3prime(i)^2 - R33 * cosd(Theta3(i)) * Theta3Doubleprime(i);
    
    
    % PATH Analysis
    Rpprime(i) = sqrt(Xprime(i) ^ 2 + Yprime(i) ^ 2);
    
    Utx(i) = Xprime(i) / Rpprime(i);
    Uty(i) = Yprime(i) / Rpprime(i);
    Uny(i) = Xprime(i) / Rpprime(i);
    Unx(i) = -Yprime(i) / Rpprime(i);
    
    Rhop(i) = (Rpprime(i) ^ 3) / (Xprime(i) * Ydprime(i) - Yprime(i) * Xdprime(i));
    
    Xcc(i) = Xp(i) + Rhop(i) * (-Yprime(i) / Rpprime(i));
    Ycc(i) = Yp(i) + Rhop(i) * (Xprime(i) / Rpprime(i));
    
    Theta2 = Theta2 + 1;
end

Apx = Xdprime * Omega2a^2 + Xprime * Alpha2;
Apy = Ydprime * Omega2a^2 + Yprime * Alpha2;

Agx = XgDoubleprime * Omega2a^2 + Xgprime * Alpha2;
Agy = YgDoubleprime * Omega2a^2 + Ygprime * Alpha2;

%% Dynamic Force Analysis
Theta2 = 0;

for i = 1:360
    A = [-R6 * sind(Theta6(i)), -R6 * cosd(Theta6(i)); R5 * sind(Theta5(i)), R5 * cosd(Theta6(i))];
    B = [Ig6 * Alpha6a(i); Ig5 * Alpha5a(i) + m5 * (-R55 * cosd(Theta5(i)) * Apy(i) + R55 * sind(Theta5(i)) * Apx(i)) + g * R55 * cosd(Theta5(i))];
    X = linsolve(A,B);
    
    F65x(i) = X(1);
    F65y(i) = X(2);
    
    F45x(i) = m5 * Apx(i) - F65x(i);
    F45y(i) = m5 * Agy(i) + g * m5 - F65y(i);
    
    A = [R4 * sind(Theta4(i)), R4 * cosd(Theta4(i)); -R3 * sind(Theta3(i)), -R3 * cosd(Theta3(i))];
    B = [Ig4 * Alpha4a(i) + F45y(i) * R44 * cosd(Theta4(i)) + F45x(i) * R44 * sind(Theta4(i)); ...
        Ig3 * Alpha3a(i) + m3 * (-R33 * cosd(Theta3(i)) * Agy(i) + R33 * sind(Theta3(i)) * Agx(i)) + g * m3 * R33 * cosd(Theta3(i))];
    X = linsolve(A,B);
    
    F34x(i) = X(1);
    F34y(i) = X(2);
    
    F14x(i) = m4 * A4x + F45x(i) - F34x(i);
    F14y(i) = m4 * A4y + g * m4 + F45y(i) - F34y(i);
    
    F16x(i) = m6 * A6x + F65x(i);
    F16y(i) = m6 * A6y + g * m6 + F65y(i);
    
    F23x(i) = m3 * Agx(i) + F34x(i);
    F23y(i) = m3 * Agy(i) + g * m3 + F34y(i);
    
    F12x(i) = m2 * A2x + F23x(i);
    F12y(i) = m2 * A2y + g * m2 + F23y(i);
    
    T12(i) = Ig2 * Alpha2 + F23x(i) * R2 * sind(Theta2) + F23y(i) * R2 * cosd(Theta2);
    
    Theta2 = Theta2 + 1;
end

%% Static Force Analysis
% Theta2 = 0;
% 
% for i = 1:360
%     A = [-R6 * sind(Theta6(i)), -R6 * cosd(Theta6(i)); R5 * sind(Theta5(i)), R5 * cosd(Theta6(i))];
%     B = [0; g * R55 * cosd(Theta5(i))];
%     X = linsolve(A,B);
%     
%     F65x(i) = X(1);
%     F65y(i) = X(2);
%     
%     F45x(i) = -F65x(i);
%     F45y(i) = g * m5 - F65y(i);
%     
%     A = [R4 * sind(Theta4(i)), R4 * cosd(Theta4(i)); -R3 * sind(Theta3(i)), -R3 * cosd(Theta3(i))];
%     B = [F45y(i) * R44 * cosd(Theta4(i)) + F45x(i) * R44 * sind(Theta4(i)); ...
%         g * m3 * R33 * cosd(Theta3(i))];
%     X = linsolve(A,B);
%     
%     F34x(i) = X(1);
%     F34y(i) = X(2);
%     
%     F14x(i) = F45x(i) - F34x(i);
%     F14y(i) = g * m4 + F45y(i) - F34y(i);
%     
%     F16x(i) = F65x(i);
%     F16y(i) = g * m6 + F65y(i);
%     
%     F23x(i) = F34x(i);
%     F23y(i) = g * m3 + F34y(i);
%     
%     F12x(i) = F23x(i);
%     F12y(i) = g * m2 + F23y(i);
%     
%     T12(i) = F23x(i) * R2 * sind(Theta2) + F23y(i) * R2 * cosd(Theta2);
%     
%     Theta2 = Theta2 + 1;
% end

%% Plots
Theta2 = 1:360;

figure(1)
hold on
plot(Theta2,F65x)
plot(Theta2,F65y)
plot(Theta2,F34x)
plot(Theta2,F34y)
plot(Theta2,F16x)
plot(Theta2,F16y)
plot(Theta2,F12x)
plot(Theta2,F12y)
plot(Theta2,F45x)
plot(Theta2,F45y)
plot(Theta2,F14x)
plot(Theta2,F14y)
plot(Theta2,F23x)
plot(Theta2,F23y)
legend('F65x', 'F65y', 'F34x', 'F34y', 'F16x', 'F16y', 'F12x', 'F12y', 'F45x', 'F45y', 'F14x', 'F14y', 'F23x', 'F23y')
title('Reaction Forces Vs. Input Posture')
ylabel('Force (N)')
xlabel('Input Posture (deg)')

figure(2)
plot(Theta2,T12)
legend('T12')
title('Dynamic Torque @ w2 = 50 rad/s Vs. Input Posture')
ylabel('Torque (Nm)')
xlabel('Input Posture (deg)')

figure(3)
hold on
plot(Theta2,Apx)
plot(Theta2, Apy)
legend('Apx', 'Apy')
title('Acceleration of Point P Vs. Input Posture')
xlabel('Input Posture (deg)')
ylabel('Acceleration (m/s^2)')

%% Export Data
writematrix([transpose(Theta3) transpose(Theta4) transpose(Theta5) transpose(Theta6)],'Postures.txt')
writematrix([transpose(Theta3prime) transpose(Theta4prime) transpose(Theta5prime) transpose(Theta6prime)],'FirstOrderCoefficients.txt')
writematrix([transpose(Theta3Doubleprime) transpose(Theta4Doubleprime) transpose(Theta5Doubleprime) transpose(Theta6Doubleprime)],'SecondOrderCoefficients.txt')

writematrix([Omega3a', Omega4a', Omega5a', Omega6a',],'LinkVelocityA')
writematrix([Omega3b', Omega4b', Omega5b', Omega6b',],'LinkVelocityB')
writematrix([Alpha3a', Alpha4a', Alpha5a', Alpha6a',],'LinkAccelerationA')
writematrix([Alpha3b', Alpha4b', Alpha5b', Alpha6b',],'LinkAccelerationB')

writematrix([transpose(Xp) transpose(Yp)],'Position.txt')
writematrix([transpose(Xprime) transpose(Yprime)],'PointPFirstOrderCoefficients.txt')
writematrix([transpose(Xdprime) transpose(Ydprime)],'PointPSecondOrderCoefficients.txt')

writematrix([transpose(Utx) transpose(Uty) transpose(Unx) transpose(Uny) transpose(Rhop) transpose(Xcc) transpose(Ycc)],'PathAnalysis.txt')

writematrix([Xg' Yg' Xgprime' Ygprime' XgDoubleprime' YgDoubleprime' Agx' Agy'],'Link3Analysis');
writematrix([F65x' F65y' F34x' F34y' F16x' F16y' F12x' F12y' F45x' F45y' F14x' F14y' F23x' F23y'],'Dynamic Forces')
writematrix([T12'], 'Dynamic Torque')