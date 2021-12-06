%% ME 452 Fall 2021
% Project 2
% Nicholas Uhlarik
% Deliverable 2
% 11 November 2021
% Updated 4 December 2021

close all
clear all
clc

%% Constants

%Global
resolution = 0.01; %deg
g = 9.810; %m/s^2

%Spring
k = 5000; %N/m
Rs0 = 150 / 1000; %m

%Damper
C = 350; %Ns/m
h = 101.6 / 1000; %m

% Link 1
R1 = 203.20 / 1000; %m
R11 = 101.60 / 1000; %m
Theta1 = 0; %deg

% Link 2
R2 = 57.15 / 1000; %m
m2 = 2; %kg
Ig2 = 0.0067; %kgm^2
Theta2 = (1:1:360)'; %deg
Theta2Short = Theta2(1:10:end,:);
empty = zeros(size(Theta2));
Theta2prime = empty + 1;
Theta2Doubleprime = empty;
Omega2a = 25; %rad/s
Omega2b = 50; %rad/s
Alpha2 = 0; %rad/s^2
A2x = 0; %m/s^2
A2y = 0; %m/s^2
XG2prime = empty;
YG2prime = empty;
XG2Doubleprime = empty;
YG2Doubleprime = empty;

% Link 3
R3 = 184.15 / 1000; %m
R33 = R3 / 2; %m
m3 = 5.5; %kg
Ig3 = 0.0433; %kgm^2
Theta3 = empty; %deg
Theta3prime = empty;
Theta3Doubleprime = empty;
Xg = empty; %m
Yg = empty; %m
Xgprime = empty;
Ygprime = empty;
XgDoubleprime = empty;
YgDoubleprime = empty;
DeltaTheta3star = 1; %rad
Det_34 = empty;

% Link 4
R4 = 177.8 / 1000; %m
R44 = 127.00 / 1000; %m
m4 = 7.5; %kg
Ig4 = 0.2426; %kgm^2
Theta4 = empty;
Theta4prime = empty;
Theta4Doubleprime = empty;
A4x = 0; %m/s^2
A4y = 0; %m/s^2
DeltaTheta4star = 1; %rad
XG4prime = empty;
YG4prime = empty;
XG4Doubleprime = empty;
YG4Doubleprime = empty;

% Link 5
R5 = 50.80 / 1000; %m
R55 = 25.40 / 1000; %m
m5 = 1.5; %kg
Ig5 = 0.0009; %kgm^2
Theta5 = empty;
Theta5prime = empty;
Theta5Doubleprime = empty;
Xp = empty; %m
Yp = empty; %m
Xprime = empty;
Yprime = empty;
Xdprime = empty;
Ydprime = empty;
Rpprime = empty;
DeltaTheta5star = 1; %rad
Utx = empty;
Uty = empty;
Unx = empty;
Uny = empty;
Rhop = empty;
Xcc = empty; %m
Ycc = empty; %m
Det_56 = empty;

%Link 6
R6 = 127.00 / 1000; %m
m6 = 6; %kg
Ig6 = 0.0634; %kgm^2
Theta6 = empty;
Theta6prime = empty;
Theta6Doubleprime = empty;
A6x = 0; %m/s^2
A6y = 0; %m/s^2
DeltaTheta6star = 1; %rad
XG6prime = empty;
YG6prime = empty;
XG6Doubleprime = empty;
YG6Doubleprime = empty;

% Dynamic Forces
F65x = empty; %N
F65y = empty; %N
F34x = empty; %N
F34y = empty; %N
F16x = empty; %N
F16y = empty; %N
F12x = empty; %N
F12y = empty; %N
F45x = empty; %N
F45y = empty; %N
F14x = empty; %N
F14y = empty; %N
F23x = empty; %N
F23y = empty; %N
T12 = empty; %Nm

% Static Forces
RF65x = empty; %N
RF65y = empty; %N
RF34x = empty; %N
RF34y = empty; %N
RF16x = empty; %N
RF16y = empty; %N
RF12x = empty; %N
RF12y = empty; %N
RF45x = empty; %N
RF45y = empty; %N
RF14x = empty; %N
RF14y = empty; %N
RF23x = empty; %N
RF23y = empty; %N
RT12 = empty; %Nm

% Spring
R7prime = empty;
Theta7prime = empty;

% Damper
DeltaRc = empty;

%% Postures

% Initial Posture Guesses
Theta3star = 60; %rad
Theta4star = 110; %rad
Theta5star = 320; %rad
Theta6star = 40; %rad

for i = 1:length(Theta2)
    
    % Links 3 & 4 Using Newton Raphson
    iteration = 0;
    while ((abs(DeltaTheta3star) > resolution) || (abs(DeltaTheta4star) > resolution))
        
        iteration = iteration + 1;
        
        epsilon_x = R2 * cosd(Theta2(i)) + R3 * cosd(Theta3star) - R4 * cosd(Theta4star) - R1 * cosd(Theta1);
        epsilon_y = R2 * sind(Theta2(i)) + R3 * sind(Theta3star) - R4 * sind(Theta4star) - R1 * sind(Theta1);
        
        DeltaTheta3star = (((epsilon_x * R4 * cosd(Theta4star)) + (epsilon_y * R4 * sind(Theta4star))) / (R3 * R4 * sind(Theta3star - Theta4star))) * (180/pi); %deg
        DeltaTheta4star = (((epsilon_y * R3 * sind(Theta3star)) + (epsilon_x * R3 * cosd(Theta3star))) / (R3 * R4 * sind(Theta3star - Theta4star))) * (180/pi); %deg
        
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
    
    Det_34(i) = R3 * R4 * sind(Theta3(i) - Theta4(i)) * (180 / pi); %deg
    
    % Links 5 & 6 Using Newton Raphson
    iteration = 0;
    while (abs(DeltaTheta5star) > resolution || abs(DeltaTheta6star) > resolution)
        
        iteration = iteration + 1;
        
        epsilon_x2 = R6 * cosd(Theta6star) - R5 * cosd(Theta5star) - R44 * cosd(Theta4(i)) - R11 * cosd(Theta1);
        epsilon_y2 = R6 * sind(Theta6star) - R5 * sind(Theta5star) - R44 * sind(Theta4(i)) - R11 * sind(Theta1);
        
        DeltaTheta5star = ((((-epsilon_x2) * R6 * cosd(Theta6star)) + (-epsilon_y2 * R6 .* sind(Theta6star))) / (R5 * R6 * sind(Theta5star - Theta6star))) * (180 / pi); %deg
        DeltaTheta6star = ((((-epsilon_y2) * R5 * sind(Theta5star)) + (-epsilon_x2 * R5 .* cosd(Theta6star))) / (R5 * R6 * sind(Theta5star - Theta6star))) * (180 / pi);
        
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
    
    Det_56(i) = R5 * R6 * sind(Theta5(i) - Theta6(i)); %deg
    
    Theta5(i) = Theta5star;
    Theta6(i) = Theta6star;
    
    DeltaTheta5star = 1;
    DeltaTheta6star = 1;
end

swing4 = max(Theta4) - min(Theta4);

%% Kinematic Coefficients

for i = 1:length(Theta2)
    %Link 3 and 4
    A = [-R3 * sind(Theta3(i)), R4 * sind(Theta4(i)); R3 * cosd(Theta3(i)), -R4 * cosd(Theta4(i))];
    B = [R2 * sind(Theta2(i)); -R2 * cosd(Theta2(i))];
    
    X = linsolve(A,B);
    Theta3prime(i) = X(1);
    Theta4prime(i) = X(2);
    
    B = [R2 * cosd(Theta2(i)) + R3 * cosd(Theta3(i)) * Theta3prime(i)^2 - R4 * cosd(Theta4(i)) * Theta4prime(i)^2; R2 * sind(Theta2(i)) + R3 * sind(Theta3(i)) * Theta3prime(i)^2 - R4 * sind(Theta4(i)) * Theta4prime(i)^2];
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

%% Kinematics of Points in Mechanism

for i = 1:length(Theta2)
    % Kinematics of Point P
    Xp(i) = R1 * cosd(Theta1) + R44 * cosd(Theta4(i)) + R55 * cosd(Theta5(i));
    Yp(i) = R1 * sind(Theta1) + R44 * sind(Theta4(i)) + R55 * sind(Theta5(i));
    
    Xprime(i) = -R44 * sind(Theta4(i)) * Theta4prime(i) - R55 * sind(Theta5(i)) * Theta5prime(i);
    Yprime(i) = R44 * cosd(Theta4(i)) * Theta4prime(i) + R55 * cosd(Theta5(i)) * Theta5prime(i);
    
    Xdprime(i) = -R44 * (cosd(Theta4(i)) * Theta4prime(i)^2 + sind(Theta4(i)) * Theta4Doubleprime(i)) ...
        -R55 * (cosd(Theta5(i)) * Theta5prime(i)^2 + sind(Theta5(i)) * Theta5Doubleprime(i));
    Ydprime(i) = -R44 * (sind(Theta4(i)) * Theta4prime(i)^2 - cosd(Theta4(i)) * Theta4Doubleprime(i)) ...
        -R55 * (sind(Theta5(i)) * Theta5prime(i)^2 - cosd(Theta5(i)) * Theta5Doubleprime(i));
    
    % Kinematics of Point G
    Xg(i) = R1 * cosd(Theta1) + R4 * cosd(Theta4(i)) - R33 * cosd(Theta3(i));
    Yg(i) = R1 * sind(Theta1) + R4 * sind(Theta4(i)) - R33 * sind(Theta3(i));
    
    Xgprime(i) = -R4 * sind(Theta4(i)) * Theta4prime(i) + R33 * sind(Theta3(i)) * Theta3prime(i);
    Ygprime(i) = R4 * cosd(Theta4(i)) * Theta4prime(i) - R33 * cosd(Theta3(i)) * Theta3prime(i);
    
    XgDoubleprime(i) = -R4 * cosd(Theta4(i)) * Theta4prime(i)^2 - R4 * sind(Theta4(i)) * Theta4Doubleprime(i) ...
        + R33 * cosd(Theta3(i)) * Theta3prime(i)^2 + R33 * sind(Theta3(i)) * Theta3Doubleprime(i);
    YgDoubleprime(i) = -R4 * sind(Theta4(i)) * Theta4prime(i)^2 + R4 * cosd(Theta4(i)) * Theta4Doubleprime(i) ...
        + R33 * sind(Theta3(i)) * Theta3prime(i)^2 - R33 * cosd(Theta3(i)) * Theta3Doubleprime(i);
    
    
    % Point P Path Analysis
    Rpprime(i) = sqrt(Xprime(i) ^ 2 + Yprime(i) ^ 2);
    
    Utx(i) = Xprime(i) / Rpprime(i); %m/s
    Uty(i) = Yprime(i) / Rpprime(i); %m/s
    Uny(i) = Xprime(i) / Rpprime(i); %m/s
    Unx(i) = -Yprime(i) / Rpprime(i); %m/s
    
    Rhop(i) = (Rpprime(i) ^ 3) / (Xprime(i) * Ydprime(i) - Yprime(i) * Xdprime(i)); %m
    
    Xcc(i) = Xp(i) + Rhop(i) * (-Yprime(i) / Rpprime(i)); %m
    Ycc(i) = Yp(i) + Rhop(i) * (Xprime(i) / Rpprime(i)); %m
    
end

% Minimum and Maximum Displacements of Point P
[Xmin_P, XminPos] = min(Xp);
[Xmax_P, XmaxPos] = max(Xp);
[Ymin_P, YminPos] = min(Yp);
[Ymax_P, YmaxPos] = max(Yp);

Xmin_Theta2 = Theta2(XminPos);
Xmax_Theta2 = Theta2(XmaxPos);
Ymin_Theta2 = Theta2(YminPos);
Ymax_Theta2 = Theta2(YmaxPos);

deltaXp = Xmax_P - Xmin_P;
deltaYp = Ymax_P - Ymin_P;

% Velocity of Points
Vpx = Xprime * Omega2a;
Vpy = Yprime * Omega2a;
Vgx = Xgprime * Omega2a;
Vgy = Ygprime * Omega2a;

% Acceleration of Points
Apx = Xdprime * Omega2a^2 + Xprime * Alpha2;
Apy = Ydprime * Omega2a^2 + Yprime * Alpha2;

Agx = XgDoubleprime * Omega2a^2 + Xgprime * Alpha2;
Agy = YgDoubleprime * Omega2a^2 + Ygprime * Alpha2;

%% Dynamic Force Analysis

for i = 1:length(Theta2)
    A = [-R6 * sind(Theta6(i)), -R6 * cosd(Theta6(i)); R5 * sind(Theta5(i)), R5 * cosd(Theta6(i))];
    B = [Ig6 * Alpha6a(i); Ig5 * Alpha5a(i) + m5 * (-R55 * cosd(Theta5(i)) * Apy(i) + R55 * sind(Theta5(i)) * Apx(i)) + m5 * g * R55 * cosd(Theta5(i))];
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
    
    T12(i) = Ig2 * Alpha2 + F23x(i) * R2 * sind(Theta2(i)) + F23y(i) * R2 * cosd(Theta2(i));

end

%% Static Force Analysis

for i = 1:length(Theta2)
    A = [-R6 * sind(Theta6(i)), -R6 * cosd(Theta6(i)); R5 * sind(Theta5(i)), R5 * cosd(Theta6(i))];
    B = [0; m5 * g * R55 * cosd(Theta5(i))];
    X = linsolve(A,B);
    
    RF65x(i) = X(1);
    RF65y(i) = X(2);
    
    RF45x(i) = -RF65x(i);
    RF45y(i) = g * m5 - RF65y(i);
    
    A = [R4 * sind(Theta4(i)), R4 * cosd(Theta4(i)); -R3 * sind(Theta3(i)), -R3 * cosd(Theta3(i))];
    B = [RF45y(i) * R44 * cosd(Theta4(i)) + RF45x(i) * R44 * sind(Theta4(i)); ...
        g * m3 * R33 * cosd(Theta3(i))];
    X = linsolve(A,B);
    
    RF34x(i) = X(1);
    RF34y(i) = X(2);
    
    RF14x(i) = RF45x(i) - RF34x(i);
    RF14y(i) = g * m4 + RF45y(i) - RF34y(i);
    
    RF16x(i) = RF65x(i);
    RF16y(i) = g * m6 + RF65y(i);
    
    RF23x(i) = RF34x(i);
    RF23y(i) = g * m3 + RF34y(i);
    
    RF12x(i) = RF23x(i);
    RF12y(i) = g * m2 + RF23y(i);
    
    RT12(i) = RF23x(i) * R2 * sind(Theta2(i)) + RF23y(i) * R2 * cosd(Theta2(i));

end

%% Kinetic Energy & Equivalent Moments of Intertia
A2 = m2 * (XG2prime.^2 + YG2prime.^2) + Ig2 * Theta2prime.^2; %kg*m^2
A3 = m3 * (Xgprime.^2 + Ygprime.^2) + Ig3 * Theta3prime.^2; %kg*m^2
A4 = m4 * (XG4prime.^2 + YG4prime.^2) + Ig4 * Theta4prime.^2; %kg*m^2
A5 = m5 * (Xprime.^2 + Yprime.^2) + Ig5 * Theta5prime.^2; %kg*m^2
A6 = m6 * (XG6prime.^2 + YG6prime.^2) + Ig6 * Theta6prime.^2; %kg*m^2

Ieq = A2 + A3 + A4 + A5 + A6; %kg*m^2

KE = 0.5 * Ieq * Omega2a^2; %Joule

B2 = m2 * (XG2prime .* XG2Doubleprime + YG2prime .* YG2Doubleprime) + Ig2 * Theta2prime .* Theta2Doubleprime; %kg*m^2
B3 = m3 * (Xgprime .* XgDoubleprime + Ygprime .* YgDoubleprime) + Ig3 * Theta3prime .* Theta3Doubleprime; %kg*m^2
B4 = m4 * (XG4prime .* XG4Doubleprime + YG4prime .* YG4Doubleprime) + Ig4 * Theta4prime .* Theta4Doubleprime; %kg*m^2
B5 = m5 * (Xprime .* Xdprime + Yprime .* Ydprime) + Ig5 * Theta5prime .* Theta5Doubleprime; %kg*m^2
B6 = m6 * (XG6prime .* XG6Doubleprime + YG6prime .* YG6Doubleprime) + Ig6 * Theta6prime .* Theta6Doubleprime; %kg*m^2

Btot = B2 + B3 + B4 + B5 + B6; %kg*m^2

dEdt = Ieq * Omega2a * Alpha2 + Btot * Omega2a ^ 3; %Watt

%% Spring Potential Energy
R7 = sqrt(Xp.^2 + Yp.^2);
Theta7 = acosd(Xp ./ R7);

for i = 1:length(Theta2)
    A = [cosd(Theta7(i)), -R7(i)*sind(Theta7(i)); sind(Theta7(i)), R7(i)*cosd(Theta7(i))];
    B = [-R55*sind(Theta5(i))*Theta5prime(i)-R44*sind(Theta4(i))*Theta4prime(i); R55*cosd(Theta5(i))*Theta5prime(i)+R44*cosd(Theta4(i))*Theta4prime(i)];
    
    x = linsolve(A,B);
    
    R7prime(i) = x(1);
    Theta7prime(i) = x(2);
end

Us = 0.5 * k * (R7 - Rs0).^2; %Joules
dUsdt = k * (R7 - Rs0) .* R7prime * Omega2a;

%% Gravitational Potential Energy
Ug3 = m3 * g * Yg; %Joules
Ug5 = m5 * g * Yp; %Joules
Ugtot = Ug3 + Ug5; %Joules
dUgdt = (m3 * g * Ygprime + m5 * g * Yprime) * Omega2a; %Watts

%% Viscous Damper Dissipation Energy
R8prime = R44 * sind(Theta4) .* Theta4prime + R55 * sind(Theta5) .* Theta5prime;
for i = 1:length(Theta2)
    if (i == 1)
        DeltaRc(i) = Xp(i) - Xp(end);
    end
    if (i > 1)
        DeltaRc(i) = Xp(i) - Xp(i-1);
    end
end

Wc = C * R8prime * Omega2a .* DeltaRc; %Joules
dWcdt = C * R8prime.^2 * Omega2a^2; %Watts

%% Power Equation & Equation of Motion
Power = dEdt + dUsdt + dUgdt + dWcdt; %Watts
EOM_T12 = Power / Omega2a; %Nm

%% Postures Plots
% figure();
% subplot(2,1,1);
% hold on
% plot(Theta2, Theta3);
% plot(Theta2, Theta4);
% plot(Theta2, Theta5);
% plot(Theta2, Theta6);
% title('Link Postures Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Link Postures (deg)');
% legend('\theta_3', '\theta_4', '\theta_5', '\theta_6', 'Location', 'eastoutside');
% hold off
% Postures = horzcat(Theta3, Theta4, Theta5, Theta6);
% Postures = Postures(1:10:end,:);

%% Determinants Plots
% subplot(2,1,2);
% hold on
% plot(Theta2, Det_34);
% plot(Theta2, Det_56);
% title('Determinants Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Determinants (m^2)');
% legend('Links 3 & 4', 'Links 5 & 6', 'Location', 'eastoutside');
% hold off
% Determinants = horzcat(Det_34, Det_56);
% Determinants = Determinants(1:10:end, :);

%% Postures Plots
% figure();
% subplot(2,1,1);
% hold on
% plot(Theta2, Theta3prime);
% plot(Theta2, Theta4prime);
% plot(Theta2, Theta5prime);
% plot(Theta2, Theta6prime);
% title('Link First Order Kinemtaic Coefficients Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Kinematic Coefficients');
% legend('\theta_3''', '\theta_4''', '\theta_5''', '\theta_6''', 'Location', 'eastoutside');
% hold off
% FirstOrderKC = horzcat(Theta3prime, Theta4prime, Theta5prime, Theta6prime);
% FirstOrderKC = FirstOrderKC(1:10:end,:);
% 
% subplot(2,1,2);
% hold on
% plot(Theta2, Theta3Doubleprime);
% plot(Theta2, Theta4Doubleprime);
% plot(Theta2, Theta5Doubleprime);
% plot(Theta2, Theta6Doubleprime);
% title('Link Second Order Kinemtaic Coefficients Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Second Order Kinematic Coefficients');
% legend('\theta_3"', '\theta_4"', '\theta_5"', '\theta_6"', 'Location', 'eastoutside');
% hold off
% SecondOrderKC = horzcat(Theta3Doubleprime, Theta4Doubleprime, Theta5Doubleprime, Theta6Doubleprime);
% SecondOrderKC = SecondOrderKC(1:10:end,:);

%% Angular Velocities Plots
% figure();
% subplot(2,1,1);
% hold on
% plot(Theta2, Omega3a);
% plot(Theta2, Omega4a);
% plot(Theta2, Omega5a);
% plot(Theta2, Omega6a);
% title('Link Angular Velocities Vs. Input Angle for \omega_2 = 25 rad/s');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Angular Velocity (rad/s)');
% legend('\omega_3', '\omega_4', '\omega_5', '\omega_6', 'Location', 'eastoutside');
% hold off
% AngularVelocityA = horzcat(Omega3a, Omega4a, Omega5a, Omega6a);
% AngularVelocityA = AngularVelocityA(1:10:end,:);
% 
% subplot(2,1,2);
% hold on
% plot(Theta2, Omega3b);
% plot(Theta2, Omega4b);
% plot(Theta2, Omega5b);
% plot(Theta2, Omega6b);
% title('Link Angular Velocities Vs. Input Angle for \omega_2 = 50 rad/s');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Angular Velocity (rad/s)');
% legend('\omega_3', '\omega_4', '\omega_5', '\omega_6', 'Location', 'eastoutside');
% hold off
% AngularVelocityB = horzcat(Omega3b, Omega4b, Omega5b, Omega6b);
% AngularVelocityB = AngularVelocityB(1:10:end,:);

%% Angular Acceleration Plots
% figure();
% subplot(2,1,1);
% hold on
% plot(Theta2, Alpha3a);
% plot(Theta2, Alpha4a);
% plot(Theta2, Alpha5a);
% plot(Theta2, Alpha6a);
% title('Link Angular Accelerations Vs. Input Angle for \omega_2 = 25 rad/s');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Angular Acceleration (rad/s^2)');
% legend('\alpha_3', '\alpha_4', '\alpha_5', '\alpha_6', 'Location', 'eastoutside');
% hold off
% AngularAccelerationA = horzcat(Alpha3a, Alpha4a, Alpha5a, Alpha6a);
% AngularAccelerationA = AngularAccelerationA(1:10:end,:);
% 
% subplot(2,1,2);
% hold on
% plot(Theta2, Alpha3b);
% plot(Theta2, Alpha4b);
% plot(Theta2, Alpha5b);
% plot(Theta2, Alpha6b);
% title('Link Angular Accelerations Vs. Input Angle for \omega_2 = 50 rad/s');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Angular Acceleration (rad/s^2)');
% legend('\alpha_3', '\alpha_4', '\alpha_5', '\alpha_6', 'Location', 'eastoutside');
% hold off
% AngularAccelerationB = horzcat(Omega3b, Omega4b, Omega5b, Omega6b);
% AngularAccelerationB = AngularAccelerationB(1:10:end,:);

%% Kinematics of Point P Plots
% % Position of Point P
% figure();
% plot(Xp, Yp)
% title('Position of Point P Through Full Travel');
% xlabel('X Position (m)')
% ylabel('Y Position (m)')
% P_pos = horzcat(Xp, Yp);
% P_pos = P_pos(1:10:end, :);
% 
% % First Order Kinematic Coefficients of Point P
% figure();
% subplot(2,1,1);
% hold on
% plot(Theta2, Xprime)
% plot(Theta2, Yprime)
% title('First Order Kinematic Coefficients of Point P Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Kinematic Coefficients');
% legend('X''', 'Y''', 'Location', 'eastoutside');
% hold off
% PFirstOrderKC = horzcat(Xprime, Yprime);
% PFirstOrderKC = PFirstOrderKC(1:10:end, :);
% 
% % Second Order Kinematic Coefficients of Point P
% subplot(2,1,2);
% hold on
% plot(Theta2, Xdprime)
% plot(Theta2, Ydprime)
% title('Second Order Kinematic Coefficients of Point P Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Second Order Kinematic Coefficients');
% legend('X"', 'Y"', 'Location', 'eastoutside');
% hold off
% PSecondOrderKC = horzcat(Xdprime, Ydprime);
% PSecondOrderKC = PSecondOrderKC(1:10:end, :);
% 
% % Velocity and Acceleration of Point P
% figure();
% subplot(2,1,1);
% hold on
% plot(Theta2, Vpx)
% plot(Theta2, Vpy)
% title('Velocity of Point P Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 306]);
% ylabel('Velocity (m/s)');
% legend('V_x', 'V_y', 'Location', 'eastoutside');
% hold off
% PVelocity = horzcat(Vpx, Vpy);
% PVelocity = PVelocity(1:10:end, :);
% 
% % Second Order Kinematic Coefficients of Point P
% subplot(2,1,2);
% hold on
% plot(Theta2, Apx)
% plot(Theta2, Apy)
% title('Acceleration of Point P Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Acceleration (m/s^2)');
% legend('A_x', 'A_y', 'Location', 'eastoutside');
% hold off
% PAcceleration = horzcat(Apx, Apy);
% PAcceleration = PAcceleration(1:10:end, :);

%% Dynamic Forces Plots
% figure();
% hold on
% plot(Theta2, F65x);
% plot(Theta2, F65y);
% plot(Theta2, F34x);
% plot(Theta2, F34y);
% plot(Theta2, F16x);
% plot(Theta2, F16y);
% plot(Theta2, F12x);
% plot(Theta2, F12y);
% plot(Theta2, F45x);
% plot(Theta2, F45y);
% plot(Theta2, F14x);
% plot(Theta2, F14y);
% plot(Theta2, F23x);
% plot(Theta2, F23y);
% title('Dynamic Link Forces Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Force (N)');
% legend('F65x', 'F65y', 'F34x', 'F34y', 'F16x', 'F16y', 'F12x', 'F12y', 'F45x', 'F45y', 'F14x', 'F14y', 'F23x', 'F23y', 'Location', 'eastoutside')
% hold off
% DynamicForces = horzcat(F65x, F65y, F34x, F34y, F16x, F16y, F12x, F12y, F45x, F45y, F14x, F14y, F23x, F23y);
% DynamicForces = DynamicForces(1:10:end,:);
% 
figure();
plot(Theta2, T12)
title('Link 2 Dynamic Torque Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Torque (Nm)');
DynamicTorque = T12;
DynamicTorque = DynamicTorque(1:10:end,:); 

%% Static Forces Plots
% figure();
% hold on
% plot(Theta2, RF65x);
% plot(Theta2, RF65y);
% plot(Theta2, RF34x);
% plot(Theta2, RF34y);
% plot(Theta2, RF16x);
% plot(Theta2, RF16y);
% plot(Theta2, RF12x);
% plot(Theta2, RF12y);
% plot(Theta2, RF45x);
% plot(Theta2, RF45y);
% plot(Theta2, RF14x);
% plot(Theta2, RF14y);
% plot(Theta2, RF23x);
% plot(Theta2, RF23y);
% title('Static Reaction Forces Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Force (N)');
% legend('RF65x', 'RF65y', 'RF34x', 'RF34y', 'RF16x', 'RF16y', 'RF12x', 'RF12y', 'RF45x', 'RF45y', 'RF14x', 'RF14y', 'RF23x', 'RF23y', 'Location', 'eastoutside')
% hold off
% StaticForces = horzcat(RF65x, RF65y, RF34x, RF34y, RF16x, RF16y, RF12x, RF12y, RF45x, RF45y, RF14x, RF14y, RF23x, RF23y);
% StaticForces = StaticForces(1:10:end,:);
% 
% figure();
% plot(Theta2, RT12)
% title('Link 2 Static Torque Vs. Input Angle');
% xlabel('Input Angle \theta_2 (deg)');
% xlim([0 360]);
% ylabel('Torque (Nm)');
% StaticTorque = RT12;
% StaticTorque = StaticTorque(1:10:end,:);

%% Kinetic Energy Plots
figure();
subplot(3,1,1)
plot(Theta2, Ieq)
title('Mechanism Equivalent Mass Moment of Inertia Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('I_E_Q (kg*m^2)');

subplot(3,1,2)
plot(Theta2, KE)
title('Mechanism Kinetic Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Kinetic Energy (Joules)');

subplot(3,1,3)
plot(Theta2, dEdt)
title('Mechanism Time Rate of Change of Kinetic Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('dE/dt (Watts)');
KineticEnergy = horzcat(Ieq, KE, dEdt);
KineticEnergy = KineticEnergy(1:10:end,:);

%% Spring Energy Plots
figure();
plot(Theta2, R7prime)
title('Spring Kinematic Coefficient Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Kinematic Coefficient');

figure();
subplot(2,1,1)
plot(Theta2, Us)
title('Spring Potential Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Potential Energy (Joules)');

subplot(2,1,2)
plot(Theta2, dUsdt)
title('Spring Time Rate of Change of Potential Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('dU_s/dt (Watts)');
Spring = horzcat(R7prime, Us, dUsdt);
Spring = Spring(1:10:end,:);

%% Gravitational Energy Plots
figure();
hold on;
plot(Theta2, Ug3)
plot(Theta2, Ug5)
plot(Theta2, Ugtot)
title('Gravitational Potential Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Gravitational Potential Energy (Joules)');
legend('Ug Link 3', 'Ug Link 5', 'Ug Total', 'Location', 'eastoutside');
hold off;

figure();
plot(Theta2, dUgdt)
title('Time Rate of Change of Gravitational Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('dU_g/dt (Watts)');
Gravity = horzcat(Ug3, Ug5, Ugtot, dUgdt);
Gravity = Gravity(1:10:end,:);

%% Damper Dissipation Energy Plots
figure();
plot(Theta2, R8prime)
title('Damper Kinematic Coefficient Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Kinematic Coefficient');

figure();
subplot(2,1,1)
plot(Theta2, Wc)
title('Viscous Damper Dissipation Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Dissipation Energy (Joules)');

subplot(2,1,2)
plot(Theta2, dWcdt)
title('Time Rate of Change of Viscous Damper Dissipation Energy Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('dW_c/dt (Watts)');
Damper = horzcat(R8prime, Wc, dWcdt);
Damper = Damper(1:10:end,:);

%% Power Equation & EOM Plots
figure();
plot(Theta2, EOM_T12)
title('Equation of Motion Torque Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Torque (Nm)');
MotorTorque = EOM_T12;
MotorTorque = MotorTorque(1:10:end,:);

figure();
hold on
plot(Theta2, Power)
plot(Theta2, dEdt)
plot(Theta2, dUsdt)
plot(Theta2, dUgdt)
plot(Theta2, dWcdt)
title('Power Contribution Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Power (Watts)');
legend('Net Power', 'Kinetic Power', 'Spring Power', 'Gravitational Power', 'Damper Power', 'Location', 'eastoutside');
hold off
MotorPower = horzcat(Power, dEdt, dUsdt, dUgdt, dWcdt);
MotorPower = MotorPower(1:10:end,:);

figure();
hold on
plot(Theta2, EOM_T12)
plot(Theta2, T12)
title('EOM Torque & Newton Euler Vs. Input Angle');
xlabel('Input Angle \theta_2 (deg)');
xlim([0 360]);
ylabel('Torque (Nm)');
legend('EOM', 'Newton Euler', 'Location', 'eastoutside');
hold off

%% Spreadsheets
Header = ["Theta2 (Deg)", "Torque (Nm)"];
writematrix([Header; [Theta2Short, DynamicTorque]], 'DynamicTorque.xls');

Header = ["Theta2 (Deg)", "I_EQ (kg*m^2)", "KE (Joules)", "dE/dt (Watt)"];
writematrix([Header; [Theta2Short, KineticEnergy]], 'KineticEnergy.xls');

Header = ["Theta2 (Deg)", "Rs'", "Us (Joules)", "dUs/dt (Watt)"];
writematrix([Header; [Theta2Short, Spring]], 'Spring.xls');

Header = ["Theta2 (Deg)", "Ug Link 3 (Joules)", "Ug Link 5 (Joules)", "Ug Total (Joules)" "dUg/dt (Watt)"];
writematrix([Header; [Theta2Short, Gravity]], 'Gravity.xls');

Header = ["Theta2 (Deg)", "Rc'", "Wc (Joules)", "dWc/dt (Watt)"];
writematrix([Header; [Theta2Short, Damper]], 'Damper.xls');

Header = ["Theta2 (Deg)", "Motor Torque (Nm)"];
writematrix([Header; [Theta2Short, MotorTorque]], 'MotorTorque.xls');

Header = ["Theta2 (Deg)", "Net Power (Watt)", "dE/dt (Watt)", "dUs/dt (Watt)" "dUg/dt (Watt)", "dWc/dt (Watt)"];
writematrix([Header; [Theta2Short, MotorPower]], 'MotorPower.xls');

Header = ["Theta2 (Deg)", "EOM Motor Torque (Nm)", "Newton Euler Motor Torque (Nm)"];
writematrix([Header; [Theta2Short, MotorTorque, DynamicTorque]], 'ComparisonTorque.xls');
