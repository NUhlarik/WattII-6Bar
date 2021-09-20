close all
clear all
clc
%%
%Constants & Initial Guess

%Angles
Theta1 = 0; %deg
Theta2 = [0:10:350]; %deg
Theta3star = 30; %deg
Theta4star = 120; %deg
Theta5star = 0; %deg
Theta6star = 45; %deg
Theta11 = 0; %deg

%Radii
R1 = 203.20; %mm
R2 = 57.15; %mm
R3 = 184.15; %mm
R4 = 177.80; %mm
R5 = 50.80; %mm
R6 = 127.00; %mm
R11 = 101.60; %mm
R44 = 127.00; %mm

resolution = 0.01; %deg
iteration = 0;

DeltaTheta3star = 1;
DeltaTheta4star = 1;
DeltaTheta5star = 1;
DeltaTheta6star = 1;

%%
% Section 1 Newton-Raphson
while ((abs(DeltaTheta3star) > resolution) | (abs(DeltaTheta4star) > resolution))
    
    iteration = iteration + 1;
    
    epsilon_x = R2 .* cosd(Theta2) + R3 .* cosd(Theta3star) - R4 .* cosd(Theta4star) - R1 .* cosd(Theta1);
    epsilon_y = R2 .* sind(Theta2) + R3 .* sind(Theta3star) - R4 .* sind(Theta4star) - R1 .* sind(Theta1);
    
    DeltaTheta3star = (((epsilon_x .* R4 .* cosd(Theta4star)) + (epsilon_y .* R4 .* sind(Theta4star))) ./ (R3 .* R4 .* sind(Theta3star - Theta4star))) .* (180/pi); %deg
    DeltaTheta4star = (((epsilon_y .* R3 .* sind(Theta3star)) + (epsilon_x .* R3 .* cosd(Theta3star))) ./ (R3 .* R4 .* sind(Theta3star - Theta4star))) .* (180/pi); %deg
    
    Theta3star = Theta3star + DeltaTheta3star; %deg
    Theta4star = Theta4star + DeltaTheta4star; %deg
    
    if (iteration > 10)
        disp('Error')
        break
    end
end

Det_1 = R3 .* R4 .* sind(Theta3star - Theta4star) .* (180 / pi); %deg

%fprintf('%.2f\n', Theta3star)
%fprintf('%.2f\n', Theta4star)
%fprintf('%.2f\n', Det_1)

% Section 2 Newton-Raphson
iteration = 0;
Theta44 = Theta4star; %deg

while ((abs(DeltaTheta5star) > resolution) | (abs(DeltaTheta6star) > resolution))
    
    iteration = iteration + 1;
    
    epsilon_x2 = R6 .* cosd(Theta6star) - R5 .* cosd(Theta5star) - R44 .* cosd(Theta44) - R11 .* cosd(Theta11);
    epsilon_y2 = R6 .* sind(Theta6star) - R5 .* sind(Theta5star) - R44 .* sind(Theta44) - R11 .* sind(Theta11);
    
    DeltaTheta5star = ((((-epsilon_x2) .* R6 .* cosd(Theta6star)) - (epsilon_y2 .* R6 .* sind(Theta6star))) ./ (R5 .* R6 .* sind(Theta5star - Theta6star))) .* (180 / pi); %deg
    DeltaTheta6star = ((((-epsilon_y2) .* R5 .* sind(Theta5star)) - (epsilon_x2 .* R5 .* cosd(Theta6star))) ./ (R5 .* R6 .* sind(Theta5star - Theta6star))) .* (180 / pi);
    
    Theta5star = Theta5star + DeltaTheta5star; %deg
    Theta6star = Theta6star + DeltaTheta6star; %deg
    
    if (iteration > 10)
        disp('Error')
        break
    end
end

Det_2 = R5 .* R6 .* sind(Theta5star - Theta6star); %deg

%fprintf('%.2f\n', Theta5star)
%fprintf('%.2f\n', Theta6star)
%fprintf('%.2f\n', Det_2)

%%
%Kinematic Coefficients

Theta3prime = (R2 .* sind(Theta2 - Theta4star)) ./ (R3 .* sind(Theta4star - Theta3star)); %rad
Theta4prime = (R2 .* sind(Theta2 - Theta3star)) ./ (R4 .* sind(Theta4star - Theta3star)); %rad

% fprintf('%.4f\n',Theta3prime)
% fprintf('%.4f\n',Theta4prime)

% Theta5prime = (R44 .* Theta4prime .* R6 .* sind(Theta4star - Theta6star)) ./ (R5 .* R6 .* sind(Theta5star - Theta6star)); %rad
% Theta6prime = (R44 .* Theta4prime .* R5 .* sind(Theta5star - Theta6star)) ./ (R5 .* R6 .* sind(Theta5star - Theta6star)); %rad

for c = 1:36
    a11 = (R5 * sind(Theta5star(c)));
    a12 = (-R5 * cosd(Theta5star(c)));
    a21 = (-R6 * sind(Theta6star(c)));
    a22 = (R6 * cosd(Theta6star(c)));
    
    b11 = (-R44 .* sind(Theta4star(c)) .* Theta4prime(c));
    b12 = (-R44 .* cosd(Theta4star(c)) .* Theta4prime(c));
    
    A = [a11 a21; a12 a22];
    B = [b11; b12];
    
    X = linsolve(A,B);
    
    Theta5prime(c) = X(1,:);
    Theta6prime(c) = X(2,:);

end
fprintf('%.4f\n',Theta5prime)
% fprintf('%.4f\n',Theta6prime)

% plot(Theta2, Theta5prime)