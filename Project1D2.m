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
R1 = 203.20; %mm
R2 = 57.15; %mm
R3 = 184.15; %mm
R4 = 177.8; %mm
R5 = 50.80; %mm
R6 = 127.00; %mm
R11 = 101.60; %mm
R44 = 127.00; %mm
R55 = 25.40; %mm

resolution = 0.01; %deg

Theta3 = [];
Theta4 = [];
Theta5 = [];
Theta6 = [];

DeltaTheta3star = 1;
DeltaTheta4star = 1;
DeltaTheta5star = 1;
DeltaTheta6star = 1;

%%
for i = 1:360
    %% Section 1 Newton-Raphson
    iteration = 0;
    while ((abs(DeltaTheta3star) > resolution) | (abs(DeltaTheta4star) > resolution))
        
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
    
    %% Section 2 Newton-Raphson
    iteration = 0;
    Theta44 = Theta4star; %deg
    
    while (abs(DeltaTheta5star) > resolution | abs(DeltaTheta6star) > resolution)
        
        iteration = iteration + 1;
        
        epsilon_x2 = R6 .* cosd(Theta6star) - R5 .* cosd(Theta5star) - R44 .* cosd(Theta44) - R11 .* cosd(Theta11);
        epsilon_y2 = R6 .* sind(Theta6star) - R5 .* sind(Theta5star) - R44 .* sind(Theta44) - R11 .* sind(Theta11);
        
        DeltaTheta5star = ((((-epsilon_x2) .* R6 .* cosd(Theta6star)) + (-epsilon_y2 .* R6 .* sind(Theta6star))) ./ (R5 .* R6 .* sind(Theta5star - Theta6star))) .* (180 / pi); %deg
        DeltaTheta6star = ((((-epsilon_y2) .* R5 .* sind(Theta5star)) + (-epsilon_x2 .* R5 .* cosd(Theta6star))) ./ (R5 .* R6 .* sind(Theta5star - Theta6star))) .* (180 / pi);
        
        Theta5star = Theta5star + DeltaTheta5star; %deg
        Theta6star = Theta6star + DeltaTheta6star; %deg
        
%         if (abs(DeltaTheta5star) < resolution & abs(DeltaTheta6star) < resolution)
%             fprintf('Theta2 = %d\tDT5 = %.4f\tDT6 = %.4f\tEx = %.4f\tEy = %.4f\n', Theta2, DeltaTheta5star, DeltaTheta6star, epsilon_x2, epsilon_y2)
%         end
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
    disp(Theta5star)
    Det_2(i) = R5 .* R6 .* sind(Theta5star - Theta6star); %deg
    
    Theta5(i) = Theta5star;
    Theta6(i) = Theta6star;
    
    DeltaTheta5star = 1;
    DeltaTheta6star = 1;
    
%     fprintf('Theta2 = %.2f\tTheta3 = %.2f\tTheta4 = %.2f\tTheta5 = %.2f\tTheta6 = %.2f\tIterations = %d\n',Theta2, Theta3(i), Theta4(i), Theta5(i), Theta6(i), iteration)
    
    %% Kinematic Coefficients of Links
    Theta3prime(i) = (R2 .* sind(Theta2 - Theta4(i))) ./ (R3 .* sind(Theta4(i) - Theta3(i))); %rad
    Theta4prime(i) = (R2 .* sind(Theta2 - Theta3(i))) ./ (R4 .* sind(Theta4(i) - Theta3(i))); %rad
    
    A = [R5 * sind(Theta5(i)), -R6 * sind(Theta6(i)); -R5 * cosd(Theta5(i)), R6 * cosd(Theta6(i))];
    B = [-R44 * sind(Theta4(i)) * Theta4prime(i); R44 * cosd(Theta4(i)) * Theta4prime(i)];
    X = linsolve(A,B);
    Theta5prime(i) = X(1);
    Theta6prime(i) = X(2);
%     Theta5prime(i) = (R44 .* Theta4prime(i) .* sind(Theta4(i) - Theta6(i))) ./ (R5 .* sind(-Theta5(i) + Theta6(i))); %rad
%     Theta6prime(i) = (R44 .* Theta4prime(i) .* sind(Theta5(i) - Theta6(i))) ./ (R6 .* sind(-Theta5(i) + Theta6(i))); %rad
    
    SecondOrdDet1 = R3 * R4 * sind(Theta3(i) - Theta4(i));
    SecondOrdDet2 = R5 * R6 * sind(Theta5(i) - Theta6(i));
    
    Theta3Doubleprime(i) = (-R2 * R4 * cosd(Theta4(i) - Theta2) - R3 * R4 * cosd(Theta4(i) - Theta3(i)) * Theta3prime(i) ^ 2 + R4 ^ 2 * Theta4prime(i) ^ 2) / SecondOrdDet1;
    Theta4Doubleprime(i) = (-R2 * R3 * cosd(Theta3(i) - Theta2) - R3 * R4 * cosd(Theta3(i) - Theta4(i)) * Theta4prime(i) ^ 2 - R3 ^ 2 * Theta3prime(i) ^ 2) / SecondOrdDet1;
    
    Theta5Doubleprime(i) = (-R44 * R6 * sind(Theta4(i) - Theta6(i)) * Theta4Doubleprime(i) - R44 * R6 * cosd(Theta4(i) - Theta6(i)) * Theta4prime(i) ^ 2 ...
        - R5 * R6 * cosd(Theta5(i) - Theta6(i)) * Theta5prime(i) ^ 2 + R6 ^ 2 * Theta6prime(i) ^ 2) / SecondOrdDet2;
    
    Theta6Doubleprime(i) = (-R44 * R5 * sind(Theta4(i) - Theta5(i)) * Theta4Doubleprime(i) - R44 * R5 * cosd(Theta4(i) - Theta5(i)) * Theta4prime(i) ^ 2 ...
        + R5 * R6 * cosd(Theta6(i) - Theta5(i)) * Theta6prime(i) ^ 2 - R5 ^ 2 * Theta5prime(i) ^ 2) / SecondOrdDet2;
    
    %% Position of Point P
    Xp(i) = R1 * cosd(Theta1) + R44 * cosd(Theta4(i)) + R55 * cosd(Theta5(i));
    Yp(i) = R1 * sind(Theta1) + R44 * sind(Theta4(i)) + R55 * sind(Theta5(i));
    
    %% Kinematic Coefficients of Point P
    Xprime(i) = -R44 * sind(Theta4(i)) * Theta4prime(i) - R55 * sind(Theta5(i)) * Theta5prime(i);
    Yprime(i) = R44 * cosd(Theta4(i)) * Theta4prime(i) + R55 * cosd(Theta5(i)) * Theta5prime(i);
    
    Xdprime(i) = -R44 * sind(Theta4(i)) * Theta4Doubleprime(i) - R44 * cosd(Theta4(i)) * Theta4prime(i) ^ 2 - R55 * sind(Theta5(i)) * Theta5Doubleprime(i) ...
        + R55 * cosd(Theta5(i)) * Theta5prime(i) ^ 2;
    Ydprime(i) = R44 * cosd(Theta4(i)) * Theta4Doubleprime(i) - R44 * sind(Theta4(i)) * Theta4prime(i) ^ 2 + R55 * cosd(Theta5(i)) * Theta5Doubleprime(i) ...
        - R55 * sind(Theta5(i)) * Theta5prime(i) ^ 2;
    
    %% Path Analysis
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

writematrix([transpose(Theta3) transpose(Theta4) transpose(Theta5) transpose(Theta6)],'Postures.txt')
writematrix([transpose(Theta3prime) transpose(Theta4prime) transpose(Theta5prime) transpose(Theta6prime)],'FirstOrderCoefficients.txt')
writematrix([transpose(Theta3Doubleprime) transpose(Theta4Doubleprime) transpose(Theta5Doubleprime) transpose(Theta6Doubleprime)],'SecondOrderCoefficients.txt')

writematrix([transpose(Xp) transpose(Yp)],'Position.txt')
writematrix([transpose(Xprime) transpose(Yprime)],'PointPFirstOrderCoefficients.txt')
writematrix([transpose(Xdprime) transpose(Ydprime)],'PointPSecondOrderCoefficients.txt')

writematrix([transpose(Utx) transpose(Uty) transpose(Unx) transpose(Uny) transpose(Rhop) transpose(Xcc) transpose(Ycc)],'PathAnalysis.txt')
