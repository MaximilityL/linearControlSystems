clear all;
close all;
clc

%% Q1 - LTI Controller

syms s;

plantNum = sym2poly(2.25 * (s + 1) * (s - 2))';
plantDen = sym2poly(s * (s^2 - 9))';


wantedPoles = [0, -5, -10, 2 + 2*i, 2 - 2*i]; % 5 Poles that we want

deltaS = sym2poly((s - wantedPoles(1)) * (s - wantedPoles(2)) * (s - wantedPoles(3)) * (s - wantedPoles(4)) * (s - wantedPoles(5)))';

Ms = getSylvesterMatrix(plantNum, plantDen);
Ms2 = getMs2FromMs(Ms);

theta = Ms2 \ deltaS;

plantNum = plantNum';
plantDen = plantDen';
theta = theta';

controllerNum = theta((length(theta) / 2) + 1:length(theta));
controllerDen = theta(1:length(theta) / 2);

P = tf(plantNum, plantDen);
C = tf(controllerNum, controllerDen);

step(P*C);











%% Helper Functions

function sylvesterMatrix = getSylvesterMatrix(num, den)
    
    augmentedNum = [zeros(length(den) - length(num)) ; num];
    n = length(den) - 1;
    % sylvesterMatrix = zeros(2 * n + 1, 2 * n + 2);
    Mden = zeros(2 * n + 1, n + 1);
    Mnum = Mden;

    denTempVector = zeros(1,length(den))';
    numTempVector = zeros(1,length(augmentedNum))';

    for i=1:(2 * n + 1)
        if (i > length(den))
            denValue = 0;
            numValue =0;
        else
            denValue = den(i);
            numValue = augmentedNum(i);
        end

        denTempVector = [denValue; denTempVector];
        denTempVector(end) = [];

        numTempVector = [numValue; numTempVector];
        numTempVector(end) = [];

        Mden(i,:) = denTempVector;
        Mnum(i,:) = numTempVector;
    end
    
    sylvesterMatrix = [Mden, Mnum];
end

function Ms2 = getMs2FromMs(Ms)
    numOfCol = length(Ms(1,:));
    n = (numOfCol / 2) + 1;
    Ms2 = Ms;
    Ms2(1,:) = [];
    Ms2(:,n) = [];
    Ms2(:,1) = [];
end




