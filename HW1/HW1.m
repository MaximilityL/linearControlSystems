clear all;
close all;
clc

%% Q1 - LTI Controller

syms s;

plantNum = sym2poly(2.25 * (s + 1) * (s - 2))';
plantDen = sym2poly(s * (s^2 - 9))';


plantNum = sym2poly(1*s^3 + 2*s^2 + 3*s +4)';
plantDen = plantNum;

wantedPoles = [0, -5, -10, 2 + 2*i, 2 - 2*i]; % 5 Poles that we want

deltaS = sym2poly((s - wantedPoles(1)) * (s - wantedPoles(2)) * (s - wantedPoles(3)) * (s - wantedPoles(4)) * (s - wantedPoles(5)));

M = transpose(getSylvesterMatrix(plantNum, plantDen));













%% Helper Functions

function sylvesterMatrix = getSylvesterMatrix(num, den)
    
    augmentedDen = [zeros(length(num) - length(den)) ; den];
    n = length(den) - 1;
    sylvesterMatrix = zeros(2 * n + 1, 2 * n + 2);
    Mden = zeros(2 * n + 1, n + 1);
    Mnum = Mden;

    denTempVector = zeros(1,length(augmentedDen))';
    numTempVector = zeros(1,length(num))';

    for i=1:(2 * n + 1)
        if (i > length(den))
            denValue = 0;
            numValue =0;
        else
            denValue = augmentedDen(i);
            numValue = num(i);
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



