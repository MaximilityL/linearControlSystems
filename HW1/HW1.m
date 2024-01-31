clear all;
close all;
clc

%% Q1 - LTI Controller

syms s;

plantNum = sym2poly(2.25 * (s + 1) * (s - 2))';
plantDen = sym2poly(s * (s^2 - 9))';

wantedPoles = [0, -5, -10, 2 + 2*i, 2 - 2*i]; % 5 Poles that we want

deltaS = sym2poly((s - wantedPoles(1)) * (s - wantedPoles(2)) * (s - wantedPoles(3)) * (s - wantedPoles(4)) * (s - wantedPoles(5)));

M = transpose(getSylvesterMatrix(plantDen, plantNum));













%% Helper Functions

function sylvesterMatrix = getSylvesterMatrix(num, den)
    
    augmentedDen = [zeros(length(num) - length(den)) ; den];
    
    M = zeros(2 * length(den), 2 * length(den));
    Mden = zeros(2 * length(den), length(den));
    Mnum = Mden;
    
    for j=1:(length(den))
        Mden(1:length(den),j) = [sen)-(j-1))];


      sylvesterMatrix = 0;
    end
end




