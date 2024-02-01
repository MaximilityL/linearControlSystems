clear all
close all
clc

set(0,'DefaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultLineLineWidth', 1.2);

%% Q1 - LTI Controller

clear all
close all
clc

syms s;

plantNum = sym2poly(2.25 * (s + 1) * (s - 2))';
plantDen = sym2poly(s * (s^2 - 9))';


wantedPoles = [-10 , -20, -5 + i, -5 - i, -2 + 3*i , -2 - 3*i]; % 5 Poles that we want

deltaS = sym2poly((s - wantedPoles(1)) * (s - wantedPoles(2)) * (s - wantedPoles(3)) * (s - wantedPoles(4)) * (s - wantedPoles(5)) * (s - wantedPoles(6)))';

Ms = getSylvesterMatrix(plantNum, plantDen);
Ms1 = getMs1FromMs(Ms);
Ms2 = getMs2FromMs(Ms);
Ms3 = getMs3FromMs(Ms);

theta = Ms3 \ deltaS;

plantNum = plantNum';
plantDen = plantDen';
theta = theta';

controllerNum = theta(floor((length(theta)) / 2) + 1:length(theta));
controllerDen = theta(1:floor(length(theta) / 2));

P = tf(plantNum, plantDen);
C = tf(controllerNum, [controllerDen, 0]);



T = (P * C) / (1 + (P * C));
S = 1 / (1 + (P * C));
Tc = C / (1 + (P * C));
Td = P / (1 + (P * C));


t = 0:0.001:10;

Kr = 1;
Kd = 0.2;
startTime = 2;
delayTime = 0.1;
r = zeros(length(t),1);
d = r;
smallSinus = 0.1*sin(t);
n = awgn(smallSinus,10,'measured');

r(t > startTime) = Kr;
d(t > (startTime + delayTime)) = - Kd;

yr = lsim(T,r,t);
yd = lsim(Td,d,t);
yn = lsim(T,n,t);

y = yr + yd - yn;

figure
plot(t,yn);

L = P*C;
sys = feedback(L,1);
figure
step (T);

yr = lsim(sys,r,t);

figure
plot(t,yd);
%% Problem 5 - SVD

clear all
close all
clc

display('Problem 5 - SVD')
% Provided plant
P = [12, 19;
     8, 21];

% perform a singular value decomposition
[U,S,Vprime] = svd(P);

disp('U is: ');
disp(U);
disp('S is: ');
disp(S);
disp('Vprime is: ');
disp(Vprime);

disp('---------------');

Uprime = inv(U);
Sinverse = inv(S);
V = inv(Vprime);

disp('V is: ');
disp(V);
disp('Sinverse is: ');
disp(Sinverse);
disp('Uprime is: ');
disp(Uprime);

disp('---------------');


VTilda = reverseOrderOfVectors(V);
SinverseTilda = ArrangeDiagonalDescending(Sinverse);
UprimeTilda = inv(reverseOrderOfVectors(U));

disp('VTilda is: ');
disp(VTilda);
disp('SinverseTilda is: ');
disp(SinverseTilda);
disp('UprimeTilda is: ');
disp(UprimeTilda);

disp('---------------');
disp('');

umax = calcOutPutVectorOfSVDForIndex(VTilda, SinverseTilda, 1);
disp('umax is: ');
disp(umax);
disp(['The magnitude is: ', num2str(norm(umax)), ', corresponding singular value is: ', num2str(SinverseTilda(1, 1))]);
disp('');
umin = calcOutPutVectorOfSVDForIndex(VTilda, SinverseTilda, length(VTilda));
disp('');
disp('umin is: ');
disp(umin);
disp(['The magnitude is: ', num2str(norm(umin)), ', corresponding singular value is: ', num2str(SinverseTilda(length(VTilda), length(VTilda)))]);
disp('');


ymax = calcOutPutVectorOfSVDForIndex(U, S, 1);
disp('');
disp('ymax is: ');
disp(ymax);
disp(['The magnitude is: ', num2str(norm(ymax)), ', corresponding singular value is: ', num2str(S(1, 1))]);
disp('');


ymin = calcOutPutVectorOfSVDForIndex(U, S, length(U));
disp('');
disp('ymin is: ');
disp(ymin);
disp(['The magnitude is: ', num2str(norm(ymin)), ', corresponding singular value is: ', num2str(S(length(U), length(U)))]);
disp('');









%% Helper Functions

function DescendingDiagonalMatrix = ArrangeDiagonalDescending(inputMatrix)
    % Check if the input is a square matrix
    [rows, cols] = size(inputMatrix);
    if rows ~= cols
        error('Input matrix must be a square matrix.');
    end

    % Extract the diagonal elements
    diagonalElements = diag(inputMatrix);

    % Sort the diagonal elements in descending order
    sortedDiagonal = sort(diagonalElements, 'descend');

    % Create a new diagonal matrix with sorted elements
    DescendingDiagonalMatrix = diag(sortedDiagonal);
end

function reversedMatrix = reverseOrderOfVectors(inputMatrix)
    % Reverse the order of vectors
    reversedMatrix = inputMatrix(:, end:-1:1);
end

function outPutVector = calcOutPutVectorOfSVDForIndex(U,S,i)
    u = U(:,i);
    s = S(i,i);

    outPutVector = u * s;
end

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

function Ms1 = getMs1FromMs(Ms)
    Ms1 = Ms;
    numOfCol = length(Ms(1,:));
    n = (numOfCol / 2) + 1;
    Ms1(:,n) = [];
end

function Ms3 = getMs3FromMs(Ms)
    Ms3 = Ms;
    numOfCol = length(Ms(1,:));
    n = (numOfCol / 2) + 1;
    Ms3(:,n - 1) = [];
end

% function deltaS = getCharectaristicVector