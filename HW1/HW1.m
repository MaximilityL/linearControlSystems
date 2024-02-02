%% Q1 - LTI Controller

clear all
close all
clc

syms s;

plantNum = sym2poly(2.25 * (s + 1) * (s - 2))';
plantDen = sym2poly(s * (s^2 - 9))';


wantedPoles = [-2 , -20, -5 + i, -5 - i, -2 + 3*i , -2 - 3*i]; % 6 Poles that we want

Csym = getSymbolicControllerFromPlantWithWantedPoles()

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

controllerDen = [controllerDen, 0];
C = tf(controllerNum, controllerDen);


plantNumSym = poly2sym(plantNum', s);
plantDenSym = poly2sym(plantDen', s);
PSym = plantNumSym / plantDenSym;

controllerNumSym = poly2sym(controllerNum, s);
controllerDenSym = poly2sym(controllerDen, s);

CSym = controllerNumSym / controllerDenSym;

TSym = (CSym * PSym) / (1 + (CSym * PSym));
SSym = 1 / (1 + (CSym * PSym));
TcSym = CSym / (1 + (CSym * PSym));
TdSym = PSym / (1 + (CSym * PSym));

T = sym2tf(TSym);
S = sym2tf(SSym);
Tc = sym2tf(TcSym);
Td = sym2tf(TdSym);

t = 0:0.001:10;

Kr = 10;
Kd = 5;
startTime = 2;
delayTime = 0.2;
r = zeros(length(t),1);
d = r;

mu = 0; % bias
sigma = 0.5; % std

n = normrnd(mu, sigma,[1, length(t)]);

r(t > startTime) = Kr;
d(t > (startTime + delayTime)) = - Kd;

figure
lsim(T,r,t);
title('y - reference response')
grid
yr = lsim(T,r,t);

figure
lsim(Td,d,t);
title('y -disturbance response')
grid
yd = lsim(Td,d,t);

figure
lsim(T,n,t);
title('y - noise response')
grid
yn = lsim(T,n,t);

y = yr + yd - yn;

figure
plot(t,y,t,r);
grid
legend('Output', 'Reference');
title("y = yr + yd - yn")
xlabel('Time [sec]');
ylabel('Amplitude')

printResponsePlot('Entire Output', t, r, 'Reference Signal' , y, 'OutPut Signal');



figure
lsim(Tc,r,t);
title('u - reference response')
grid
ur = lsim(Tc,r,t);

figure
lsim(T,d,t);
title('u - disturbance response')
grid
ud = lsim(T,d,t);

figure
lsim(Tc,n,t);
title('u - noise response')
grid
un = lsim(Tc,n,t);

u = ur - ud - un;

figure
plot(t,u);
grid
title("u = ur - ud - un")
xlabel('Time [sec]');
ylabel('Amplitude')
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

function G = sym2tf(g)
    %SYM2TF - Symbolic Transfer Function to Numerical Transfer Function
    %conversion
    %
    % Syntax:  G = sym2tf(g)
    %
    % Inputs:
    %    g - Symbolic Transfer Function representation
    %
    % Outputs:
    %    G - Numeric Transfer Function representation
    %
    % Example: 
    %    syms p
    %    g11=(p + 2)/(p^2 + 2*p + 1);
    %    g12=(p - 1)/(p^2 + 5*p + 6);
    %    g21=(p - 1)/(p^2 + 3*p + 2);
    %    g22=(p + 2)/(p + 1);
    %    g=[g11 g12; g21 g22];
    %    G=sym2tf(g)
    %
    % See also: tf2sym, ss2sym
    %
    % Author: Oskar Vivero Osornio
    % email: oskar.vivero@gmail.com
    % Created: February 2006; 
    % Last revision: 25-March-2006;
    % May be distributed freely for non-commercial use, 
    % but please leave the above info unchanged, for
    % credit and feedback purposes
    %------------- BEGIN CODE -------------
    [n,m]=size(g);
    for i=1:n
        for j=1:m
            [num,den]=numden(g(i,j));
            num_n=sym2poly(num);
            den_n=sym2poly(den);
            G(i,j)=tf(num_n,den_n);
        end
    end
    %------------- END OF CODE --------------
end

function printResponsePlot(titleName, time, input, inputName, output, outputName)

    % Create a figure and set it to full screen
    hFig = figure('Position', get(0, 'Screensize'));

    plot(time, input, 'Color', [0.5, 0.5, 0.5]); % Gray color for input
    hold on;
    plot(time, output, 'Color', [0, 0, 1]);        % Blue color for output
    hold off;

    grid;
    fontSize = 20;
    legend({inputName, outputName}, 'Interpreter', 'latex', 'FontSize', fontSize);
    title(titleName, 'Interpreter', 'latex', 'FontSize', fontSize);
    xlabel('Time [sec]', 'Interpreter', 'latex', 'FontSize', fontSize);
    ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', fontSize);

    % Set the DPI (dots per inch) for the exported image
    dpi = 2000;

    % Save the plot as a JPG with the title as the filename
    filename = strrep(titleName, ' ', '_'); % Replace spaces with underscores
    print(hFig, filename, '-djpeg', ['-r', num2str(dpi)]);

end


