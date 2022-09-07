%% init

isCallingFromOutside = 0; %must be 1 when executing an optimization
if ~isCallingFromOutside, clear all; isCallingFromOutside = 0; disp('IS NOT CALLING FROM OUTSIDE'); end

% settings when not calling from outside
if ~isCallingFromOutside
    isEvaluateComsol = 0;
end

isEnablePlotting = 1;

positionOfAerodynamicForce = [0.088; 0.018];
positionOfAerodynamicForceOnChordLine = [0.085; -0.0195];
distanceVerticalBeamsGap = 30e-3;

positionOfFlapAxis = [0.163; -0.035];

if 1 %use values from optimization
    isCallingFromOutside = 1; %same approach
    
    %30.09.2020
    x0 = [-0.96875,0.16625,0.08890625,0.1675,-0.14,0.61,0.1203125,0.15015625,-0.15,0.96,-0.0796875,-36.9375,0.045,0.27,5]';

    %x -> variables
    x = x0;
    
    i = 1;

    %... general

    angleOfAttackWrtMain = deg2rad(x(i)); i = i + 1;

    %... main element

    main_x1  = x(i); i = i + 1;
    main_ch1 = x(i); i = i + 1;
    main_th1 = x(i); i = i + 1;
    main_tw1 = x(i); i = i + 1;
    main_x2  = x(i); i = i + 1;
    main_ch2 = x(i); i = i + 1;
    main_th2 = x(i); i = i + 1;
    main_tw2 = x(i); i = i + 1;

    %... flap element

    flap_positionX = x(i); i = i + 1;
    flap_positionY = x(i); i = i + 1;
    flap_rotation = deg2rad(x(i)); i = i + 1;

    flap_maximumChamber = x(i); i = i + 1;
    flap_positionOfMaximumChamber = x(i); i = i + 1;

    flap_gurneyLengthPercentage = x(i); i = i + 1;
end



%% cost initialization

cost = 0;

penaltyOffset = 1e3;
penaltyMultiplyer = 1e3;



%% airfoil parameters

% general

if ~isCallingFromOutside
    angleOfAttackWrtMain = deg2rad(-5); %OPTIMIZE
end
chord = 0.15; %chord of main element



% main element

main_yTE = 0.5e-3 / chord; %assuming we can manufacture 0.5mm trailing edge

if ~isCallingFromOutside
    main_x1  = 0.18; %OPTIMIZE
    main_ch1 = 0.10; %OPTIMIZE
    main_th1 = 0.16; %OPTIMIZE
    main_tw1 = -0.15; %OPTIMIZE
    main_x2  = 0.60; %OPTIMIZE
    main_ch2 = 0.12; %OPTIMIZE
    main_th2 = 0.16; %OPTIMIZE
    main_tw2 = -0.15; %OPTIMIZE
end



% flap element

flap_size = 1/2; %keep fix
if ~isCallingFromOutside
    flap_positionX = 0.95; %OPTIMIZE
    flap_positionY = -0.08; %OPTIMIZE
    flap_rotation = deg2rad(-35); %OPTIMIZE
end

flap_thickness = 0.15; %keep fix
if ~isCallingFromOutside
    flap_maximumChamber = 0.05; %OPTIMIZE
    flap_positionOfMaximumChamber = 0.39; %OPTIMIZE
end

if ~isCallingFromOutside
    flap_gurneyLengthPercentage = 3; %OPTIMIZE
end
flap_gurneyLength = flap_gurneyLengthPercentage/100 / flap_size;
flap_gurneyThickness = 0.5/100 / flap_size;
flap_gurneyAdditionalAngle = deg2rad(-30*0);
flap_isClosedGurney = 0;



% slat element

hasSlat = 0;

slat_size = 0.2;
slat_positionX = -0.1;
slat_positionY = -0.05;
slat_rotation = deg2rad(60);

slat_thickness = 0.15;
slat_maximumChamber = 0.15;
slat_positionOfMaximumChamber = 0.4;

slat_gurneyLength = 0;
slat_gurneyThickness = 0;
slat_gurneyAdditionalAngle = 0;
slat_isClosedGurney = 0;



%% derive global positions of airfoil ensemble

% references and helpers

xShift = 0.25 * 0;
R = @(phi) [cos(phi) -sin(phi); sin(phi) cos(phi)];



% main elment

ss = 0 : 1/(200 - 1) : 1;
[main_x, main_y, main_x_spline, main_y_spline] = airfoil2(ss, main_yTE, main_x1, main_ch1, main_th1, main_tw1, main_x2, main_ch2, main_th2, main_tw2);
main_xy = [main_x', main_y']';
mainSpline_xy = [main_x_spline, main_y_spline]';

% ... move, scale, rotate for AoA

for i = 1 : size(main_xy, 2)
    main_xy(:, i) = R(-angleOfAttackWrtMain) * chord * (main_xy(:, i) - [xShift; 0]);
end
for i = 1 : size(mainSpline_xy, 2)
    mainSpline_xy(:, i) = R(-angleOfAttackWrtMain) * chord * (mainSpline_xy(:, i) - [xShift; 0]);
end



% flap element

[flap_x, flap_y] = NACA(flap_thickness, flap_maximumChamber, flap_positionOfMaximumChamber, flap_gurneyLength, flap_gurneyThickness, flap_gurneyAdditionalAngle, flap_isClosedGurney);
flap_xy = [flap_x, flap_y]';

% ... scale

flap_xy = flap_size * flap_xy;

% ... rotate

for i = 1 : size(flap_xy, 2)
    flap_xy(:, i) = R(flap_rotation) * flap_xy(:, i);
end

% ... move, scale, rotate for AoA

for i = 1 : size(flap_xy, 2)
    flap_xy(:, i) = flap_xy(:, i) + [flap_positionX; flap_positionY];
end

% ... rotate AoA

for i = 1 : size(flap_xy, 2)
    flap_xy(:, i) = R(-angleOfAttackWrtMain) * chord * (flap_xy(:, i) - [xShift; 0]);
end



% slat element

if hasSlat
    % airfoil
    
    [slat_x, slat_y] = NACA(slat_thickness, slat_maximumChamber, slat_positionOfMaximumChamber, slat_gurneyLength, slat_gurneyThickness, slat_gurneyAdditionalAngle, slat_isClosedGurney);
    slat_xy = [slat_x, slat_y]';

    % ... scale slat

    slat_xy = slat_size * slat_xy;

    % ... rotate slat

    for i = 1 : size(slat_xy, 2)
        slat_xy(:, i) = R(slat_rotation) * slat_xy(:, i);
    end

    % ... move slat

    for i = 1 : size(slat_xy, 2)
        slat_xy(:, i) = slat_xy(:, i) + [slat_positionX; slat_positionY];
    end

    % ... move, scale, rotate for AoA

    for i = 1 : size(slat_xy, 2)
        slat_xy(:, i) = R(-angleOfAttackWrtMain) * chord * (slat_xy(:, i) - [xShift; 0]);
    end
end



%% compute minimum distance between main and flap, and check

% compute

distanceBetweenMainAndFlapMinimum = inf;
for i = 1 : size(main_xy, 2)
    for j = 1 : size(flap_xy, 2)
        distance_ij = norm(main_xy(:, i) - flap_xy(:, j));
        if distance_ij < distanceBetweenMainAndFlapMinimum
            distanceBetweenMainAndFlapMinimum = distance_ij;
        end
    end
end



% check

distanceBetweenMainAndFlapMinimallyDesired = 10e-3;

if distanceBetweenMainAndFlapMinimum < distanceBetweenMainAndFlapMinimallyDesired
    disp('Penalty due to distanceBetweenMainAndFlapMinimum.');
    cost = cost + penaltyOffset + penaltyMultiplyer * (distanceBetweenMainAndFlapMinimum - distanceBetweenMainAndFlapMinimallyDesired)^2;
end



%% compute thicknesses of main element and checks

% compute upper and lower surface of main element

xMinimumOfMain = main_x(1);
indexOfXMinimumOfMain = 1;
for i = 2 : length(main_x)
    if xMinimumOfMain > main_x(i)
        indexOfXMinimumOfMain = i;
        xMinimumOfMain = main_x(i);
    end
end

mainUpper_xy = main_xy(:, indexOfXMinimumOfMain : -1 : 1);
mainLower_xy = main_xy(:, indexOfXMinimumOfMain + 1 : end);



% compute thickness

thicknessAtIndexUpper = nan(size(mainUpper_xy(1, :)));
for indexUpper = 1 : length(mainUpper_xy(1, :))
    for indexLower = 1 : length(mainLower_xy(1, :)) - 1
        if mainLower_xy(1, indexLower) <= mainUpper_xy(1, indexUpper) && mainUpper_xy(1, indexUpper) <= mainLower_xy(1, indexLower + 1)
            
            t = (mainUpper_xy(1, indexUpper) - mainLower_xy(1, indexLower)) / (mainLower_xy(1, indexLower + 1) - mainLower_xy(1, indexLower));
            newPoint = mainLower_xy(:, indexLower) + t * (mainLower_xy(:, indexLower + 1) - mainLower_xy(:, indexLower));
            yLowerAtXUpper = newPoint(2);
            
            thicknessAtIndexUpper(indexUpper) = mainUpper_xy(2, indexUpper) - yLowerAtXUpper;
            break;
        end
    end
end


% extract important thicknesses

thicknessMaximum = max(thicknessAtIndexUpper); %note: can be constrained > some value

indexUppderAtThicknessMaximum = find(thicknessMaximum == thicknessAtIndexUpper);

thicknessMinimum = min(thicknessAtIndexUpper); %note: should be constrained > 0
thicknessMinimumBehindPointOfMaximumThickness = min(thicknessAtIndexUpper(indexUppderAtThicknessMaximum : end)); %note: should be constrained > trailing edge thickness



% checks

thicknessPercentageMaximumOfMainElementMinimallyDesired = 0.20;
thicknessPercentageMaximumOfMainElementMaximallyDesired = 0.25;

if thicknessMaximum < thicknessPercentageMaximumOfMainElementMinimallyDesired * chord
    disp('Penalty due to thicknessPercentageMaximumOfMainElementMinimallyDesired.');
    cost = cost + penaltyOffset + penaltyMultiplyer * (thicknessMaximum - thicknessPercentageMaximumOfMainElementMinimallyDesired * chord)^2;
end

if thicknessMaximum > thicknessPercentageMaximumOfMainElementMaximallyDesired * chord
    disp('Penalty due to thicknessPercentageMaximumOfMainElementMaximallyDesired.');
    cost = cost + penaltyOffset + penaltyMultiplyer * (thicknessMaximum - thicknessPercentageMaximumOfMainElementMaximallyDesired * chord)^2;
end

if thicknessMinimum <= 0
    disp('Penalty due to thicknessMinimum.');
    cost = cost + penaltyOffset + penaltyMultiplyer * (thicknessMinimum - 0)^2;
end

%notes: this seems to have some numerical problems -- should be covered by monotony check anyways
% if thicknessMinimumBehindPointOfMaximumThickness + eps < main_yTE * chord
%     disp('Penalty due to thicknessMinimumBehindPointOfMaximumThickness.');
%     cost = cost + penaltyOffset + penaltyMultiplyer * (thicknessMinimumBehindPointOfMaximumThickness - main_yTE * chord)^2;
% end



%% thickness at special points

% points of interest

thicknessMinimumAtChordPercentages = [
    %chord [%]  %thickness minimum desired [%]  %computed thickness [%]
%   0.40        0.2                             nan %NOT NECESSARY (?)
    0.70        0.1                             nan
];



% compute

for index = 1 : size(thicknessMinimumAtChordPercentages, 1)
    xValue = (thicknessMinimumAtChordPercentages(index, 1) - xShift) * chord;
    for indexUpper = 1 : length(mainUpper_xy(1, :)) - 1
        if mainUpper_xy(1, indexUpper) <= xValue && xValue <= mainUpper_xy(1, indexUpper + 1)
            
            t = (xValue - mainUpper_xy(1, indexUpper)) / (mainUpper_xy(1, indexUpper + 1) - mainUpper_xy(1, indexUpper));
            thicknessAtX = thicknessAtIndexUpper(indexUpper) + t * (thicknessAtIndexUpper(indexUpper + 1) - thicknessAtIndexUpper(indexUpper));
            
            thicknessMinimumAtChordPercentages(index, 3) = thicknessAtX / chord;
            break;
        end
    end
end



% add panalty

for index = 1 : size(thicknessMinimumAtChordPercentages, 1)
    if thicknessMinimumAtChordPercentages(index, 2) > thicknessMinimumAtChordPercentages(index, 3)
        disp('Penalty due to thicknessMinimumAtChordPercentages.');
        cost = cost + penaltyOffset + penaltyMultiplyer * (thicknessMinimumAtChordPercentages(index, 2) - thicknessMinimumAtChordPercentages(index, 3))^2;
    end
end



%% compute if thickness is big enough to fit servo well

% servo dimensions, including wall of wing extrusions and servo 3D printed housing

thicknessServo = 21e-3 + 4e-3;
widthServo = 41e-3 + 4e-3;



% further settings

%chordPercentageMaximumWhenServoEndsDesired = 1/3; %seems way to conservative
chordPercentageMaximumWhenServoEndsDesired = 0.40;



% compute

widthAvailableForServo = 0;
chordPercentageMinimumWhenServoEnds = nan;

servo_xy = nan(5, 2);

for indexUpper1 = 1 : length(mainUpper_xy(1, :))
    if thicknessAtIndexUpper(indexUpper1) > thicknessServo
        for indexUpper2 = indexUpper1 : length(mainUpper_xy(1, :))
            if thicknessAtIndexUpper(indexUpper2) < thicknessServo
                widthAvailableForServo = mainUpper_xy(1, indexUpper2) - mainUpper_xy(1, indexUpper1);
                chordPercentageMinimumWhenServoEnds = (mainUpper_xy(1, indexUpper1) + widthServo) / chord + xShift;
                
                servo_xy = [
                    0   0
                    widthServo   0
                    widthServo   thicknessServo
                    0   thicknessServo
                    0   0
                ]';
                for i = 1 : 5
                    servo_xy(1, i) = servo_xy(1, i) + mainUpper_xy(1, indexUpper1);
                end
                
                break;
            end
        end
        
        break;
    end
end



% penalty

if widthAvailableForServo <= 0
    disp('Penalty due to widthAvailableForServo.');
    cost = cost + penaltyOffset;
elseif chordPercentageMinimumWhenServoEnds > chordPercentageMaximumWhenServoEndsDesired
    disp('Penalty due to chordPercentageMinimumWhenServoEnds.');
    cost = cost + penaltyOffset + penaltyMultiplyer * (chordPercentageMinimumWhenServoEnds - chordPercentageMaximumWhenServoEndsDesired)^2;
end



%% compute chamber line

chamberYAtIndexUpper = mainUpper_xy(2, :) - thicknessAtIndexUpper / 2;



%% compute ensemble chord line

positionChordLineLeadingEdge = [mainUpper_xy(1, 1); chamberYAtIndexUpper(2)]; %note: chamberYAtIndexUpper(1) is nan
positionChordLineTrailingEdge = flap_xy(:, end);
positionChordLineQuarter = positionChordLineLeadingEdge + (positionChordLineTrailingEdge - positionChordLineLeadingEdge) / 4;

chordOfEnsemble = norm(positionChordLineTrailingEdge - positionChordLineLeadingEdge);

temp = positionChordLineTrailingEdge - positionChordLineLeadingEdge;
angleOfIncidenceOfEnsemble = -atan2(temp(2), temp(1));



%% compute monotony behind point of maximum thickness

% compute

monotonyBehindPointOfMaximumThickness = 0;
for indexUpper = indexUppderAtThicknessMaximum : length(thicknessAtIndexUpper) - 1
    changeOfThickness = thicknessAtIndexUpper(indexUpper + 1) - thicknessAtIndexUpper(indexUpper);
    if changeOfThickness > 0
        monotonyBehindPointOfMaximumThickness = monotonyBehindPointOfMaximumThickness + changeOfThickness;
    end
end



% check (note: desired is to have negative monotony)

if monotonyBehindPointOfMaximumThickness > 0
    disp('Penalty due to monotonyBehindPointOfMaximumThickness.');
    cost = cost + penaltyOffset + penaltyMultiplyer * (monotonyBehindPointOfMaximumThickness - 0)^2;
end







%% plot

% plot ensemble, main element, and main element's thickness

if isEnablePlotting
    figure(1); clf;

    subplot(2, 2, 1);
    plot(mainUpper_xy(1, :), mainUpper_xy(2, :), 'r.-'); hold on;
    plot(mainLower_xy(1, :), mainLower_xy(2, :), 'b.-'); hold on;
    plot(mainUpper_xy(1, :), chamberYAtIndexUpper, 'g.-'); hold on;
    plot([positionChordLineLeadingEdge(1) positionChordLineTrailingEdge(1)], [positionChordLineLeadingEdge(2) positionChordLineTrailingEdge(2)], 'c'); hold on;
    plot(positionChordLineQuarter(1), positionChordLineQuarter(2), 'co');
    plot(positionOfAerodynamicForceOnChordLine(1), positionOfAerodynamicForceOnChordLine(2), 'cx');
    plot(flap_xy(1, :), flap_xy(2, :), 'r.-'); hold on;
    if hasSlat, plot(slat_xy(1, :), slat_xy(2, :), 'r.-'); hold on; end
    plot(positionOfAerodynamicForce(1), positionOfAerodynamicForce(2), 'kx');
    plot(positionOfAerodynamicForce(1) - distanceVerticalBeamsGap/2, positionOfAerodynamicForce(2), 'k*');
    plot(positionOfAerodynamicForce(1) + distanceVerticalBeamsGap/2, positionOfAerodynamicForce(2), 'k*');
    plot(positionOfFlapAxis(1), positionOfFlapAxis(2), 'ko');
    grid on;
    axis equal;
    title('ensemble')

    subplot(2, 2, 2);
    %plot((main_x_spline - xShift) * chord, main_y_spline * chord, 'ko-'); hold on; %TODO: rotate spline or rotate airfoil back
    plot(mainSpline_xy(1, :), mainSpline_xy(2, :), 'ko-'); hold on;
    plot(mainUpper_xy(1, :), mainUpper_xy(2, :), 'r.-'); hold on;
    plot(mainLower_xy(1, :), mainLower_xy(2, :), 'b.-'); hold on;
    plot(mainUpper_xy(1, :), chamberYAtIndexUpper, 'g.-'); hold on;
    plot(positionOfAerodynamicForce(1), positionOfAerodynamicForce(2), 'kx');
    plot(positionOfAerodynamicForce(1) - distanceVerticalBeamsGap/2, positionOfAerodynamicForce(2), 'k*');
    plot(positionOfAerodynamicForce(1) + distanceVerticalBeamsGap/2, positionOfAerodynamicForce(2), 'k*');
    grid on;
    axis equal;
    title('main element')

    subplot(2, 2, 4);
    plot(mainUpper_xy(1, :), thicknessAtIndexUpper, 'k.-'); hold on;
    for i = 1 : size(thicknessMinimumAtChordPercentages, 1)
        stem((thicknessMinimumAtChordPercentages(i, 1) - xShift) * chord, thicknessMinimumAtChordPercentages(i, 2) * chord, 'k')
    end
    plot(servo_xy(1, :), servo_xy(2, :), 'k-', 'LineWidth', 2);
    plot([1, 1] * chordPercentageMaximumWhenServoEndsDesired * chord, [0, thicknessServo], 'r-', 'LineWidth', 2);
    plot(positionOfAerodynamicForce(1), 0, 'kx');
    plot(positionOfAerodynamicForce(1) - distanceVerticalBeamsGap/2, 0, 'k*');
    plot(positionOfAerodynamicForce(1) + distanceVerticalBeamsGap/2, 0, 'k*');
    grid on;
    axis equal;
    title('main element thickness')
end



% add Makani airfoil

if 0
    main_Makani = [
       9.9791376912378293e-01	   2.7777777777775459e-03	
       8.4337814866326666e-01	   5.8333333333333182e-02	
       7.1113042806366855e-01	   1.0138888888888875e-01	
       5.5801653531138906e-01	   1.4583333333333320e-01	
       4.1603693401329001e-01	   1.8749999999999994e-01	
       3.0890511512903723e-01	   2.0138888888888878e-01	
       2.5188533456961831e-01	   1.9999999999999990e-01	
       2.1575104311543811e-01	   1.9027777777777760e-01	
       1.4488100757224542e-01	   1.6805555555555540e-01	
       7.5463606861381524e-02	   1.2361111111111095e-01	
       3.3843300880853001e-02	   8.6111111111110916e-02	
       6.1427909133052910e-04	   3.1944444444444220e-02	
      -2.0591871426363807e-03	  -6.9444444444446973e-03	
       2.8647040642868150e-02	  -4.5833333333333615e-02	
       7.5985164580435763e-02	  -6.3888888888889050e-02	
       1.6224694792149594e-01	  -7.5000000000000233e-02	
       2.3736671302735279e-01	  -8.0555555555555769e-02	
       3.5698114665430397e-01	  -8.1944444444444819e-02	
       4.4040720135991346e-01	  -7.3611111111111321e-02	
       5.3354195642095492e-01	  -5.5555555555555691e-02	
       6.4334337814866305e-01	  -2.9166666666666924e-02	
       7.4064673157162697e-01	  -9.7222222222223820e-03	
       8.5047906042342736e-01	   5.5555555555554248e-03	
       9.2418482460207063e-01	   8.3333333333330817e-03	
       9.9791763251429433e-01	   1.3888888888886619e-03	
    ]';
    flap_Makani = [
       1.2766998918250654e+00	  -2.2083333333333349e-01	
       1.2083565136764023e+00	  -1.5138888888888902e-01	
       1.1525575645186215e+00	  -9.1666666666666785e-02	
       1.0995827538247565e+00	  -4.7222222222222471e-02	
       1.0383132437026732e+00	  -2.0833333333333565e-02	
       9.9937413073713477e-01	  -2.2222222222222421e-02	
       9.8831324370267348e-01	  -4.5833333333333615e-02	
       1.0092682738371193e+00	  -7.9166666666666829e-02	
       1.0816875289754289e+00	  -1.1388888888888901e-01	
       1.1527391438726624e+00	  -1.5694444444444483e-01	
       1.2265646731571627e+00	  -1.9722222222222235e-01	
       1.2753129346314322e+00	  -2.2222222222222243e-01	
    ]';

    if isEnablePlotting
        subplot(2, 2, 1);
        plot((main_Makani(1, :) - xShift) * chord, main_Makani(2, :) * chord, 'b-.'); hold on;
        plot((flap_Makani(1, :) - xShift) * chord, flap_Makani(2, :) * chord, 'b-.'); hold on;
    end
end



% add previously optimized

if 1
    main_other = [
        0.11249,-0.0017134
        0.11116,-0.0006362
        0.10979,0.00043816
        0.10835,0.001509
        0.10686,0.0025755
        0.10532,0.003637
        0.10373,0.0046929
        0.10209,0.0057424
        0.10041,0.0067849
        0.098679,0.0078196
        0.096905,0.0088459
        0.09509,0.009863
        0.093236,0.01087
        0.091343,0.011867
        0.089414,0.012852
        0.08745,0.013826
        0.085453,0.014787
        0.083425,0.015734
        0.081367,0.016668
        0.079282,0.017586
        0.07717,0.01849
        0.075035,0.019377
        0.072876,0.020247
        0.070697,0.0211
        0.068498,0.021935
        0.066282,0.02275
        0.06405,0.023546
        0.061804,0.024322
        0.059545,0.025077
        0.057276,0.02581
        0.054998,0.02652
        0.052713,0.027208
        0.050422,0.027871
        0.048128,0.02851
        0.045831,0.029124
        0.043534,0.029712
        0.041238,0.030274
        0.038946,0.030808
        0.036658,0.031314
        0.034377,0.031791
        0.032104,0.032239
        0.029841,0.032657
        0.02759,0.033044
        0.025352,0.0334
        0.023129,0.033723
        0.020923,0.034014
        0.018735,0.034271
        0.016568,0.034494
        0.014422,0.034681
        0.0123,0.034833
        0.010204,0.034949
        0.0081348,0.035028
        0.0060944,0.035069
        0.0040844,0.035071
        0.0021067,0.035034
        0.000163,0.034957
        -0.0017451,0.03484
        -0.0036157,0.034682
        -0.0054473,0.034481
        -0.007238,0.034238
        -0.0089862,0.033952
        -0.01069,0.033621
        -0.012348,0.033246
        -0.013958,0.032825
        -0.015519,0.032358
        -0.017029,0.031845
        -0.018486,0.031283
        -0.019889,0.030674
        -0.021236,0.030018
        -0.022529,0.029318
        -0.023768,0.028577
        -0.024953,0.027797
        -0.026084,0.02698
        -0.027162,0.026131
        -0.028188,0.02525
        -0.02916,0.024342
        -0.030081,0.023408
        -0.030949,0.022452
        -0.031766,0.021475
        -0.032532,0.020481
        -0.033247,0.019472
        -0.033911,0.018451
        -0.034525,0.01742
        -0.03509,0.016383
        -0.035604,0.015342
        -0.03607,0.014299
        -0.036486,0.013257
        -0.036854,0.012219
        -0.037174,0.011188
        -0.037446,0.010166
        -0.037671,0.0091552
        -0.037848,0.0081593
        -0.037978,0.0071805
        -0.038062,0.0062215
        -0.038099,0.005285
        -0.038091,0.0043736
        -0.038037,0.00349
        -0.037938,0.0026368
        -0.037794,0.0018167
        -0.037605,0.0010324
        -0.037372,0.00028641
        -0.037096,-0.00042024
        -0.036775,-0.0010882
        -0.03641,-0.0017183
        -0.036002,-0.0023112
        -0.03555,-0.0028675
        -0.035055,-0.0033882
        -0.034516,-0.0038737
        -0.033934,-0.004325
        -0.033309,-0.0047428
        -0.03264,-0.0051276
        -0.031929,-0.0054804
        -0.031174,-0.0058018
        -0.030377,-0.0060925
        -0.029536,-0.0063532
        -0.028653,-0.0065848
        -0.027728,-0.0067879
        -0.02676,-0.0069632
        -0.02575,-0.0071114
        -0.024697,-0.0072334
        -0.023602,-0.0073298
        -0.022465,-0.0074013
        -0.021286,-0.0074487
        -0.020065,-0.0074727
        -0.018802,-0.0074741
        -0.017497,-0.0074535
        -0.016151,-0.0074117
        -0.014763,-0.0073493
        -0.013334,-0.0072673
        -0.011863,-0.0071662
        -0.010351,-0.0070467
        -0.0087976,-0.0069097
        -0.0072033,-0.0067558
        -0.005568,-0.0065858
        -0.0038924,-0.0064005
        -0.002178,-0.0062006
        -0.00042656,-0.0059871
        0.0013603,-0.0057607
        0.003181,-0.0055224
        0.0050337,-0.005273
        0.006917,-0.0050134
        0.0088291,-0.0047444
        0.010768,-0.0044668
        0.012733,-0.0041815
        0.014722,-0.0038893
        0.016733,-0.0035912
        0.018764,-0.003288
        0.020814,-0.0029804
        0.022882,-0.0026694
        0.024965,-0.0023559
        0.027062,-0.0020406
        0.029171,-0.0017244
        0.031291,-0.0014082
        0.03342,-0.0010928
        0.035557,-0.00077911
        0.037699,-0.00046794
        0.039845,-0.00016015
        0.041993,0.00014341
        0.044142,0.00044188
        0.046291,0.00073442
        0.048436,0.0010202
        0.050578,0.0012983
        0.052713,0.0015679
        0.054841,0.0018282
        0.05696,0.0020782
        0.059068,0.0023172
        0.061163,0.0025444
        0.063244,0.0027587
        0.06531,0.0029595
        0.067358,0.0031458
        0.069387,0.0033168
        0.071395,0.0034716
        0.07338,0.0036095
        0.075342,0.0037294
        0.077278,0.0038307
        0.079186,0.0039123
        0.081066,0.0039736
        0.082914,0.0040136
        0.084731,0.0040314
        0.086513,0.0040263
        0.08826,0.0039974
        0.08997,0.0039437
        0.09164,0.0038646
        0.09327,0.003759
        0.094858,0.0036262
        0.096402,0.0034653
        0.0979,0.0032755
        0.099351,0.0030559
        0.10075,0.0028057
        0.1021,0.0025239
        0.1034,0.0022098
        0.10465,0.0018625
        0.10584,0.0014812
        0.10697,0.0010649
        0.10805,0.00061289
        0.10906,0.00012427
        0.11001,-0.00040182
        0.1109,-0.00096622
        0.11172,-0.0015698
        0.11248,-0.0022134
    ]';
    flap_other = [
        0.16458,-0.059564
        0.16457,-0.059547
        0.16454,-0.059495
        0.1645,-0.059408
        0.16443,-0.059287
        0.16435,-0.059132
        0.16424,-0.058943
        0.16412,-0.05872
        0.16398,-0.058464
        0.16382,-0.058175
        0.16365,-0.057853
        0.16345,-0.057499
        0.16323,-0.057114
        0.163,-0.056698
        0.16274,-0.056251
        0.16247,-0.055775
        0.16218,-0.055269
        0.16186,-0.054736
        0.16153,-0.054175
        0.16118,-0.053587
        0.16081,-0.052974
        0.16042,-0.052335
        0.16002,-0.051673
        0.15959,-0.050987
        0.15914,-0.05028
        0.15868,-0.049551
        0.15819,-0.048802
        0.15769,-0.048034
        0.15717,-0.047247
        0.15662,-0.046444
        0.15607,-0.045625
        0.15549,-0.044791
        0.1549,-0.043943
        0.15428,-0.043083
        0.15366,-0.042211
        0.15301,-0.041329
        0.15235,-0.040438
        0.15167,-0.039539
        0.15098,-0.038633
        0.15027,-0.037722
        0.14954,-0.036807
        0.14881,-0.035888
        0.14805,-0.034968
        0.14729,-0.034047
        0.14651,-0.033127
        0.14572,-0.032208
        0.14492,-0.031293
        0.14411,-0.030382
        0.14328,-0.029477
        0.14245,-0.028578
        0.14161,-0.027688
        0.14076,-0.026806
        0.1399,-0.025936
        0.13904,-0.025076
        0.13817,-0.02423
        0.13729,-0.023398
        0.13641,-0.02258
        0.13551,-0.021765
        0.13459,-0.020959
        0.13365,-0.020177
        0.13272,-0.01942
        0.13177,-0.018688
        0.13083,-0.017984
        0.12988,-0.017306
        0.12894,-0.016657
        0.12799,-0.016037
        0.12705,-0.015446
        0.12611,-0.014885
        0.12518,-0.014354
        0.12426,-0.013854
        0.12334,-0.013386
        0.12244,-0.012948
        0.12154,-0.012541
        0.12066,-0.012166
        0.1198,-0.011822
        0.11895,-0.01151
        0.11812,-0.011228
        0.1173,-0.010977
        0.11651,-0.010757
        0.11574,-0.010567
        0.11499,-0.010407
        0.11426,-0.010276
        0.11357,-0.010174
        0.11289,-0.010101
        0.11225,-0.010055
        0.11163,-0.010037
        0.11104,-0.010045
        0.11048,-0.01008
        0.10995,-0.010141
        0.10945,-0.010227
        0.10899,-0.010337
        0.10856,-0.010471
        0.10816,-0.010629
        0.10779,-0.01081
        0.10746,-0.011014
        0.10716,-0.01124
        0.10689,-0.011487
        0.10666,-0.011756
        0.10646,-0.012046
        0.1063,-0.012357
        0.10617,-0.012684
        0.10608,-0.013023
        0.10603,-0.013373
        0.10601,-0.013735
        0.10603,-0.014107
        0.10608,-0.01449
        0.10617,-0.014883
        0.10629,-0.015287
        0.10645,-0.0157
        0.10664,-0.016122
        0.10687,-0.016554
        0.10712,-0.016994
        0.10741,-0.017443
        0.10773,-0.0179
        0.10808,-0.018364
        0.10846,-0.018836
        0.10887,-0.019315
        0.10931,-0.019801
        0.10978,-0.020294
        0.11027,-0.020793
        0.11079,-0.021298
        0.11134,-0.021809
        0.11191,-0.022326
        0.11251,-0.022849
        0.11312,-0.023377
        0.11376,-0.02391
        0.11443,-0.024449
        0.11511,-0.024994
        0.11581,-0.025543
        0.11653,-0.026098
        0.11727,-0.026658
        0.11803,-0.027224
        0.1188,-0.027795
        0.11959,-0.028372
        0.12039,-0.028954
        0.12121,-0.029542
        0.12204,-0.030137
        0.12288,-0.030737
        0.12373,-0.031343
        0.1246,-0.031956
        0.12547,-0.032575
        0.12635,-0.033201
        0.12725,-0.033846
        0.12818,-0.034507
        0.12912,-0.035169
        0.13006,-0.035831
        0.13102,-0.036494
        0.13197,-0.037157
        0.13293,-0.03782
        0.1339,-0.038482
        0.13486,-0.039143
        0.13583,-0.039803
        0.13679,-0.040461
        0.13776,-0.041117
        0.13872,-0.041771
        0.13968,-0.042421
        0.14063,-0.043068
        0.14157,-0.043711
        0.14251,-0.044349
        0.14345,-0.044983
        0.14437,-0.04561
        0.14528,-0.046232
        0.14618,-0.046847
        0.14707,-0.047454
        0.14795,-0.048053
        0.14881,-0.048644
        0.14966,-0.049225
        0.15049,-0.049796
        0.15131,-0.050357
        0.15211,-0.050907
        0.15289,-0.051445
        0.15365,-0.05197
        0.15439,-0.052483
        0.15511,-0.052982
        0.15581,-0.053466
        0.15649,-0.053936
        0.15715,-0.05439
        0.15778,-0.054828
        0.15838,-0.05525
        0.15897,-0.055654
        0.15952,-0.056041
        0.16005,-0.05641
        0.16056,-0.056761
        0.16103,-0.057092
        0.16148,-0.057404
        0.1619,-0.057696
        0.16229,-0.057969
        0.16265,-0.05822
        0.16298,-0.058451
        0.16328,-0.058661
        0.16355,-0.058849
        0.16379,-0.059016
        0.16396,-0.059135
        0.16139,-0.062827
        0.16201,-0.063256
        0.16458,-0.059564
    ]';

    if isEnablePlotting
        subplot(2, 2, 1);
        plot((main_other(1, :) + 0.15 / 4), main_other(2, :), '--', 'Color', [0.5 0.5 0.5]); hold on;
        plot((flap_other(1, :) + 0.15 / 4), flap_other(2, :), '--', 'Color', [0.5 0.5 0.5]); hold on;
    end
end



% add S1223

if 1
    main_other = [
          1.00000     0.00000
          0.99838     0.00126
          0.99417     0.00494
          0.98825     0.01037
          0.98075     0.01646
          0.97111     0.02250
          0.95884     0.02853
          0.94389     0.03476
          0.92639     0.04116
          0.90641     0.04768
          0.88406     0.05427
          0.85947     0.06089
          0.83277     0.06749
          0.80412     0.07402
          0.77369     0.08044
          0.74166     0.08671
          0.70823     0.09277
          0.67360     0.09859
          0.63798     0.10412
          0.60158     0.10935
          0.56465     0.11425
          0.52744     0.11881
          0.49025     0.12303
          0.45340     0.12683
          0.41721     0.13011
          0.38193     0.13271
          0.34777     0.13447
          0.31488     0.13526
          0.28347     0.13505
          0.25370     0.13346
          0.22541     0.13037
          0.19846     0.12594
          0.17286     0.12026
          0.14863     0.11355
          0.12591     0.10598
          0.10482     0.09770
          0.08545     0.08879
          0.06789     0.07940
          0.05223     0.06965
          0.03855     0.05968
          0.02694     0.04966
          0.01755     0.03961
          0.01028     0.02954
          0.00495     0.01969
          0.00155     0.01033
          0.00005     0.00178
          0.00044    -0.00561
          0.00264    -0.01120
          0.00789    -0.01427
          0.01718    -0.01550
          0.03006    -0.01584
          0.04627    -0.01532
          0.06561    -0.01404
          0.08787    -0.01202
          0.11282    -0.00925
          0.14020    -0.00563
          0.17006    -0.00075
          0.20278     0.00535
          0.23840     0.01213
          0.27673     0.01928
          0.31750     0.02652
          0.36044     0.03358
          0.40519     0.04021
          0.45139     0.04618
          0.49860     0.05129
          0.54639     0.05534
          0.59428     0.05820
          0.64176     0.05976
          0.68832     0.05994
          0.73344     0.05872
          0.77660     0.05612
          0.81729     0.05219
          0.85500     0.04706
          0.88928     0.04088
          0.91966     0.03387
          0.94573     0.02624
          0.96693     0.01822
          0.98255     0.01060
          0.99268     0.00468
          0.99825     0.00115
          1.00000     0.00000
    ]';
    for i = 1 : size(main_other, 2)
        main_other(:, i) = R(-deg2rad(17)) * 0.21 * main_other(:, i);
    end
    main_other(1, :) = main_other(1, :) - 0.00;
    main_other(2, :) = main_other(2, :) + 0.02;

    if isEnablePlotting
        subplot(2, 2, 1);
        plot(main_other(1, :), main_other(2, :), '-.', 'Color', [.5 .5 .5]); hold on;
    end
end



% add previously optimized Ziemkewitz

if 0    
    main_Z = [0.0749943706948113,0.0749715181156618,0.0749029720916880,0.0747887764213641,0.0746290095107223,0.0744237841658753,0.0741732474240817,0.0738775803687702,0.0735369979136907,0.0731517485490964,0.0727221140450934,0.0722484091078840,0.0717309809845708,0.0711702090117858,0.0705665041027960,0.0699203081669392,0.0692320934542725,0.0685023618171615,0.0677316438791636,0.0669204980999252,0.0660695097228607,0.0651792895900182,0.0642504728056574,0.0632837172265061,0.0622797017522111,0.0612391243838642,0.0601627000112536,0.0590511578801124,0.0579052386782955,0.0567256911634086,0.0555132682323035,0.0542687223026862,0.0529927998353613,0.0516862347670834,0.0503497405404983,0.0489840002964909,0.0475896546149029,0.0461672859180859,0.0447173982302050,0.0432403903113141,0.0417365190710481,0.0402058482509242,0.0386481739135978,0.0370629117075032,0.0354489174878517,0.0338041831726390,0.0321252761194695,0.0304061783373767,0.0286354193294442,0.0267863752729027,0.0247141884925106,0.0226426409566000,0.0207952818585745,0.0190271420923706,0.0173115385692589,0.0156369557177146,0.0139973364777137,0.0123892119739093,0.0108105392342136,0.00926013979864606,0.00773739517312135,0.00624206725267252,0.00477418551116074,0.00333397249321378,0.00192179254981031,0.000538115338664885,-0.000816510933161939,-0.00214147962775386,-0.00343614114293448,-0.00469981503865476,-0.00593179971232215,-0.00713138024000746,-0.00829783482056341,-0.00943044013808581,-0.0105284758742922,-0.0115912285435205,-0.0126179947810321,-0.0136080841848664,-0.0145608217891160,-0.0154755502298170,-0.0163516316520422,-0.0171884493971063,-0.0179854095012417,-0.0187419420310812,-0.0194575022763592,-0.0201315718160598,-0.0207636594704921,-0.0213533021481448,-0.0219000655922826,-0.0224035450275076,-0.0228633656999508,-0.0232791832945303,-0.0236506841950656,-0.0239775855197565,-0.0242596347961223,-0.0244966089848450,-0.0246883121695088,-0.0248345700674692,-0.0249352152072927,-0.0249900330415010,-0.0249981235649371,-0.0249620747909169,-0.0248832740549530,-0.0247605301969664,-0.0245936411969542,-0.0243825803078603,-0.0241274089270634,-0.0238282489430724,-0.0234852707453849,-0.0230986869725940,-0.0226687487885701,-0.0221957434109963,-0.0216799923093245,-0.0211218497765601,-0.0205217017118358,-0.0198799645169105,-0.0191970840447615,-0.0184735345576875,-0.0177098176630715,-0.0169064612007570,-0.0160640180587475,-0.0151830648946405,-0.0142642007393416,-0.0133080454573958,-0.0123152380347301,-0.0112864346595676,-0.0102223065554478,-0.00912353751615700,-0.00799082108019517,-0.00682485726607690,-0.00562634876770009,-0.00439599647886314,-0.00313449417429765,-0.00184252211605505,-0.000520739270658251,0.000830226298601573,0.00220978947893821,0.00361742348065115,0.00505267984145846,0.00651521504143261,0.00800482681888910,0.00952150518850920,0.0110655066041955,0.0126374662600786,0.0142385768697073,0.0158708918703016,0.0175378843428109,0.0192456041554311,0.0210055354967476,0.0228442621034359,0.0249061112410624,0.0269686371835782,0.0288090725602918,0.0305715913480426,0.0322826962942024,0.0339538102407558,0.0355909299503175,0.0371974791913538,0.0387754637585295,0.0403260292757022,0.0418497637906175,0.0433468762284251,0.0448173085379637,0.0462608098084337,0.0476769873127873,0.0490653428968836,0.0504252997000448,0.0517562222864013,0.0530574321578573,0.0543282199490227,0.0555678551850509,0.0567755942132109,0.0579506867405846,0.0590923812897405,0.0601999298011761,0.0612725915530702,0.0623096365273810,0.0633103483213117,0.0642740266811724,0.0651999897193353,0.0660875758627110,0.0669361455718281,0.0677450828624075,0.0685137966557047,0.0692417219794565,0.0699283210377162,0.0705730841649774,0.0711755306776162,0.0717352096337074,0.0722517005106072,0.0727246138082693,0.0731535915850289,0.0735383079315011,0.0738784693872503,0.0741738153039542,0.0744241181578220,0.0746291838128983,0.0747888517352377,0.0749029951547651,0.0749715211630632,0.0749943706948113;0.000918892860580999,0.000930106832122252,0.000964186387603382,0.00102173585408332,0.00110331442413724,0.00120944315203926,0.00134060476892396,0.00149724070911413,0.00167974647477108,0.00188846581148256,0.00212368395925671,0.00238562016520791,0.00267441961467466,0.00299014492935591,0.00333276738338254,0.00370215799599834,0.00409807866999274,0.00452017355649914,0.00496796083799108,0.00544082513122078,0.00593801071949445,0.00645861582816489,0.00700158815768958,0.00756572188426337,0.00814965632818952,0.00875187647423298,0.00937071550580498,0.0100043594857756,0.0106508542810752,0.0113081147864072,0.0119739364550812,0.0126460090932922,0.0133219328196730,0.0139992360366438,0.0146753952065921,0.0153478561775452,0.0160140567641094,0.0166714502659991,0.0173175296072212,0.0179498518178959,0.0185660626821217,0.0191639215852239,0.0197413270061523,0.0202963439229260,0.0208272361275442,0.0213325103946505,0.0218109895442803,0.0222619615692772,0.0226855647458904,0.0230841905872409,0.0234793278920144,0.0238222940028976,0.0240834039172307,0.0242932407304684,0.0244590328780584,0.0245845973479159,0.0246724107085062,0.0247242528899996,0.0247414781805627,0.0247251509666320,0.0246761227281874,0.0245950804948862,0.0244825803668996,0.0243390727883650,0.0241649230459957,0.0239604288416468,0.0237258359161567,0.0234613522148504,0.0231671608083936,0.0228434316276675,0.0224903319900743,0.0221080358596073,0.0216967317772505,0.0212566294108993,0.0207879646975117,0.0202910035792833,0.0197660443663193,0.0192134187875858,0.0186334918175661,0.0180266603862610,0.0173933510935161,0.0167340170538705,0.0160491339938825,0.0153391957085238,0.0146047089542576,0.0138461878097694,0.0130641474640998,0.0122590972841005,0.0114315328472509,0.0105819263627247,0.00971071446907729,0.00881828164510276,0.00790493609878699,0.00697087234051002,0.00601610911100015,0.00504037875364043,0.00404291212578212,0.00302196835813178,0.00197360705845989,0.000887276393762775,-0.000306297620193603,-0.00139450366812106,-0.00226550960026949,-0.00302070944782519,-0.00368354904761712,-0.00426593892328930,-0.00477535120503876,-0.00521705845161450,-0.00559508996019719,-0.00591271546785802,-0.00617271810687749,-0.00637756056785846,-0.00652949184909050,-0.00663061847866822,-0.00668295323019054,-0.00668844889422339,-0.00664902174049830,-0.00656656765133233,-0.00644297293045584,-0.00628012119420245,-0.00607989737505044,-0.00584418962247138,-0.00557488972171593,-0.00527389253672016,-0.00494309489956411,-0.00458439430357723,-0.00419968770223079,-0.00379087066608134,-0.00335983710175691,-0.00290847968812418,-0.00243869113410547,-0.00195236630961420,-0.00145140524577554,-0.000937716943403259,-0.000413223870203535,0.000120133032123996,0.000660389069825452,0.00120555054916824,0.00175358992774818,0.00230244181686904,0.00285000040360108,0.00339411905307698,0.00393261320863208,0.00446326843552179,0.00498385705896261,0.00549217069118975,0.00598608616669270,0.00646371386967896,0.00692379809429040,0.00736721890723713,0.00781577627777878,0.00820909700549124,0.00851305854281966,0.00876197133997204,0.00896332531137513,0.00912087148254172,0.00923705438785120,0.00931377544266054,0.00935271455088309,0.00935548720802457,0.00932372677666243,0.00925912770393127,0.00916346619274489,0.00903860687102992,0.00888650033192328,0.00870917457621460,0.00850872239466546,0.00828728615336958,0.00804704108847795,0.00779017797713287,0.00751888587731147,0.00723533549273660,0.00694166360512314,0.00663995891681543,0.00633224955823673,0.00602049243457873,0.00570656451405396,0.00539225609565821,0.00507926603783578,0.00476919888084304,0.00446356375502216,0.00416377493454882,0.00387115387126353,0.00358693252554073,0.00331225780025622,0.00304819687914000,0.00279574327143692,0.00255582337009958,0.00232930333998306,0.00211699616504199,0.00191966869881515,0.00173804858020886,0.00157283089679494,0.00142468450119298,0.00129425791447086,0.00118218478827501,0.00108908895503103,0.00101558919999774,0.000962304121423627,0.000929858122009220,0.000918892860581001];
    flap_Z = [0.109735038379615,0.109728153919413,0.109707498288153,0.109673064801938,0.109624842546852,0.109562816714407,0.109486969064415,0.109397278508648,0.109293721807009,0.109176274366366,0.109044911130880,0.108899607551476,0.108740340621202,0.108567089962518,0.108379838952103,0.108178575868589,0.107963295048695,0.107733998037537,0.107490694719463,0.107233404416529,0.106962156942745,0.106676993603378,0.106377968129967,0.106065147543154,0.105738612937020,0.105398460180262,0.105044800531185,0.104677761165190,0.104297485615025,0.103904134125645,0.103497883926950,0.103078929429010,0.102647482345489,0.102203771751992,0.101748044086788,0.101280563101914,0.100801609772978,0.100311482176089,0.0998104953401592,0.0992989810825091,0.0987772878351072,0.0982457804680141,0.0977048401156759,0.0971548640106163,0.0966011793766565,0.0960431615273800,0.0954789777179297,0.0949090761185130,0.0943339122139638,0.0937539486664520,0.0931696552271639,0.0925815086936104,0.0919899929078561,0.0913955987896714,0.0907988243974282,0.0902001750084943,0.0896001632099606,0.0889993089897751,0.0883981398177684,0.0877971907056504,0.0871970042348449,0.0865981305410096,0.0860011272442680,0.0854065593145494,0.0848149988619949,0.0842270248431227,0.0836432226743545,0.0830641837455630,0.0824905048274978,0.0819227873682617,0.0813616366754268,0.0808076609818671,0.0802614703949396,0.0797236757302216,0.0791948872326085,0.0786757131891556,0.0781667584395975,0.0776686227919663,0.0771818993521563,0.0767071727776006,0.0762450174664456,0.0757959956946973,0.0753606557147665,0.0749395298296391,0.0745331324575438,0.0741419582024633,0.0737664799461427,0.0734071469773790,0.0730643831743340,0.0727385852553960,0.0724301211137308,0.0721393282501109,0.0718665123179045,0.0716119457932459,0.0713758667824124,0.0711584779773066,0.0709599457686997,0.0707803995255480,0.0706199310472586,0.0704785941942796,0.0703581803442001,0.0702604568761970,0.0701854056914472,0.0701329987770402,0.0701031980159757,0.0700959548900177,0.0701112100801938,0.0701488929710736,0.0702089210662302,0.0702911993234795,0.0703956194195753,0.0705220589550192,0.0706703806104926,0.0708404312671391,0.0710320411034944,0.0712450226822866,0.0714791700405892,0.0717342577969093,0.0720100402887264,0.0723062507537644,0.0726226005678756,0.0729587785518564,0.0733144503587899,0.0736892579526378,0.0740828191877941,0.0744947274981665,0.0749245517030914,0.0753718359360302,0.0758360997005463,0.0763168380565541,0.0768135219382735,0.0773255986037424,0.0778524922141579,0.0783936045397512,0.0789483157873783,0.0795159855435546,0.0800959538252847,0.0806875422297778,0.0812900551730003,0.0819027812060225,0.0825249943972839,0.0831559557682448,0.0837949147694158,0.0844411107834813,0.0850937746421499,0.0857521301434881,0.0864153955568141,0.0870827851027464,0.0877535103967041,0.0884267818450385,0.0891018099840181,0.0897778067530746,0.0904539866950311,0.0911295680774479,0.0918037739307139,0.0924734298294635,0.0931351106848242,0.0937925237130187,0.0944449581443419,0.0950917126676949,0.0957320959125506,0.0963654268544376,0.0969910351498594,0.0976082614078879,0.0982164574068285,0.0988149862653248,0.0994032225780480,0.0999805525266584,0.100546373977043,0.101100096573906,0.101641141843608,0.102168943315757,0.102682946673381,0.103182609940671,0.103667403716208,0.104136811458347,0.104590329828031,0.105027469092795,0.105447753594113,0.105850722278575,0.106235929291686,0.106602944631420,0.106951354857044,0.107280763847165,0.107590793599552,0.107881085064007,0.108151298998436,0.108401116837388,0.108630241561607,0.108838398556707,0.109025336448819,0.109190827905108,0.109306535914844,0.107440376977111,0.107868882134537,0.109735038379615;-0.0361335113627648,-0.0361224550636948,-0.0360893038939496,-0.0360341109787756,-0.0359569646703293,-0.0358579882647921,-0.0357373396102979,-0.0355952106095782,-0.0354318266222641,-0.0352474457727479,-0.0350423581704005,-0.0348168850497419,-0.0345713778388646,-0.0343062171649958,-0.0340218118065542,-0.0337185976013925,-0.0333970363211122,-0.0330576145213933,-0.0327008423781841,-0.0323272525193493,-0.0319373988609768,-0.0315318554569962,-0.0311112153700686,-0.0306760895708747,-0.0302271058719718,-0.0297649079013070,-0.0292901541193000,-0.0288035168821370,-0.0283056815525906,-0.0277973456582990,-0.0272792180960422,-0.0267520183791559,-0.0262164759238611,-0.0256733293689766,-0.0251233259222600,-0.0245672207255082,-0.0240057762295771,-0.0234397615696658,-0.0228699519305879,-0.0222971278913323,-0.0217220747380250,-0.0211455817344453,-0.0205684413395470,-0.0199914483619780,-0.0194187589785831,-0.0188485939098463,-0.0182797759310676,-0.0177130603263718,-0.0171492030284920,-0.0165889594198359,-0.0160330830286698,-0.0154823241229833,-0.0149374282067208,-0.0143991344251944,-0.0138681738885977,-0.0133452679245804,-0.0128311262727823,-0.0123264452360283,-0.0118319058045128,-0.0113481717707316,-0.0108758878541045,-0.0104156778551617,-0.00996814285981176,-0.00953385951454400,-0.00911337839344435,-0.00870722247759496,-0.00831588576679354,-0.00793983204256338,-0.00757949380013747,-0.00723527136550516,-0.00690753221172146,-0.00659661048652428,-0.00630280676090889,-0.00602638800570494,-0.00576758780042639,-0.00552660677575709,-0.00530361328803804,-0.00509874432108075,-0.00491210660759027,-0.00474377795948877,-0.00459380879353161,-0.00446222383584921,-0.00434902398647274,-0.00425418832255252,-0.00417767621689217,-0.00411942954663542,-0.00407937496548507,-0.00405742621173168,-0.00405348642364311,-0.00406745043342935,-0.00409920701105961,-0.00414864102967293,-0.00421563552518499,-0.00430007362394472,-0.00440184031391655,-0.00452082403683825,-0.00465691808110170,-0.00481002175769339,-0.00498004134437579,-0.00516689078634939,-0.00536786128313568,-0.00558023353083261,-0.00580388995998597,-0.00603868430938483,-0.00628444230512618,-0.00654096257501076,-0.00680801778695831,-0.00708535599719582,-0.00737270219125012,-0.00766975999830490,-0.00797621355729915,-0.00829172951128031,-0.00861595910500795,-0.00894854035965266,-0.00928910029766688,-0.00963725719052728,-0.00999262280206637,-0.0103548046005202,-0.0107234079132110,-0.0110980379989436,-0.0114783020146983,-0.0118638108550342,-0.0122541808447281,-0.0126490352675513,-0.0130480057166717,-0.0134507332549328,-0.0138568693761528,-0.0142660767615681,-0.0146780298285571,-0.0150924150717915,-0.0155089311999150,-0.0159272890737021,-0.0163472114543654,-0.0167684325732125,-0.0171906975361681,-0.0176137615787476,-0.0180373891888567,-0.0184613531162863,-0.0188854332889457,-0.0193094156567255,-0.0197330909843901,-0.0201562536150787,-0.0205787002258289,-0.0210002285960638,-0.0214206364091876,-0.0218397201063608,-0.0222572738101810,-0.0226730883344144,-0.0230869502941350,-0.0234986413286656,-0.0239079374476152,-0.0243146085081065,-0.0247184178290188,-0.0251191219457865,-0.0255164705070072,-0.0259084400401374,-0.0262942684559621,-0.0266776494967475,-0.0270582857638367,-0.0274358669221731,-0.0278100705251302,-0.0281805629654469,-0.0285470005380476,-0.0289090305997099,-0.0292662928099985,-0.0296184204376011,-0.0299650417161810,-0.0303057812340749,-0.0306402613426123,-0.0309681035684799,-0.0312889300163913,-0.0316023647493130,-0.0319080351346248,-0.0322055731458259,-0.0324946166107101,-0.0327748103982966,-0.0330458075381981,-0.0333072702674963,-0.0335588710015704,-0.0338002932266493,-0.0340312323131241,-0.0342513962498429,-0.0344605063006983,-0.0346582975858034,-0.0348445195904126,-0.0350189366054911,-0.0351813281044408,-0.0353314890609754,-0.0354692302134797,-0.0355943782814039,-0.0357067761393301,-0.0358062829543123,-0.0358758572579792,-0.0389795245914518,-0.0392371742179985,-0.0361335113627648];

    if isEnablePlotting
        subplot(2, 2, 1);
        plot((main_Z(1, :)) * 1.2, main_Z(2, :) * 1.2, 'g-.'); hold on;
        plot((flap_Z(1, :)) * 1.2, flap_Z(2, :) * 1.2, 'g-.'); hold on;
    end
end



% add servo also in ensemble

if 1    
    thicknessServo = 41e-3;
    widthServo = 21e-3;
    
    servo_xy = [
        0   0
        thicknessServo   0
        thicknessServo   widthServo
        0   widthServo
        0   0
    ]';

    angle = deg2rad(10);
    delta = [12e-3 - xShift * chord; 0e-3];

    for i = 1 : size(servo_xy, 2)
        servo_xy(:, i) = R(angle) * servo_xy(:, i) + delta;
    end

    if isEnablePlotting
        subplot(2, 2, 1);
        for i = 1 : size(servo_xy, 2)
            plot(servo_xy(1, :), servo_xy(2, :), 'k-', 'LineWidth', 2);
        end
    end
end



%% write airfoil coordinates to file

if ~isCallingFromOutside || 1
    csvwrite('Main.csv', main_xy');
    csvwrite('Flap.csv', flap_xy');
    if hasSlat, csvwrite('Slat.csv', slat_xy'); end
end



%% display cost

if ~isCallingFromOutside
    cost
end
