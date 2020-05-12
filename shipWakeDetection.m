% Code for detecting ship wakes by using the inverse problem solution via
% the code "IP_solve.m"
% We  assume that 5 ship wake structures, which are:
% 1 Turbulent wake,
% 2 arms of Narrow V wake,
% 2 arms of Kelvin wake.
% For details please refer to the reference given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some Important Variables
%       ** thetaRotate: refer to "IP_solve.m.
%
%       ** GMCImage: refer to "IP_solve.m.
%
%       ** SinMax: Maximum amplitude (A) for sine wave. Reference paper [1]
%       p.1669 --> A.sin(theta).
%
%       ** minKelvinAngleConstant: Minimum angle for Kelvin wake search 
%       range. 
%
%       ** maxKelvinAngleConstant: Maximum angle for Kelvin wake search 
%       range. 
%
%       ** indiV & indjV: Stores regional maximum points in Radon space.
%
%       ** indiVV & indjVV: Stores regional minimum points in Radon space.
%
%       ** minMaxPair: Stores possible detection results for the
%       turbulent-1st Narrow V wake pair. 
%
%       ** lineStats: Stores detection statistics for each wake given
%       above. Each column refers to: 
%           Detection Status(1 or 0) & Row index of detected point & Column
%           index for detected point & Amplitude difference
%
%       ** F: Measure index for each variable. Please refer to (35) in the
%       reference paper [1] p.1670.
%
%       ** y: Resuting image with embedded detected-validated ship wakes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICENSE
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% Copyright © Oktay Karakus,PhD 
% o.karakus@bristol.ac.uk
% University of Bristol, UK
% October, 2018
% Updated: April 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE
%
% [1] O Karakus, I Rizaev, and A Achim. "Ship Wake Detection in SAR Images 
%       via Sparse Regularisation."
%       IEEE Transactions on Geoscience and Remote Sensing, 2020.
% [2] O Karakus, and A Achim. "Ship Wake Detection in X-band SAR Images 
%       Using Sparse GMC Regularization." 
%       IEEE ICASSP 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
clc
addpath('.\source functions');
addpath('.\saved data');
addpath('.\Images');
%
load('wakeDetection_test.mat')          % Loading saved results obtained from the code 'IP_solve.m'
%% Obtaining inverse problem solution from saved data 
thetaInt = -90:0.25:(-90+thetaRotate); 
SinMax = round(length(Image)/10);                                   % Maximum amplitude (A) for sine wave. arXiv paper page p.7 A.sin(theta)
verLength = round(SinMax);
minKelvinAngleConstant = 10;                                        % Kelvin wake search minimum angle
maxKelvinAngleConstant = 20;                                        % Kelvin wake search maximum angle
% Obtain Radon transform of observed SAR image
thetaInt2 = -90:0.25:90;
[radonImage, xp] = radon(Image, thetaInt2);
radonImage = radonImage./max(radonImage(:));
% Start merging inv.problem solutions for each 45 degrees intervals
radonMergedImage = [];
radonMergedImageNormalized = [];
numRotate = 180/thetaRotate;
for rotation = 1:numRotate
    thetaInt2 = thetaInt + thetaRotate*(rotation-1);
    theta = thetaInt2;
    A  = @(x) radonT(x,theta);
    AH = @(x) radon(x,theta);
    image_GMC = GMCImage(:, :, rotation);
    [Rver, xp] = AH(image_GMC);%
    maxu = max(max(Rver(20:end-20,20:end-20))); %                   % Normalize Radon Image ignoring side effects
    Rver = min(1,max(0,Rver/maxu));
    radonMergedImage = [radonMergedImage Rver];
end
% Remove overlapping samples
if numRotate == 4
    Rver = [radonMergedImage(:, 1:182) radonMergedImage(:, 184:363) radonMergedImage(:, 365:544) radonMergedImage(:, 546:end)];
elseif numRotate == 2
    Rver = [radonMergedImage(:, 1:363) radonMergedImage(:, 365:end)];
elseif numRotate == 1
    Rver = radonMergedImage;
end
thetaInt2 = -90:0.25:90;
A  = @(x) radonT(x, thetaInt2);                                     % Redefine operator A for the whole angle range
ipSolvedImage = A(Rver);                                            % Obtain inverse problem solved SAR image by taking inv.Radon of result Rver.
ipSolvedImage = ipSolvedImage./max(ipSolvedImage(:));
[h, ~] = size(Rver);

%% Wake detection starts here
% Some filtering and blurring operations to make peaks clearer in Radon space
Len = round(find(xp == 0));
Len2 = 1/mean(diff(thetaInt2));
blurRF = imfilter(Rver./max(Rver(:)),fspecial('gaussian',3,3));%
blurRF = imdilate(blurRF,strel('disk',3));
% Detect regional maximum points in Radon space
localmx = imregionalmax(blurRF);
stats = regionprops(localmx,'Centroid');
indiV = zeros(length(stats),1);
indjV = zeros(length(stats),1);
for k = 1:length(stats)
    indiV(k) = round(stats(k).Centroid(2));
    indjV(k) = round(stats(k).Centroid(1));
end
% Remove points lie outside of two sine waves
idxI = or(xp(indiV)' < -cosd(thetaInt2(indjV))*SinMax, xp(indiV)' > cosd(thetaInt2(indjV))*SinMax)';
idxJ = or(indjV < Len2+1, indjV > size(Rver, 2)-Len2);
indiV(or(idxI, idxJ)) = [];
indjV(or(idxI, idxJ)) = [];

% Some filtering and blurring operations to make minimum points clearer in Radon space
blurRF = imfilter(Rver./max(Rver(:)),fspecial('gaussian',3,3));%
blurRF = imdilate(blurRF,strel('disk',1));
% Detect regional minimum points in Radon space
localmn = imregionalmin(blurRF);
stats = regionprops(localmn,'Centroid');
indiV2 = zeros(length(stats),1);
indjV2 = zeros(length(stats),1);
for k = 1:length(stats)
    indiV2(k) = round(stats(k).Centroid(2));
    indjV2(k) = round(stats(k).Centroid(1));
end
% Remove points lie outside of two sine waves
idxI = or(xp(indiV2)' < -cosd(thetaInt2(indjV2))*SinMax, xp(indiV2)' > cosd(thetaInt2(indjV2))*SinMax)';
idxJ = or(indjV2 < Len2+1, indjV2 > size(Rver, 2)-Len2);
indiV2(or(idxI, idxJ)) = [];
indjV2(or(idxI, idxJ)) = [];
tempRver2 = Rver((indjV2-1)*h + indiV2);
%Combine min and max points into same Vector for x and y axes separetely
indiVV = [indiV; indiV2];
indjVV = [indjV; indjV2];
clear local*
%% Turbulent & 1st-Narrow-V-Wake pair detection
% Detect min/max pairs
minMaxPair = [];
tempMinMaxPair = [];
for i = 1:length(indiV2)
    jIdx = indjV2(i);
    len = verLength;
    len2 = 4/mean(diff(thetaInt2)); % maximum angle between a Narrow V-wake and the Turbulent wake
    jInterval = jIdx-len2+1:jIdx+len2;
    tempIdx = indjV >= jIdx-len2+1 & indjV <= jIdx+len2;
    if not(sum(tempIdx) == 0)
        tempIndiV = indiV(tempIdx);
        tempIndjV = indjV(tempIdx);
        distance = abs(tempIndjV - jIdx);
        valuesForMaxPoints = Rver((tempIndjV-1)*h + tempIndiV);
        maxPointValue = max(valuesForMaxPoints);
        maxPointIdxAll = find(valuesForMaxPoints == maxPointValue);
        for j = 1:length(maxPointIdxAll)
            maxPointIdx = maxPointIdxAll(j);
            minPointValue = Rver((indjV2(i)-1)*h + indiV2(i));
            differenceMinMax = maxPointValue - minPointValue;
            blurMaxValue = blurRF((tempIndjV(maxPointIdx)-1)*h + tempIndiV(maxPointIdx));
            blurMinValue = blurRF((indjV2(i)-1)*h + indiV2(i));
            blurDifferenceMinMax = blurMaxValue - blurMinValue;
            realMaxValue = radonImage((tempIndjV(maxPointIdx)-1)*h + tempIndiV(maxPointIdx));
            realMinValue = radonImage((indjV2(i)-1)*h + indiV2(i));
            realDifferenceMinMax = realMaxValue - realMinValue;
            if abs(tempIndiV(maxPointIdx) - indiV2(i)) < len
                maxPoint = [tempIndiV(maxPointIdx) tempIndjV(maxPointIdx)];
                minMaxPair = [minMaxPair; maxPoint maxPointValue indiV2(i) indjV2(i) minPointValue distance(maxPointIdx) differenceMinMax realDifferenceMinMax blurDifferenceMinMax];
            else
                tempMaxPoint = [tempIndiV(maxPointIdx) tempIndjV(maxPointIdx)];
                tempMinMaxPair = [tempMinMaxPair; tempMaxPoint maxPointValue indiV2(i) indjV2(i) minPointValue distance(maxPointIdx) differenceMinMax realDifferenceMinMax blurDifferenceMinMax];
            end
        end
    end
end
[h, w] = size(Rver);
if isempty(minMaxPair)
    error('Turbulent and 1st Narraw-V wake pair cannot be detected!')
else % if some pairs deteced, then choose the pair which has maximum difference 
    fprintf('Turbulent and 1st Narrow V Vake pair is detected!\n\n');
    maxIdx = find(minMaxPair(:, 8) >= max(minMaxPair(:, 8)));
    turWake2(1) = 1;
    narrowVWake2(1) = 1;
    % Take x and y axis values of Turbulent and 1st Narrow V-wake
    turWake2(2:4) = minMaxPair(maxIdx(1), 4:6); % Turbulent Wake
    narrowVWake2(2:4) = minMaxPair(maxIdx(1), 1:3); % 1st Narrow V-wake
    differenceVector2 = minMaxPair(maxIdx(1), 8:10);
end

lineStats = zeros(5, 3);
lineStats(1, 1:4) = turWake2; % Tubulent Wake
lineStats(2, 1:4) = narrowVWake2; % 1st Narrow V-wake
%% 2nd Narrow V-Vake detection
% Search other side of turbulent wake for the 2nd Narroe V-wake
distanceBtwTurV = lineStats(1, 3) - lineStats(2, 3);
if distanceBtwTurV == 0
    distanceBtwTurV = 1;
end
secondVVakeThetaInt = (lineStats(1, 3)):sign(distanceBtwTurV):(lineStats(1, 3) + sign(distanceBtwTurV)*len2);
startIdxJ = min(secondVVakeThetaInt);
endIdxJ = max(secondVVakeThetaInt);
inRangeMaxPointIdxJ = indjV >= startIdxJ & indjV <= endIdxJ;
if lineStats(1, 2) - lineStats(2, 2) >= 0
    startIdxI = lineStats(1, 2);
    endIdxI = round(lineStats(1, 2)+(len/1));
else
    startIdxI = round(lineStats(1, 2)-(len/1));
    endIdxI = lineStats(1, 2);
end
inRangeMaxPointIdxI = indiV >= startIdxI & indiV <= endIdxI;
possibleSecondVVakeIdx = inRangeMaxPointIdxJ & inRangeMaxPointIdxI;
if not(sum(possibleSecondVVakeIdx) == 0)
    fprintf('2nd Narrow V Vake is detected!\n\n');
    tempIndiV = indiV(possibleSecondVVakeIdx);
    tempIndjV = indjV(possibleSecondVVakeIdx);
    valuesForMaxPoints = Rver((tempIndjV-1)*h + tempIndiV);
    maxPointValue = max(valuesForMaxPoints); % If there is more than one detection, only take the one that has maximum value
    maxPointIdx = find(valuesForMaxPoints == maxPointValue);
    lineStats(3, 1) = 1;
    lineStats(3, 2) = tempIndiV(maxPointIdx(end));
    lineStats(3, 3) = tempIndjV(maxPointIdx(end));
    lineStats(3, 4) = maxPointValue;
else
    fprintf('2nd Narrow V Vake cannot be detected!\n\n');
end
len = len*round(SinMax/verLength);
%% Turbulent and Narrov V wakes removal from the candidate lines
% Before detecting Kelvin wakes, detected points before this line should be
% removed to prevent misdetection
tempIdxI = indiV == lineStats(2, 2);
tempIdxJ = indjV == lineStats(2, 3);
tempIdx3 = not(tempIdxI & tempIdxJ);
tempIdxI2 = indiV == lineStats(3, 2);
tempIdxJ2 = indjV == lineStats(3, 3);
tempIdx2 = not(tempIdxI2 & tempIdxJ2);
tempIdx = tempIdx3 & tempIdx2;
indiV_new = indiV(tempIdx);
indjV_new = indjV(tempIdx);

%% 1st Kelvin wake detection
% Same procedure in Narrow V-wake detection but for different angle
% distance from the Turbulent wake. This part starts searching from the
% right of the Turbulent wake as 1st Kelvin wake. It can be reversed. No 
%important change on results observed.
startIdxJ = lineStats(1, 3) + minKelvinAngleConstant*len2/4;
endIdxJ = min(length(thetaInt2), round(lineStats(1, 3)+(maxKelvinAngleConstant*len2/4)));
inRangeMaxPointIdxJ = indjV_new >= startIdxJ & indjV_new <= endIdxJ;

startIdxI = round(lineStats(1, 2)-(len));
endIdxI = round(lineStats(1, 2)+(len));
inRangeMaxPointIdxI = indiV_new >= startIdxI & indiV_new <= endIdxI;
possibleFirstKelvinVakeIdx = inRangeMaxPointIdxJ & inRangeMaxPointIdxI;
if not(sum(possibleFirstKelvinVakeIdx) == 0)
    tempIndiV = indiV_new(possibleFirstKelvinVakeIdx);
    tempIndjV = indjV_new(possibleFirstKelvinVakeIdx);
    fprintf('1st Kelvin Vake is detected!\n\n');
    if length(tempIndiV) > 1
        valuesForMaxPoints = Rver((tempIndjV-1)*h + tempIndiV);
        terter = radonImage((tempIndjV-1)*h + tempIndiV);
        maxPointValue = max(valuesForMaxPoints.*terter);
        maxPointIdx = find(valuesForMaxPoints.*terter == maxPointValue);
        lineStats(4, 1) = 1;
        lineStats(4, 2) = tempIndiV(maxPointIdx);
        lineStats(4, 3) = tempIndjV(maxPointIdx);
        lineStats(4, 4) = maxPointValue;
    else
        valuesForMaxPoints = Rver((tempIndjV-1)*h + tempIndiV);
        lineStats(4, 1) = 1;
        lineStats(4, 2) = tempIndiV;
        lineStats(4, 3) = tempIndjV;
        lineStats(4, 4) = valuesForMaxPoints;
    end
else
    fprintf('1st Kelvin Vake cannot be detected!\n\n');
end

%% 2nd Kelvin Wake detection
% Now search left side of the Turbulent wake.
startIdxJ = max(0, round(lineStats(1, 3)-(maxKelvinAngleConstant*len2/4)));
endIdxJ = lineStats(1, 3)-minKelvinAngleConstant*len2/4;
inRangeMaxPointIdxJ = indjV_new >= startIdxJ & indjV_new <= endIdxJ;

if lineStats(4, 2) ==0
    startIdxI = round(lineStats(1, 2)-(len));
    endIdxI = round(lineStats(1, 2)+(len));
else
    if lineStats(1, 2) - lineStats(4, 2) >= 0
        startIdxI = lineStats(1, 2);
        endIdxI = round(lineStats(1, 2)+(len/1));
    else
        startIdxI = round(lineStats(1, 2)-(len/1));
        endIdxI = lineStats(1, 2);
    end
end

inRangeMaxPointIdxI = indiV_new >= startIdxI & indiV_new <= endIdxI;
possibleSecondKelvinVakeIdx = inRangeMaxPointIdxJ & inRangeMaxPointIdxI;
if not(sum(possibleSecondKelvinVakeIdx) == 0)
    tempIndiV = indiV_new(possibleSecondKelvinVakeIdx);
    tempIndjV = indjV_new(possibleSecondKelvinVakeIdx);
    fprintf('2nd Kelvin Vake is detected!\n\n');
    if length(tempIndiV) > 1
        valuesForMaxPoints = Rver((tempIndjV-1)*h + tempIndiV);
        maxPointValue = max(valuesForMaxPoints);
        maxPointIdx = find(valuesForMaxPoints == maxPointValue);
        lineStats(5, 1) = 1;
        lineStats(5, 2) = tempIndiV(maxPointIdx(1));
        lineStats(5, 3) = tempIndjV(maxPointIdx(1));
        lineStats(5, 4) = maxPointValue;
    else
        valuesForMaxPoints = Rver((tempIndjV-1)*h + tempIndiV);
        lineStats(5, 1) = 1;
        lineStats(5, 2) = tempIndiV;
        lineStats(5, 3) = tempIndjV;
        lineStats(5, 4) = valuesForMaxPoints;
    end
else
    fprintf('2nd Kelvin Vake cannot be detected!\n\n');
end
% Up to this points, some of the wakes detected. After this line,
% confirmation of the detected wakes will be performed. Please read end of 
% section III in arXiv paper.
% The function confirmedHalflines() takes detected points in Radon space
% and perform removing half lines and confirmation of those halg lines. 
% Please read arXiv paper Section 3. 
[F, y] = confirmedHalflines(lineStats, thetaInt2, ipSolvedImage, Rver, Image);

% Plot detected points in Radon Domain
figure(2);
imagesc(thetaInt2, xp, Rver); colormap gray; hold on;
plot(thetaInt2, cosd(thetaInt2)*SinMax, 'LineWidth', 2)
plot(thetaInt2, -cosd(thetaInt2)*SinMax, 'LineWidth', 2)
plot(thetaInt2, (thetaInt2)*0, 'k:', 'LineWidth', 2)
shape = {'x', 'o', 'd', '>', '>'};
colour = {'y', 'g', 'g', 'r', 'r'};
for i = 1:5
    if lineStats(i, 1) == 1
        plot(thetaInt2(lineStats(i, 3)), xp(lineStats(i, 2)), [colour{i} shape{i}]);
    end
end
hold off

% Plot detected wakes on the image
figure(1)
imshowpair(Image, y, 'montage')
drawnow