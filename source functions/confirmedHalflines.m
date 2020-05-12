function [F, y] = confirmedHalflines(lineStats, theta, ipSolvedImage, Rver, Image)

% Code for half line confirmation after detection procedure. Please refer
% to the reference paper p.1670, last part of Section 3. 
%
% INPUTS
%       ** lineStats: 5x4 matrix which includes detection statistics for 
%       each wake.      
%
%       ** theta: Angles for Radon transform.
% 
%       ** ipSolvedImage: GMC reconstructed SAR image.
%
%       ** Rver: GMC reconstructed image.
%
%       ** Image: Ship-centred-masked SAR Image. Please refer to the 
%       reference document [1].
%
%
% OUTPUTS
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
[h, w] = size(Rver);
indI = lineStats(1, 2);
indJ = lineStats(1, 3);
linesmatV = zeros(h,w);
linesmatV((indJ-1)*h + indI) = 1;
linesmatV = iradon(linesmatV, theta);
if size(linesmatV,1)>length(ipSolvedImage)
    linesmatV = linesmatV(2:end-1,2:end-1);
end
lineSize2 = 1;
lineSize = 1;
linesmatV = imdilate(linesmatV/max(linesmatV(:)),strel('disk',lineSize));
linesmatV = imresize(linesmatV, size(ipSolvedImage));
linesmatV = max(0, linesmatV./(max(linesmatV(:))));
idx = not(linesmatV == 0);
half = length(ipSolvedImage)/2;
idxUp = [idx(1:half, :); zeros(half, 2*half)];
idxDown = [zeros(half, 2*half); idx(half+1:end, :)];
upImage = ipSolvedImage(idxUp == 1);
downImage = ipSolvedImage(idxDown == 1);
meanUp = mean(upImage);
meanDown = mean(downImage);
if meanDown < meanUp
    TurbulentHalfLine = idxDown;
    State = 1;
else
    TurbulentHalfLine = idxUp;
    State = 0;
end

linesmatV2 = zeros(h,w); % 1st narrow V wake
linesmatV3 = zeros(h,w); % 2nd narrow V wake
linesmatV4 = zeros(h,w); % 1st Kelvin Wake
linesmatV5 = zeros(h,w); % 2nd Kelvin wake
missingLines = lineStats(2:5, 1) == 0;

if missingLines(1) && missingLines(2)
    linesmatV2 = zeros(size(ipSolvedImage));
    linesmatV3 = zeros(size(ipSolvedImage));
elseif missingLines(1)
    indiVVV = lineStats(3, 2);
    indjVVV = lineStats(3, 3);
    linesmatV3((indjVVV-1)*h + indiVVV) = 1;
    linesmatV3 = iradon(linesmatV3, theta);
    if size(linesmatV3,1)>length(ipSolvedImage)
        linesmatV3 = linesmatV3(2:end-1,2:end-1);
    end
    linesmatV3 = imdilate(linesmatV3/max(linesmatV3(:)),strel('disk',lineSize));
    linesmatV3 = imresize(linesmatV3, size(ipSolvedImage));
elseif missingLines(2)
    indiVVV = lineStats(2, 2);
    indjVVV = lineStats(2, 3);
    linesmatV2((indjVVV-1)*h + indiVVV) = 1;
    linesmatV2 = iradon(linesmatV2, theta);
    if size(linesmatV2,1)>length(ipSolvedImage)
        linesmatV2 = linesmatV2(2:end-1,2:end-1);
    end
    linesmatV2 = imdilate(linesmatV2/max(linesmatV2(:)),strel('disk',lineSize));
    linesmatV2 = imresize(linesmatV2, size(ipSolvedImage));
else
    indiVVV = lineStats(2, 2);
    indjVVV = lineStats(2, 3);
    linesmatV2((indjVVV-1)*h + indiVVV) = 1;
    linesmatV2 = iradon(linesmatV2, theta);
    if size(linesmatV2,1)>length(ipSolvedImage)
        linesmatV2 = linesmatV2(2:end-1,2:end-1);
    end
    linesmatV2 = imdilate(linesmatV2/max(linesmatV2(:)),strel('disk',lineSize));
    linesmatV2 = imresize(linesmatV2, size(ipSolvedImage));
    
    indiVVV = lineStats(3, 2);
    indjVVV = lineStats(3, 3);
    linesmatV3((indjVVV-1)*h + indiVVV) = 1;
    linesmatV3 = iradon(linesmatV3, theta);
    if size(linesmatV3,1)>length(ipSolvedImage)
        linesmatV3 = linesmatV3(2:end-1,2:end-1);
    end
    linesmatV3 = imdilate(linesmatV3/max(linesmatV3(:)),strel('disk',lineSize));
    linesmatV3 = imresize(linesmatV3, size(ipSolvedImage));
end

if missingLines(3) && missingLines(4)
    linesmatV4 = zeros(size(ipSolvedImage));
    linesmatV5 = zeros(size(ipSolvedImage));
elseif missingLines(3)
    indiVVV = lineStats(5, 2);
    indjVVV = lineStats(5, 3);
    linesmatV5((indjVVV-1)*h + indiVVV) = 1;
    linesmatV5 = iradon(linesmatV5, theta);
    if size(linesmatV5,1)>length(ipSolvedImage)
        linesmatV5 = linesmatV5(2:end-1,2:end-1);
    end
    linesmatV5 = imdilate(linesmatV5/max(linesmatV5(:)),strel('disk',lineSize));
    linesmatV5 = imresize(linesmatV5, size(ipSolvedImage));
elseif missingLines(4)
    indiVVV = lineStats(4, 2);
    indjVVV = lineStats(4, 3);
    linesmatV4((indjVVV-1)*h + indiVVV) = 1;
    linesmatV4 = iradon(linesmatV4, theta);
    if size(linesmatV4,1)>length(ipSolvedImage)
        linesmatV4 = linesmatV4(2:end-1,2:end-1);
    end
    linesmatV4 = imdilate(linesmatV4/max(linesmatV4(:)),strel('disk',lineSize));
    linesmatV4 = imresize(linesmatV4, size(ipSolvedImage));
else
    indiVVV = lineStats(4, 2);
    indjVVV = lineStats(4, 3);
    linesmatV4((indjVVV-1)*h + indiVVV) = 1;
    linesmatV4 = iradon(linesmatV4, theta);
    if size(linesmatV4,1)>length(ipSolvedImage)
        linesmatV4 = linesmatV4(2:end-1,2:end-1);
    end
    linesmatV4 = imdilate(linesmatV4/max(linesmatV4(:)),strel('disk',lineSize));
    linesmatV4 = imresize(linesmatV4, size(ipSolvedImage));
    
    indiVVV = lineStats(5, 2);
    indjVVV = lineStats(5, 3);
    linesmatV5((indjVVV-1)*h + indiVVV) = 1;
    linesmatV5 = iradon(linesmatV5, theta);
    if size(linesmatV5,1)>length(ipSolvedImage)
        linesmatV5 = linesmatV5(2:end-1,2:end-1);
    end
    linesmatV5 = imdilate(linesmatV5/max(linesmatV5(:)),strel('disk',lineSize));
    linesmatV5 = imresize(linesmatV5, size(ipSolvedImage));
end
linesmatV2 = max(0, linesmatV2./(max(linesmatV2(:))));
linesmatV3 = max(0, linesmatV3./(max(linesmatV3(:))));
linesmatV4 = max(0, linesmatV4./(max(linesmatV4(:))));
linesmatV5 = max(0, linesmatV5./(max(linesmatV5(:))));
if State == 1
    linesmatV2(1:half, :) = 0;
    linesmatV3(1:half, :) = 0;
    linesmatV4(1:half, :) = 0;
    linesmatV5(1:half, :) = 0;
    NarrowV1HalfLine = linesmatV2;
    NarrowV2HalfLine = linesmatV3;
    KelvinHalfLine1 = linesmatV4;
    KelvinHalfLine2 = linesmatV5;
else
    linesmatV2(half+1:end, :) = 0;
    linesmatV3(half+1:end, :) = 0;
    linesmatV4(half+1:end, :) = 0;
    linesmatV5(half+1:end, :) = 0;
    NarrowV1HalfLine = linesmatV2;
    NarrowV2HalfLine = linesmatV3;
    KelvinHalfLine1 = linesmatV4;
    KelvinHalfLine2 = linesmatV5;
end

turHalfIntensity = ipSolvedImage(not(TurbulentHalfLine == 0));
nar1HalfIntensity = ipSolvedImage(not(NarrowV1HalfLine == 0));
nar2HalfIntensity = ipSolvedImage(not(NarrowV2HalfLine == 0));
kel1HalfIntensity = ipSolvedImage(not(KelvinHalfLine1 == 0));
kel2HalfIntensity = ipSolvedImage(not(KelvinHalfLine2 == 0));

if not(sum(nar1HalfIntensity)) == 0
    Val = 1;
    FI2 = mean(nar1HalfIntensity(nar1HalfIntensity<Val))/mean(ipSolvedImage(:)) - 1;
else
    FI2 = -1;
end
if not(sum(nar2HalfIntensity)) == 0
    Val = 1;
    FI3 = mean(nar2HalfIntensity(nar2HalfIntensity<Val))/mean(ipSolvedImage(:)) - 1;
else
    FI3 = -1;
end

if not(sum(kel1HalfIntensity)) == 0
    Val = 1;
    FI4 = mean(kel1HalfIntensity(kel1HalfIntensity<Val))/mean(ipSolvedImage(:)) - 1;
else
    FI4 = -1;
end

if not(sum(kel2HalfIntensity)) == 0
    Val = 1;
    FI5 = mean(kel2HalfIntensity(kel2HalfIntensity<Val))/mean(ipSolvedImage(:)) - 1;
else
    FI5 = -1;
end
FI1 = mean(turHalfIntensity)/mean(ipSolvedImage(:)) - 1;

y = Image;
rangecol = 1+(0:length(Image)-1);
y((1:length(Image)),rangecol) = Image;
y = repmat(y,[1 1 3]);
T = zeros(5, 1);
if FI1 < 0 
    y((1:length(Image)),rangecol,2) = y((1:length(Image)),rangecol,2) + imdilate(TurbulentHalfLine,strel('disk',lineSize2));
    y((1:length(Image)),rangecol,1) = y((1:length(Image)),rangecol,1) + imdilate(TurbulentHalfLine,strel('disk',lineSize2));
    T(1) = 1;
end
if FI2 > 0
    y((1:length(Image)),rangecol,2) = y((1:length(Image)),rangecol,2) + imdilate(NarrowV1HalfLine,strel('disk',lineSize2));
    T(2) = 1;
end
if FI3 > 0.1
    y((1:length(Image)),rangecol,2) = y((1:length(Image)),rangecol,2) + imdilate(NarrowV2HalfLine,strel('disk',lineSize2));
    T(3) = 1;
end
if FI4 > 0.1
    y((1:length(Image)),rangecol,1) = y((1:length(Image)),rangecol,1) + imdilate(KelvinHalfLine1,strel('disk',lineSize2));
    T(4) = 1;
end
if FI5 > 0.1
    y((1:length(Image)),rangecol,1) = y((1:length(Image)),rangecol,1) + imdilate(KelvinHalfLine2,strel('disk',lineSize2));
    T(5) = 1;
end
F = [FI1;FI2;FI3;FI4;FI5];
F = [F T];