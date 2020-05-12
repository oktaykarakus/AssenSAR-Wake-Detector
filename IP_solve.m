% Code for inverse problem solution to promote linear structures for ship 
% wake detection in SAR images in terms of the GMC penalty function.
%
% Some Important Variables
%       ** thetaRotate: Inverse problem will be solved for thetaRotate 
%       angle interval in Radon domain. The full range is 180 degrees.
%       However for computational reasons, the range between -90 to 90 
%       degrees into "numRotate" sections. Default usage is 45.
%
%       ** GMCImage: Stores reconstructed SAR images for each theta
%       interval.
%
%       ** lambda: Regularisation constant for the GMC penalty function.
%
%       ** gamma: Convexity parammeter for the GMC penalty.
%
%       ** Image: Ship-centred-masked SAR Image. Please refer to the 
%       reference document [1].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
addpath('.\source functions');
addpath('.\saved data');
addpath('.\Images');
%% Parameter Initialization
thetaRotate = 45;
numRotate = 180/thetaRotate;               
thetaInt = -90:0.25:(-90+thetaRotate);       
normType = 'GMC';                                                           % Norm type. For storing only.
lambda = 50;                                                                
gamma = 0.9;                                                                                     
load('testImage.mat');                                                      % Loading the test image. It returns an image within a variable "Image"
%% Inverse problem solution
GMCImage = zeros([size(Image) numRotate]);
for rotation = 1:numRotate
    fprintf('Processing... %d/%d\n', rotation, numRotate)
    theta = thetaInt + thetaRotate*(rotation-1);
    A  = @(x) radonT(x,theta);                                              % Inverse Radon Transform operator - C
    AH = @(x) radon(x,theta);                                               % Radon Transform Operator - R
    temp = GMC_regularisation(Image, A, AH, 1, lambda, gamma);                     % Variable "temp" is the solution in Radon space
    image_GMC = A(temp);                                                    % Tranform into image domain.
    image_GMC = image_GMC/max(image_GMC(:));                                % Normalization between 0 and 1.
    GMCImage(:, :, rotation) = image_GMC;                                   % Storing results.
end
%% Storing processed data for future use
cd('.\saved data')
filename = 'wakeDetection_test.mat';
save(filename)
cd('.\..')