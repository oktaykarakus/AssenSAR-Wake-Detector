function Ax = radonT(x, theta)

% Code for inverse Radon transform. It includes removing the border pixels.
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
%
% inverse radon (line -> image)
Ax = iradon(x,theta);%,'Hamming'); % filtering if want super smooth
Ax = max(0,Ax(2:end-1,2:end-1));