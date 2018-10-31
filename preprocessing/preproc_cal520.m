function data_index = preproc_cal520(data)
%preproc_cal520 contains the preprocessing routine for the Cal520 signal.
%
% The signal may contain an outlier (low value) at the first point.
% Exclude such a point from fitting.
%
% Required variables:
% ===================
% data				vector of measured data
% 
% Provided variables:
% ===================
% data_index		indices of datapoints to be used for fitting
%
% HINTS:
% 1.	Due to parallelization errors, "data" must only be the vector of
%		the current trace, and not an array of all traces.
% 2.	If an empty vector is returned, it will automatically be
%		populated with all indices of the whole track.
%
% Copyright © 2018 Daniel Woschée <daniel.woschee@physik.lmu.de>
% Faculty of Physics / Ludwig-Maximilians-Universität München
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

% Set threshold for excluding first point
thresh = 5;

% Get difference between first and second point
init_diff = data(2) - data(1);

% Get difference between 
consec_diff = abs(data(3) - data(2));

% Compare to threshold and return corresponding index array
if init_diff > thresh * consec_diff
	data_index = 2:length(data);
else
	data_index = [];
end
