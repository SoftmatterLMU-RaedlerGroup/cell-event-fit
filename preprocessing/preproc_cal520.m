function data_index = preproc_cal520(data, mf, trace_ind)
%preproc_cal520 contains the preprocessing routine for the Cal520 signal.
%
% The signal may contain an outlier (low value) at the first point.
% Exclude such a point from fitting.
%
% Required variables:
% ===================
% data				vector of measured data
% mf				matfile handle for further information
% trace_ind			index of trace in `mf`
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
% Copyright © 2018-2019 Daniel Woschée <daniel.woschee@physik.lmu.de>
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

file_ind = mf.file_ind(trace_ind, 1);
data_index = 1:mf.data_len(file_ind, 1);

% Check for measurement 20181130
% In measurement 20181130, a systematic transient dip occurs
% for all traces during time points 12 and 13. Exclude those points.
file_name = mf.name_F(mf.name_F_ind(file_ind, 1), 1:mf.name_F_len(file_ind, 1));
if regexpi(file_name, '^\d+_\d+_20181130$')
	data_index = [data_index(1:11) data_index(14:end)];
end

% Set threshold for excluding first point
thresh = 5;

% Get difference between first and second point
init_diff = data(data_index(2)) - data(data_index(1));

% Get difference between 
consec_diff = abs(data(data_index(3)) - data(data_index(2)));

% Compare to threshold and return corresponding index array
if init_diff > thresh * consec_diff
	data_index = data_index(2:end);
end
