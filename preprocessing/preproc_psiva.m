function data_index = preproc_psiva(data, ~, ~)
%preproc_psiva contains the preprocessing routine for the pSIVA signal.
%
% For pSIVA, only fit to the minimum after the first peak if
% there are two large peaks
%
% Required variables:
% ===================
% data				vector of measured data
% 
% Provided variables:
% ===================
% data_index		indices of datapoints to be used for fitting
%
% The goal of this preprocessing is to only fit up to the minimum
% after the first peak (if there are several peaks).
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

% Get the amplitude of the signal
amplitude = max(data) - min(data);

% Define minimum prominence (peaks below this level will be ignored)
min_peak_prom = amplitude * 0.7;

% Get indices and prominences of the two most prominent peaks
% (actually we get all peaks, but we only need the two largest)
[~,peak_loc,peak_wid,peak_prom] = findpeaks(data, ...
							'SortStr', 'descend', ...
							'WidthReference', 'halfprom');

if peak_prom(1) >= min_peak_prom && ...
		( peak_prom(2) >= .15 * peak_prom(1) || peak_wid(2) >= .15 * peak_wid(1) )
	% There are two large peaks => find the minimum in between
	if peak_loc(1) < peak_loc(2)
		[~,latest_t] = min(data( peak_loc(1) : peak_loc(2)));
	else
		[~,latest_t] = min(data( peak_loc(2) : peak_loc(1)));
	end

	% Fit the data only until to the minimum after first peak
	data_index = 1:latest_t;
else
	data_index = [];
end