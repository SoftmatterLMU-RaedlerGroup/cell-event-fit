function varargout = postproc_LATE_decay(D, F, R)
% This function contains the postprocessing routine for Calcium marker.
%
% Required arguments:
% ===================
% D					structure, containing these fields:
% D.data					array of measured data
% D.index					index of current trace in datafile
% R					structure, containing these fields:
% R.data_sim				array of simulated data
% R.max_val					maximum value
% R.max_ind					index of maximum
% R.min_val					minimum value
% R.min_ind					index of minimum
% R.amplitude				amplitude (R.max_val - R.min_val)
% R.amp_max					file-wide maximum value of amplitude
% R.params					vector of fit parameters
% F					structure, containing these fields:
% F.t_sim					time vector of simulated data (in hours)
% F.fullpath				absolute path to the initial file

%
% Return value:
% =============
% varargout{1}			onset time of the signal
% varargout{2}			structure with model-specific information
% varargout{3}			derivative of the rising edge
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

% Allocate memory
isGoodTrace = true;
t_onset = NaN;

% Define threshold values
assume_const_delta = 1E-2;
max_half_width = 4;
min_relative_peak_height = .1;
onset_slope_multiplicator = 1 / 20;
onset_height_relative_threshold = .05;

%% Calculate derivative
deriv = zeros(1, size(R.data_sim,1));
deriv(1) = (R.data_sim(2) - R.data_sim(1)) / (F.t_sim(2) - F.t_sim(1));
deriv(end) = (R.data_sim(end) - R.data_sim(end-1)) / (F.t_sim(end) - F.t_sim(end-1));
for i = 2:length(deriv)-1
	deriv(i) = (R.data_sim(i+1) - R.data_sim(i-1)) / (F.t_sim(i+1) - F.t_sim(i-1));
end

%% Find peak and minima
% Start at time 0 and climb down the decay
idx_min_before = 1;
while true
	if idx_min_before == length(R.data_sim)
		% Whole trace is monotonously decreasing: no onset found
		t_onset = -Inf;
		isGoodTrace = false;
		break
	elseif R.data_sim(idx_min_before + 1) > R.data_sim(idx_min_before)
		% Fluorescence starts increasing: minimum found
		break
	else
		% Climb one step further down
		idx_min_before = idx_min_before + 1;
	end
end

if isGoodTrace
	% Get positions and amplitudes related to peak
	[height_peak, idx_peak] = max(R.data_sim(idx_min_before:end));
	idx_peak = idx_peak + idx_min_before - 1;
	[min_after, idx_min_after] = min(R.data_sim(idx_peak:end));
	idx_min_after = idx_min_after + idx_peak - 1;
	min_before = R.data_sim(idx_min_before);
	amp_peak = height_peak - min_before;
	amp_decay = R.data_sim(1) - R.min_val;

	% Avoid indexing errors upon bad peak/minimum positions
	if idx_peak <= 1
		isGoodTrace = false;
	elseif idx_min_before <= 1
		isGoodTrace = false;
	end

	% DEBUG
	%fprintf('trace %d: min_before = (%fh | %f) (%d)\n', D.index, F.t_sim(idx_min_before), min_before, isGoodTrace)
end

if isGoodTrace
	%% Calculate FWHM
	% Half maximum before peak
	mid_amp_before = min_before + 0.5 * amp_peak;
	mid_amp_before_idx = idx_peak - 1;
	if mid_amp_before_idx > idx_min_before
		while R.data_sim(mid_amp_before_idx - 1) >= mid_amp_before
			mid_amp_before_idx = mid_amp_before_idx - 1;
% 			if mid_amp_before_idx <= 1
% 				fprintf('index %d: mid_amp_before_idx = %d', D.index, mid_amp_before_idx) % DEBUG
% 			end
			if mid_amp_before_idx == idx_min_before
				break
			end
		end
	end
	if abs(R.data_sim(mid_amp_before_idx - 1) - mid_amp_before) < ...
			abs(R.data_sim(mid_amp_before_idx) - mid_amp_before)
		mid_amp_before_idx = mid_amp_before_idx - 1;
	end
	%mid_amp_before = R.data_sim(mid_amp_before_idx);

	% Half maximum after peak
	if idx_min_after > idx_peak
		mid_amp_after = 0.5 * (height_peak + min_after);
		mid_amp_after_idx = idx_peak + 1;
		if mid_amp_after_idx < idx_min_after
			while R.data_sim(mid_amp_after_idx + 1) >= mid_amp_after
				mid_amp_after_idx = mid_amp_after_idx + 1;
				if mid_amp_after_idx == idx_min_after
					break
				end
			end
		end
		if mid_amp_after_idx < idx_min_after && ...
				abs(R.data_sim(mid_amp_after_idx + 1) - mid_amp_after) < ...
				abs(R.data_sim(mid_amp_after_idx) - mid_amp_after)
			mid_amp_after_idx = mid_amp_after_idx + 1;
		end
		%mid_amp_after = R.data_sim(mid_amp_after_idx);
	else
		mid_amp_after_idx = idx_peak;
	end

	% Get FHWM
	if mid_amp_after_idx > idx_peak && mid_amp_after_idx < length(F.t_sim)
		fwhm = F.t_sim(mid_amp_after_idx) - F.t_sim(mid_amp_before_idx);
	else
		fwhm = Inf;
	end

	%% Derivative at rising edge
	rising_edge_xoff = F.t_sim(mid_amp_before_idx);
	rising_edge_yoff = R.data_sim(mid_amp_before_idx);
	rising_edge_scale = deriv(mid_amp_before_idx);
% 	rising_edge_scale = (R.data_sim(mid_amp_before_idx+1)-R.data_sim(mid_amp_before_idx-1)) / ...
% 		(F.t_sim(mid_amp_before_idx+1) - F.t_sim(mid_amp_before_idx-1));
end

if ~isGoodTrace
	rising_edge_xoff = NaN;
	rising_edge_yoff = NaN;
	rising_edge_scale = NaN;
end

%% Populate status matrix

% Calculate onset time
if isGoodTrace
	if amp_peak < min_relative_peak_height * amp_decay
		% Peak is too small in relation to bleaching decay
		t_onset = NaN;

	elseif min(deriv) >= -assume_const_delta || max(deriv) <= assume_const_delta
		% Whole trace is actually flat
		t_onset = NaN;

	elseif fwhm > max_half_width
		% Too wide peak, unexplained behaviour
		t_onset = Inf;

	else
% 		% Set t_onset when ascending edge reaches threshold
% 		thresh = onset_height_relative_threshold * amp_peak + min_before;
% 		idx_onset = idx_min_before;
% 		while R.data_sim(idx_onset + 1) < thresh
% 			idx_onset = idx_onset + 1;
% 		end
% 		t_onset = F.t_sim(idx_onset);

		% Set t_onset when slope reaches threshold
		onset_slope = onset_slope_multiplicator * (R.data_sim(1) - min_before);
		idx_onset = find(deriv(idx_min_before:mid_amp_before_idx) <= onset_slope, 1, 'last') + idx_min_before - 1;
		if length(idx_onset) < 1 || idx_onset == mid_amp_before_idx
			idx_onset = idx_min_before;
		end
		t_onset = F.t_sim(idx_onset);
	end
end

%% Return values

% Onset time
varargout{1} = t_onset;

% Fit type
if nargout > 1
	varargout{2}.type = int8(4); % type: min between peaks
	varargout{2}.rising_edge_xoff = rising_edge_yoff;
	varargout{2}.rising_edge_yoff = rising_edge_xoff;
	varargout{2}.rising_edge_scale = 1 / rising_edge_scale;
end

% Initial derivative
if nargout > 2
	if isfinite(t_onset)
		varargout{3} = 1 / rising_edge_scale;
	else
		varargout{3} = NaN;
	end
end
