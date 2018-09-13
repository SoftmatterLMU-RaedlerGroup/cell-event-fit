function varargout = postproc_kink(~, F, R)
% This function contains the postprocessing routine for
% LysoTracker and CellROX.
%
% Required arguments:
% ===================
% %D					structure, containing these fields:
% %D.data					array of measured data
% %D.index					index of current trace in datafile
% R					structure, containing these fields:
% R.data_sim				array of simulated data
% R.max_val					maximum value
% R.max_ind					index of maximum
% R.min_val					minimum value
% R.min_ind					index of minimum
% R.amplitude				amplitude (R.max_val - R.min_val)
% R.params					vector of fit parameters
% F					structure, containing these fields:
% F.t_sim					time vector of simulated data (in hours)
% F.modelname				name of the model used
%
% Return value:
% =============
% varargout{1}			breakdown time of the signal
% varargout{2}			structure with model-specific information
% varargout{3}			maximum derivative of the rising edge
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

%% Define thresholds
assume_const_delta = 1E-2 / 6;
min_half_width = 1.7;
min_breakdown_amp = .1 * R.amplitude;
min_breakdown_amp_smooth = .2 * R.amplitude;
parab_threshold = 0.03;
min_decay_steepness = 0.07;

%% Set default values
parab_threshold = parab_threshold * R.amplitude;
%deriv_min_ind = NaN;
t_breakdown = NaN;
breakdown_too_smooth = false;
parabola = NaN(length(F.t_sim), 1);


%% Get parabola
% Obtain parabola coefficients (index depends on model)
switch F.modelname

	case 'ROS_Parabola'
		vertex_pos = R.params(3);
		coeff_const = R.params(1);
		coeff_lin = 0;
		coeff_quad = - R.params(2);
		decay_steepness = R.params(5);

	case 'Lyso'
		vertex_pos = R.params(3);
		coeff_const = R.params(1);
		coeff_lin = 0;
		coeff_quad = - R.params(2);
		decay_steepness = R.params(5);

	case 'Lyso1'
		vertex_pos = R.params(3);
		coeff_const = R.params(4);
		coeff_lin = R.params(5);
		coeff_quad = - R.params(6);
		decay_steepness = R.params(8);

	case 'Lyso2'
		vertex_pos = R.params(1);
		coeff_const = R.params(2);
		coeff_lin = R.params(3);
		coeff_quad = - R.params(4);
		decay_steepness = R.params(6);

	case 'Lyso3'
		vertex_pos = R.params(1);
		coeff_const = R.params(2);
		coeff_lin = 0;
		coeff_quad = - R.params(3);
		decay_steepness = R.params(5);

	otherwise
		error(['Unknown model: ' F.modelname]);
end

% Calculate parabola
for i = 1:length(F.t_sim)
	x = F.t_sim(i) - vertex_pos;
	parabola(i) = coeff_const + x * ( coeff_lin + x * coeff_quad );

	% Prevent display of irrelevant values of parabola (design issue)
	if parabola(i) < R.min_val - 0.1 * R.amplitude
		parabola(i) = -Inf;
	elseif parabola(i) > R.max_val + 0.5 * R.amplitude
		parabola(i) = Inf;
	end
end

%% Get 1st derivative of flourescence signal w.r.t. time
deriv = NaN(size(R.data_sim,1), 1);

for i = 2:length(R.data_sim)
	deriv(i) = ( R.data_sim(i) - R.data_sim(i-1) ) ...
		/ ( F.t_sim(i) - F.t_sim(i-1) );
end

% Get maximum of derivative of rising edge
rising_deriv = deriv(2);
rising_ind = 2;

for i = 3:R.max_ind
	if deriv(i) > rising_deriv
		rising_deriv = deriv(i);
		rising_ind = i;
	elseif deriv(i) < 0
		break
	end
end

% Find minimum of derivative of falling edge
[~,deriv_min_ind] = min(deriv(R.max_ind+1:end));
deriv_min_ind = deriv_min_ind + R.max_ind;

% Define rising/falling edge coefficients
rising_edge_xoff = R.data_sim(rising_ind);
rising_edge_yoff = F.t_sim(rising_ind);
rising_edge_scale = 1 / deriv(rising_ind);

if ~isempty(deriv_min_ind)
	falling_edge_xoff = R.data_sim(deriv_min_ind);
	falling_edge_yoff = F.t_sim(deriv_min_ind);
	falling_edge_scale = 1 / deriv(deriv_min_ind);
else
	falling_edge_xoff = NaN;
	falling_edge_yoff = NaN;
	falling_edge_scale = NaN;
end

%% Calculate full width at half maximum
middle_amplitude = R.min_val + R.amplitude / 2;
%middle_amplitude = model.simulate(R.params(3),R.params);

% Find index next to middle amplitude in rising edge
min_diff = Inf;
for i = R.max_ind:-1:1
	temp = abs(R.data_sim(i) - middle_amplitude);

	if temp < min_diff
		mid_amp_ind1 = i;
		min_diff = temp;
	elseif isfinite(min_diff) && temp > min_diff
		break
	end
end

% Find index next to middle amplitude in falling edge
min_diff = Inf;
mid_amp_ind2 = Inf;
for i = R.max_ind:length(R.data_sim)
	temp = abs(R.data_sim(i) - middle_amplitude);

	if temp < 0.3 * R.amplitude && temp < min_diff
		mid_amp_ind2 = i;
		min_diff = temp;
	elseif isfinite(min_diff) && temp > min_diff
		break
	end
end

% Get full width at half maximum
if ~isfinite(mid_amp_ind2)
	fwhm = Inf;
else
	fwhm = F.t_sim(mid_amp_ind2) - F.t_sim(mid_amp_ind1);
end

%% Sort out badly-behaved cells; set breakdown time
% Sort out cells dying too early or too late
if max(deriv) <= assume_const_delta || R.max_ind == 1
	% Cell died before measurement
	t_breakdown = -Inf;

elseif R.data_sim(end) == R.max_val || deriv_min_ind == length(R.data_sim) ...
		|| min(deriv) >= assume_const_delta
	% Cell died to late
	t_breakdown = Inf;

elseif fwhm < min_half_width
	% Peak not wide enough
	t_breakdown = NaN;

elseif decay_steepness < min_decay_steepness
	% Sigmoid decay not sharp enough
	breakdown_too_smooth = true;

else
	% Cell dies within observed time

	% Locate kink by deviation from parabola
	for i = 1:length(F.t_sim)
		parab_devi = parabola(i) - R.data_sim(i);

		if abs(parab_devi) <= parab_threshold
			% Nothing interesting happens
			continue

		elseif parab_devi < 0
			% Fit is larger than parabola
			breakdown_too_smooth = true;
			break

		elseif i == 1 && (parabola(i) - R.data_sim(i) > 0)
			% Fit deviates from parabola already at the first time
			breakdown_too_smooth = true;
			break

		elseif i > 1 && i ~= length(R.data_sim)
			% Breakdown time found
			t_breakdown = F.t_sim(i);
			breakdown_too_smooth = false;

			% Ensure that we have no too strong ascent any more
			if R.data_sim(i) - R.data_sim(end) < min_breakdown_amp
				% Too strong ascent after breakdown
				t_breakdown = -Inf;

			elseif i <= R.max_ind
				% Breakdown was detected before maximum
				% Go on to extrapolation
				breakdown_too_smooth = true;
			end

			break
		else
			continue
		end
	end
end

if breakdown_too_smooth
	%% Sigmoid function not sharp enough
	% Use extrapolation method instead

	% Test if breakdown amplitude is large enough
	if R.max_val - min(R.data_sim(R.max_ind:end)) > min_breakdown_amp_smooth

		% Calculate breakdown time by extrapolation
		t_breakdown = ( R.max_val - falling_edge_xoff ) ...
			* falling_edge_scale + falling_edge_yoff;
	end	
end

%% Return values

% Breakdown time
varargout{1} = t_breakdown;

% Fit type
if nargout > 1
	varargout{2} = struct( ...
		'parabola', parabola, ...
		'rising_edge_xoff', rising_edge_xoff, ...
		'rising_edge_yoff', rising_edge_yoff, ...
		'rising_edge_scale', rising_edge_scale, ...
		'falling_edge_xoff', falling_edge_xoff, ...
		'falling_edge_yoff', falling_edge_yoff, ...
		'falling_edge_scale', falling_edge_scale ...
		);

	if ~breakdown_too_smooth
		% Used parabola kink algorithm
		varargout{2}.type = int8(2); % type: kink_parabola

	else
		% Used alternative algorithm
		varargout{2}.type = int8(-2); % type: nokink_parabola
	end
end

% Linear parameter of parabola in normal form
if nargout > 2
	varargout{3} = rising_deriv;
end
