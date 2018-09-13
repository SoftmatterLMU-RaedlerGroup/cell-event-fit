function varargout = postproc_LATE_extrapol(D, F, R)
% This function contains the postprocessing routine for PI, Caspase, and pSIVA.
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
ascent_time = @(x) NaN;
mid_amp_deriv = NaN;

% Define threshold values
assume_const_delta = 1E-2;
min_half_width = 2.5;

%% Calculate FWHM
middle_amplitude = R.min_val + R.amplitude / 2;
mid_amp_ind = NaN;
%middle_amplitude = model.simulate(R.params(3),R.params);

%% Get 1st derivative of flourescence signal w.r.t. time
deriv = zeros(1,size(R.data_sim,1));
deriv(1) = NaN;

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


% Define rising edge coefficients
rising_edge_xoff = R.data_sim(rising_ind);
rising_edge_yoff = F.t_sim(rising_ind);
rising_edge_scale = 1 / deriv(rising_ind);
%rising_edge = @(y) (y - R.data_sim(rising_ind)) / deriv(rising_ind) + F.t_sim(rising_ind);

% Find index next to middle amplitude in rising edge
min_diff = Inf;
for i = R.max_ind:-1:1
	temp = abs(R.data_sim(i) - middle_amplitude);

	if temp < min_diff
		mid_amp_ind = i;
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
	fwhm = F.t_sim(mid_amp_ind2) - F.t_sim(mid_amp_ind);
end

%% Populate status matrix

% Calculate onset time
if fwhm < min_half_width
	% No significant signal
	t_onset = NaN;
else

	% Calculate temporal derivative
	start_deriv = deriv(2);
	end_deriv = deriv(end);

	if R.max_ind <= R.min_ind || start_deriv < 0 || start_deriv > assume_const_delta || max(deriv) < assume_const_delta
		% Cell died before measurement
		t_onset = -Inf;

	elseif end_deriv > assume_const_delta || min(deriv) > assume_const_delta
		% Cell died to late; reaches maximum after measurement
		t_onset = Inf;

	else
		% Cell dies within observed time

		% Extrapolate ascent
		% DEBUG OUTPUT
		if (~(isfinite(mid_amp_ind) && mid_amp_ind > 1))
			disp(['Invalid mid_amp_ind for cell ' num2str(D.index) ' in file ' F.fullpath ': mid_amp_ind=' num2str(mid_amp_ind)]);
		end

		% Define coefficients for rising edge
		mid_amp_deriv = ...
			( F.t_sim(mid_amp_ind + 1) - F.t_sim(mid_amp_ind - 1) ) ...
			/ ( R.data_sim( mid_amp_ind + 1 ) - R.data_sim( mid_amp_ind - 1 ) );

		rising_edge_xoff = middle_amplitude;
		rising_edge_yoff = F.t_sim(mid_amp_ind);
		rising_edge_scale = mid_amp_deriv;

		%ascent_time = @(Y) ( Y - middle_amplitude ) * mid_amp_deriv  ...
		%	+ F.t_sim( mid_amp_ind );
		ascent_time = @(Y) ( Y - rising_edge_xoff ) * rising_edge_scale  ...
			+ rising_edge_yoff;

		% Calculate onset time
		t_onset = ascent_time(R.data_sim(1));

	end
end

%% Return values

% Onset time
varargout{1} = t_onset;

% Fit type
if nargout > 1
	varargout{2}.type = int8(1); % type: extrapol_rising
	varargout{2}.rising_edge_xoff = rising_edge_xoff;
	varargout{2}.rising_edge_yoff = rising_edge_yoff;
	varargout{2}.rising_edge_scale = rising_edge_scale;
end

% Initial derivative
if nargout > 2
	if isfinite(t_onset)
		varargout{3} = 1 / mid_amp_deriv;
	else
		varargout{3} = 1 / deriv(rising_ind);
	end
end
