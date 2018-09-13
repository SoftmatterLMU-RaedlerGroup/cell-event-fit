function stat_mat = postproc_ros_kink(D, F, output_info)
% This file contains the postprocessing routine for CellROX signal.
%
% Required arguments:
% ===================
% D					structure, containing these fields:
% D.data					array of measured data
% D.data_sim				array of simulated data
% D.max_val					maximum value
% D.max_ind					index of maximum
% D.min_val					minimum value
% D.min_ind					index of minimum
% D.amplitude				amplitude (D.max_val - D.min_val)
% D.amp_max					file-wide maximum value of amplitude
% D.params					vector of fit parameters
% D.logPost					log-likelihood for parameter set D.params
% D.index					index of current trace in datafile
% F					structure, containing these fields:
% F.t						time vector of measured data (in hours)
% F.t_sim					time vector of simulated data (in hours)
% F.t_unit					time unit (for plot)
% F.target_dir				directory in which to save data
% F.name					name of the initial file
% F.format_cell_number		formatted string of cell number
% output_info		structure, containing these fields:
% output_info.debug			boolean whether debug figure should be created
% output_info.marker		marker model name
% output_info.time_now		formatted string of current time
%
% Return value:
% =============
% stat_mat			status (line) vector; contains onset times and amplitudes
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
noise = D.amp_max * .1;
assume_const_delta = 1E-2;
%slope_threshold = 0.01;
parab_threshold = 0.01;

%% Populate status matrix
% Set default values
parab_threshold = parab_threshold * D.amplitude;
%norm_deriv_min_val = NaN;
deriv_min_ind = NaN;
t_breakdown = NaN;
switch_to_extrapolate = 0;
parabola = zeros(length(F.t_sim), 1);
parabola(:) = NaN;
fit_type = 0;
init_slope = NaN;

% Evaluate differences to parabola
for i = 1:length(F.t_sim)
	parabola(i) = D.params(1) + D.params(2) ...
		* ( D.params(7) - D.params(5) * (F.t_sim(i) - D.params(6))^2 );

	if parabola(i) < D.min_val - parab_threshold
		parabola(i) = NaN;
	end
end

% Calculate onset time
if D.amplitude < noise
	% Cell does not die at all
	t_breakdown = NaN;
else
	% Calculate temporal derivative
	start_deriv = D.data_sim(2) - D.data_sim(1);
	end_deriv = D.data_sim(length(F.t_sim)) - D.data_sim(length(F.t_sim)-1);
	init_slope = start_deriv / (F.t_sim(2) - F.t_sim(1));

	if start_deriv <= assume_const_delta
		% Cell died before measurement
		t_breakdown = -Inf;
	elseif ( end_deriv > assume_const_delta && ~(start_deriv > 0 && D.min_ind ~= 1 && D.min_ind ~= size(D.data_sim,1)) ) ...
			|| end_deriv < -assume_const_delta || D.data_sim(size(D.data_sim,1)) == D.max_val
		% Cell died to late
		t_breakdown = Inf;
	else
		% Cell dies within observed time
		%first_val = D.data_sim(1);

		% Get 1st derivative of flourescence signal w.r.t. time
		deriv = zeros(1,size(D.data_sim,1));
		deriv(1) = NaN;

		for i = 2:size(D.data_sim,1)
			deriv(i) = ( D.data_sim(i) - D.data_sim(i-1) ) ...
				/ ( F.t_sim(i) - F.t_sim(i-1) );
		end

% FIRST TRY; DISCARDED
% 		% Find the time when slope becomes non-linear
% 		slope_lin_slope = (norm_deriv(3) - norm_deriv(2)) ...
% 			/ (F.t_sim(3) - F.t_sim(2));
% 		slope_offset = norm_deriv(2);
%
% 		for i = 4:length(norm_deriv);
% 			if abs( slope_offset + (i-1) * slope_lin_slope - norm_deriv(i) ) >= slope_threshold
% 				t_breakdown = F.t_sim(i);
%				break
% 			end
% 		end
%
% SECOND TRY; ALSO DISCARDED
% 		% Find the first kink in the derivative
% 		for i = deriv_min_ind:-1:3
% 			if abs( (norm_deriv(i) - norm_deriv(i-1)) ...
% 					/ (F.t_sim(i) - F.t_sim(i-1)) ) <= slope_threshold
% 				t_breakdown = F.t_sim(i);
% 				disp(['Cell ' num2str(c) ': t_breakdown = ' num2str(t_breakdown) '; F.t_sim(i) = ' num2str(F.t_sim(i))])
% 				break
% 			end
% 		end
%
% THIRD TRY
		% Locate kink by deviation from parabola
		for i = 1:length(F.t_sim)
			parab_devi = abs(parabola(i) - D.data_sim(i));

			if isnan(parab_devi) || parab_devi <= parab_threshold
				% Nothing interesting happens
				continue;

			elseif deriv(i) >= 0 && i ~= D.max_ind
				% Too flat breakdown, switch to extrapolation mode
				switch_to_extrapolate = 1;
				break;

			else
				if i > D.max_ind && D.data_sim(i-1) <= D.min_val + noise
					% Signal is too low, still no breakdown, switch
					% to extrapolation mode
					switch_to_extrapolate = 1;
					break;
				elseif i > 1
					% Breakdown time found
					t_breakdown = F.t_sim(i);
					fit_type = 1;
					%disp(['Cell ' num2str(D.index) ': t_breakdown=' num2str(t_breakdown) ...
					%	' data_sim=' num2str(D.data_sim(i)) ' parab=' num2str(parabola(i)) ]);
					break;
				else
					continue;
				end

			end
		end


		if switch_to_extrapolate
			% Get normalized 1st derivative
			%norm_deriv = deriv / D.amplitude;

			% Find slope minimum
			[~,deriv_min_ind] = min(deriv);
			%norm_deriv_min_val = norm_deriv(deriv_min_ind);

			% Set up extrapolation function
			breakdownfunc = @(y) (y - D.data_sim(deriv_min_ind)) / deriv(deriv_min_ind) ...
				+ F.t_sim(deriv_min_ind);

			% Calculate onset time by extrapolating from slope minimum
			t_breakdown = breakdownfunc(D.max_val);
			fit_type = -1;
		end
	end
end

% Write status matrix row:
%	1st column: cell number
%	2nd column: breakdown time
%	3rd column: absolute amplitude
%	4th column: relative amplitude
%	5th column: log-likelihood
%	6th column: fit_type (1: kink, 0: no fit, -1: extrapolation)
%	7th column: initial slope
stat_mat = [ D.index t_breakdown D.amplitude D.amplitude/D.amp_max D.logPost fit_type init_slope ];

%% Plot debug figures
if output_info.debug
	fh = figure;
	plot(F.t_sim,D.data_sim,'r-', F.t,D.data,'b-');

	% Minimum amplitude level:
	line(xlim, [D.min_val+noise D.min_val+noise], 'LineStyle','-.', 'Color',[1,.5,0]);

	if isfinite(t_breakdown)
		hold on;

		if switch_to_extrapolate
			%plot(F.t_sim,norm_deriv,'g-');
			% Extrapolation:
			line(breakdownfunc(ylim), ylim, 'LineStyle','-', 'Color',[1,.5,0]);
			% Minimum of slope:
			line([F.t_sim(deriv_min_ind), F.t_sim(deriv_min_ind)], ylim, 'LineStyle',':', 'Color',[1,.5,0]);
		else
			plot(F.t_sim,parabola,'m-');
		end

		% Onset time:
		line([t_breakdown t_breakdown], ylim, 'LineStyle','-', 'Color',[1,.5,0]);
	end

	%uistack(p, 'top');
	title([output_info.marker ' Fluorescence in Cell ' num2str(D.index)]);
	xlabel(['Time ' F.t_unit]); xlim([0,max(F.t)]);
	ylabel('Fluorescence Intensity [a.u.]');
	set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
	print(fh, '-dpdf', fullfile(F.target_dir, [F.name '_DEBUG_' num2str(D.index, F.format_cell_number) '_' output_info.time_now '.pdf']));
	close(fh);
end
