%% postproc_ros.m
% This script contains the postprocessing routine for CellROX signal.
%
% Required variables:
% ===================
% data					array of measured data
% t						time vector of measured data (in hours)
% data_sim				array of simulated data
% t_sim					time vector of simulated data (in hours)
% max_val				vector of maximum values of all traces
% max_ind				vector of indices of maxima in all traces
% min_val				vector of minimum values of all traces
% min_ind				vector of indices of minima in all traces
% amplitude				vector of amplitudes of all traces
% outmode.debug			boolean whether debug figure should be created
% outmode.states		boolean whether state values should be created
% model.marker			marker model name
% t_unit				time unit (for plot)
% target_dir			directory in which to save data
% filename				name of the initial file
% format_cell_number	formatted string of cell number
% time_now				formatted string of current time
% 
% Provided variables:
% ===================
% stat_mat				status matrix; contains onset times and amplitudes
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

%% Allocate memory
stat_mat = zeros(size(data,2), 5);

%% Define thresholds
amp_max = max(amplitude);
noise = amp_max * .05;
assume_const_delta = 1E-2;
slope_threshold = -50;

%% Populate status matrix
parfor c = 1:size(data,2)

	% Set default values
	norm_deriv_min_val = NaN;
	deriv_min_ind = NaN;

	% Calculate onset time
	if amplitude(c) < noise
		% Cell does not die at all
		t_breakdown = NaN;
	else
		% Calculate temporal derivative
		start_deriv = data_sim(2,c) - data_sim(1,c);
		end_deriv = data_sim(length(t_sim),c) - data_sim(length(t_sim)-1,c);

		if start_deriv <= assume_const_delta
			% Cell died before measurement
			t_breakdown = -Inf;
		elseif ( end_deriv > assume_const_delta && ~(start_deriv > 0 && min_ind(c) ~= 1 && min_ind(c) ~= size(data_sim,1)) ) ...
				|| end_deriv < -assume_const_delta || data_sim(size(data_sim,1),c) == max_val(c)
			% Cell died to late
			t_breakdown = Inf;
		else
			% Cell dies within observed time
			first_val = data_sim(1,c);

			% Get 1st derivative of flourescence signal w.r.t. time
			deriv = zeros(1,size(data_sim,1));
			deriv(1) = NaN;

			for i = 2:size(data_sim,1)
				deriv(i) = ( data_sim(i,c) - data_sim(i-1,c) ) ...
					/ ( t_sim(i) - t_sim(i-1) );
			end
			
			% Get normalized 1st derivative
			norm_deriv = deriv / amplitude(c);

			% Find slope minimum
			[~,deriv_min_ind] = min(deriv);
			norm_deriv_min_val = norm_deriv(deriv_min_ind);

			% Set up extrapolation function
			breakdownfunc = @(y) (y - data_sim(deriv_min_ind,c)) / deriv(deriv_min_ind) ...
				+ t_sim(deriv_min_ind);
			
			% Calculate onset time by extrapolating from slope minimum
			t_breakdown = breakdownfunc(max_val(c));
		end
	end

	% Write status matrix row:
	%	1st column: cell number
	%	2nd column: onset time
	%	3rd column: normalized minimal slope
	%	4th column: absolute amplitude
	%	5th column: relative amplitude
	stat_mat(c,:) = [ c t_breakdown norm_deriv_min_val amplitude(c) amplitude(c)/amp_max ];

	%% Plot debug figures
	if outmode.debug
		fh = figure;
		p = plot(t_sim,data_sim(:,c),'r-', t,data(:,c),'b-');
			% Minimum amplitude level:
			line(xlim, [min_val(c)+noise min_val(c)+noise], 'LineStyle','-.', 'Color',[1,.5,0]);
			if isfinite(t_breakdown)
				hold on;
				plot(t_sim,norm_deriv,'g-');
				% Onset time:
				line([t_breakdown t_breakdown], ylim, 'LineStyle','-', 'Color',[1,.5,0]);
				% Minimum of slope:
				line([t_sim(deriv_min_ind), t_sim(deriv_min_ind)], ylim, 'LineStyle',':', 'Color',[1,.5,0]);
				% Extrapolation:
				line(breakdownfunc(ylim), ylim, 'LineStyle','-', 'Color',[1,.5,0]);
			end
		uistack(p, 'top');
		title([model.marker ' Fluorescence in Cell ' num2str(c)]);
		xlabel(['Time ' t_unit]); xlim([0,max(t)]);
		ylabel('Fluorescence Intensity [a.u.]');
		set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
		print(fh, '-dpdf', fullfile(target_dir, [filename '_DEBUG_' num2str(c, format_cell_number) '_' time_now '.pdf']));
		close(fh);
	end
end

%% Export status matrix
if outmode.states
	csvwrite(fullfile(target_dir, [filename '_ALL_STATE_' time_now '.csv']), stat_mat);
end
