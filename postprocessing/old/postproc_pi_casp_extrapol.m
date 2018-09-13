%% postproc_pi_casp_extrapol.m
% This script contains the postprocessing routine for PI and Caspase signal.
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
stat_mat = zeros(size(data,2), 4);

%% Define thresholds
amp_max = max(amplitude);
noise = amp_max * .05;
assume_const_delta = 1E-2;

%% Populate status matrix
for c = 1:size(data,2)

	% Calculate onset time
	if amplitude(c) < noise
		% Cell does not die at all
		t_onset = NaN;
	else
		% Calculate temporal derivative
		start_deriv = data_sim(2,c) - data_sim(1,c);
		end_deriv = data_sim(length(t_sim),c) - data_sim(length(t_sim)-1,c);

		if max_ind(c) <= min_ind(c) || start_deriv < 0 || start_deriv > assume_const_delta
			% Cell died before measurement
			t_onset = -Inf;
		elseif end_deriv > assume_const_delta
			% Cell died to late; reaches maximum after measurement
			t_onset = Inf;
		else
			% Cell dies within observed time
			first_val = data_sim(1,c);

			% Extrapolate ascent
			%middle_amplitude = min_val(c) + amplitude(c) / 2;
			middle_amplitude = model.simulate(params(c,3),params(c,:));
			[~,mid_amp_ind] = min(abs(data_sim(min_ind(c):max_ind(c),c) - middle_amplitude));
			mid_amp_val1 = data_sim(mid_amp_ind-1,c);
			mid_amp_val2 = data_sim(mid_amp_ind+1,c);
			mid_amp_time = t_sim(mid_amp_ind);

			ascent_time = @(Y) (Y-middle_amplitude) / (mid_amp_val2-mid_amp_val1) * 2 * t_sim_res + mid_amp_time;

			% Calculate onset time
			t_onset = ascent_time(first_val);
		end
	end

	% Write status matrix row:
	%	 1st column: cell number
	%	 2nd column: onset time
	%	 3rd column: absolute amplitude
	%	 4th column: relative amplitude
	stat_mat(c,:) = [ c amplitude(c) amplitude(c)/amp_max t_onset ];

	%% Plot debug figures
	if outmode.debug
		fh = figure;
		p = plot(t_sim,data_sim(:,c),'r-', t,data(:,c),'b-');
			% Minimum amplitude level:
			line(xlim, [min_val(c)+noise min_val(c)+noise], 'LineStyle', '-.', 'Color', [1,.5,0]);
			if isfinite(t_onset)
				% Horizontal line:
				line(xlim, [first_val first_val], 'LineStyle', '-', 'Color', [1,.5,0]);
				% Ascent line:
				yl = ylim;
				line([ascent_time(yl(1)) ascent_time(yl(2))], ylim, 'LineStyle', '-', 'Color', [1,.5,0]);
				% Onset time:
				line([t_onset t_onset], ylim, 'LineStyle', ':', 'Color', [1,.5,0]);
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
