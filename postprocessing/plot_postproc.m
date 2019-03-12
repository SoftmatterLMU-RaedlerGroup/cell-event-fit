function plot_postproc(info)
% This function plots the postprocessing results in a figure.
%
% Required arguments:
% ===================
% info					structure containing information needed for plot
% info.t					time vector
% info.t_sim				simulated time vector
% info.t_unit				time unit of info.t/info.t_sim as a string (for plot)
% info.data					vector of measured data to be plotted
% info.data_sim				vector of simulated data of current trace
% info.noise				noise level of source data file
% info.t_event				time of event (onset, breakdown)
% info.fit_type				number of fit type used
% info.rising_edge_xoff		abscissa offset of rising edge interpolation
% info.rising_edge_yoff		ordinate offset of rising edge interpolation
% info.rising_edge_scale	scaling of rising edge interpolation
% info.falling_edge_xoff	abscissa offset of falling edge interpolation
% info.falling_edge_yoff	ordinate offset of falling edge interpolation
% info.falling_edge_scale	scaling of falling edge interpolation
% info.fit_parabola			parabola to be plotted (model-dependent)
% info.marker				string of marker name
% info.index				index of current trace
% info.target_dir			string of directory to save plots in
% info.name					name of source data file (without extension)
% info.format_cell_number	formatstring for plot file numbering
% info.time_now				time fingerprint of current session
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

%% Define extrapolation functions
rising_edge = @(x) (x - info.rising_edge_xoff) * info.rising_edge_scale ...
					+ info.rising_edge_yoff;

falling_edge = @(x) (x - info.falling_edge_xoff) * info.falling_edge_scale ...
					+ info.falling_edge_yoff;

%% Plot debug figures
% Plot measured data
fh = figure('Visible','off');
ax = axes(fh);
hold(ax, 'on');
plot(ax, info.t, info.data, 'b-');

% Minimum amplitude level:
line(ax, ax.XLim, [info.noise info.noise], 'LineStyle', ':', 'Color', [1,.5,0]);


switch info.fit_type

	case 1 % extrapol_rising

		if isfinite(info.t_event)
			% Horizontal line:
			line(ax, ax.XLim, [info.data_sim(1) info.data_sim(1)], 'LineStyle', '-', 'Color', [0,1,0]);
		end

		yl = ax.YLim;

		% Rising edge:
		line(ax, [rising_edge(yl(1)) rising_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');

		if isfinite(info.t_event)
			% Onset time:
			line(ax, [info.t_event info.t_event], ax.YLim, 'LineStyle', '--', 'Color', [0,.5,0]);
		end

	case -2 % nokink_parabola
		yl = ax.YLim;

		% Plot maximum slope of the rising edge
		line(ax, [rising_edge(yl(1)) rising_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);

		if isfinite(info.t_event)
			% Plot maximum slope of the falling edge
			line(ax, [falling_edge(yl(1)) falling_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);
		end

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');

		% Pure parabola
		plot(ax, info.t_sim, info.fit_parabola, 'LineStyle','-.', 'Color',[0,1,0]);

		if isfinite(info.t_event)
			% Onset time:
			line(ax, [info.t_event info.t_event], ax.YLim, 'LineStyle','--', 'Color',[0,.5,0]);
		end

	case 2 % kink_parabola
		yl = ax.YLim;

		% Plot maximum slope of the rising edge
		line(ax, [rising_edge(yl(1)) rising_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);

		if isfinite(info.t_event)
			% Plot maximum slope of the falling edge
			line(ax, [falling_edge(yl(1)) falling_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);
		end

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');

		% Pure parabola
		plot(ax, info.t_sim, info.fit_parabola, 'LineStyle','-.', 'Color',[0,1,0]);

		if isfinite(info.t_event)
			% Onset time:
			line(ax, [info.t_event info.t_event], ax.YLim, 'LineStyle','--', 'Color',[0,.5,0]);
		end

	case 3 % mRNA transfection

		if isfinite(info.t_event)
			% Horizontal line:
			line(ax, ax.XLim, [info.data_sim(1) info.data_sim(1)], 'LineStyle', '-', 'Color', [0,1,0]);
		end

		yl = ax.YLim;

		% Rising edge:
		line(ax, [rising_edge(yl(1)) rising_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');

		if isfinite(info.t_event)
			% Onset time:
			line(ax, [info.t_event info.t_event], ax.YLim, 'LineStyle', '--', 'Color', [0,.5,0]);
		end

	case 4 % min between peaks (see postproc_LATE_decay)
		yl = ax.YLim;

		% Rising edge:
		line(ax, [rising_edge(yl(1)) rising_edge(yl(2))], yl, 'LineStyle','-', 'Color',[0,1,0]);

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');

		if isfinite(info.t_event)
			% Onset time:
			line(ax, [info.t_event info.t_event], ax.YLim, 'LineStyle', '--', 'Color', [0,.5,0]);
		end

	case {5, -5} % interactive EARLY model

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');

		% Save y limits for restoring later
		yl = ax.YLim;

		% Plot parabola
		parabola = parabola_EARLY_simulate(info.t_sim, info.params, 'parabola');
		parabola((parabola > yl(2)) & (parabola < yl(1))) = NaN;
		plot(ax, info.t_sim, parabola, 'LineStyle', '-.', 'Color', [0,1,0]);

		% Display fit limits
		t_start = info.params(7);
		t_end = info.params(8);
		if isfinite(t_start) && t_start > min(info.t)
			line(ax, [t_start t_start], yl, 'LineStyle', ':', 'Color', 'k');
		end
		if isfinite(t_end) && t_end < max(info.t)
			line(ax, [t_end t_end], yl, 'LineStyle', ':', 'Color', 'k');
		end

		% Breakdown time
		if isfinite(info.t_event)
			line(ax, [info.t_event info.t_event], yl, 'LineStyle', '--', 'Color', [0,.5,0]);
			text(ax, info.t(end), yl(2), ...
				sprintf('t_{breakdown}=%.2fh', info.t_event), ...
				'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
				'FontSize', 8);
		end

		% Restore y limits
		ax.YLim = yl;

	otherwise
		% Fit type not defined
		disp([ 'Unknown fit type: ' num2str(info.fit_type) ]);

		% Plot simulation
		plot(ax, info.t_sim, info.data_sim, 'r-');
end

% Plot labels and print
title(ax, [info.marker ' Fluorescence in Cell ' num2str(info.index)]);
xlabel(ax, ['Time ' info.t_unit]);
ylabel(ax, 'Fluorescence Intensity [a.u.]');

ax.Box = 'on';
ax.XLim = [0, max(info.t)];

set(fh, 'PaperUnits','centimeters', 'PaperSize',[10 7], 'PaperPosition',[0 0 10 7])
print(fh, '-dpdf', fullfile(info.target_dir, strcat(info.name, '_DEBUG_', ...
	num2str(info.index, info.format_cell_number), '_', info.time_now, '.pdf')));
delete(fh);
