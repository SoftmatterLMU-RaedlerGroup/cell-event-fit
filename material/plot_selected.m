% plot_selected.m
% This script helps to plot special traces from different files.
% It takes as input the cell array tracesToPlot, which can
% for example be created by plot_select.m, and the variable
% verbose_plot, which switches on/off verbose plotting.
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

%% Get plot options
if ~exist('verbose_plot', 'var')
	verbose_plot = true;
end

%% Create figure with desired options
fh = figure('PaperUnits','centimeters', 'PaperSize',[25 14], ...
	'PaperPosition',[0 0 25 14], 'PaperPositionMode','manual');
clmp = colormap(fh, lines(length(tracesToPlot)));

%% Iterate over traces and plot them
for iline = 1:length(tracesToPlot)
	for itrace = 1:length(tracesToPlot{iline})
		% Get trace struct
		trace = tracesToPlot{iline}(itrace);

		% Plot data
		plot(trace.t, trace.data, '-', 'linewidth', 1, 'color', clmp(iline,:))

		hold on;

		% Plot simulated data
		plot(trace.t_sim, trace.data_sim, '-', 'linewidth', 2, 'color', clmp(iline,:));

		% Break here, if non-verbose plotting desired
		if ~verbose_plot
			continue
		end

		% Plot noise level
		line(xlim, [trace.noise, trace.noise], 'linestyle', ':', 'linewidth', 1, 'color', clmp(iline,:));

		% Plot parabola (if exists), else plot baseline
		if ~isempty(trace.parab)
			plot(trace.t_sim, trace.parab, ':', 'linewidth', 1, 'color', clmp(iline,:));
		else
			line(xlim, [trace.data_sim(1), trace.data_sim(1)], 'linestyle', '-', 'linewidth', 1, 'color', clmp(iline,:));
		end

		% Plot slope
		yl = ylim;
		line([trace.slope(yl(1)), trace.slope(yl(2))], ylim, 'linestyle', '-', 'linewidth', 1, 'color', clmp(iline,:));

		% Plot event time
		line([trace.t_event, trace.t_event], ylim, 'linestyle', '--', 'linewidth', 1, 'color', clmp(iline,:));
	end
end

% Get labelling
xlabel('Time [#]');
ylabel('Fluorescence Intensity [a.u.]');
user_title = inputdlg('Please enter a figure title:', 'Title required');
if ~isempty(user_title)
	user_title = user_title{1};
	title(user_title);
else
	user_title = '';
end

% Save figure (to be done manually)
%print(fh, 'example.eps', '-depsc2', '-cmyk')
fprintf('To save your figure, use:\nprint(fh, ''%s.eps'', ''-depsc2'', ''-cmyk'')\n', user_title);