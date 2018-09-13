% plot_select.m
% This script helps to select special traces from different matfiles.
% It creates the cell array tracesToPlot, which can be used as an input
% to plot_selected.m.
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

%% Prompt for source mat-file and open it
[name,path] = uigetfile('*.mat', 'Select the source file:');
mf = matfile( fullfile(path, name) );

%% Prompt for files to be used
sourceFiles = cell( mf.ndatafiles, 1 );

for f = 1:length(sourceFiles)
	sourceFiles{f,:} = mf.name_F(mf.name_F_ind(f,1),:);
end

fileIndices = listdlg('ListString', sourceFiles, 'ListSize', [500,300], ...
	'PromptString', 'Please select the files you want to use:');

%% Prompt for traces to be used
traceIndices = uint16([]);

for ifile = fileIndices
	tracesInFile = cell(mf.size_F(ifile, 1), 1);
	for itrace = 1:mf.size_F(ifile, 1)
		tracesInFile{itrace} = num2str(mf.index_F( mf.dataindices(ifile, itrace), 1 ));
	end

	tracesFromFile = listdlg('ListString', tracesInFile, ...
		'PromptString', {'Select traces from file:'; ...
		sourceFiles{ifile,:}});

	if isempty(tracesFromFile)
		continue
	end

	for itrace = tracesFromFile
		traceIndices = [ traceIndices mf.dataindices(ifile, itrace) ];
	end
end

%% Clear variables that are not needed any more
clear f fileIndices ifile itrace name path tracesInFile tracesFromFile sourceFiles;

%% Write relevant information to final cell array
if ~exist('tracesToPlot', 'var')
	itraces = 1;
else
	itraces = length(tracesToPlot) + 1;
end

tracesToPlot{itraces}(1,1:length(traceIndices)) = struct( ...
	't',[], 'data',[], 't_sim',[], 'data_sim',[], ...
	't_event',0, 'noise',0, 'parab',[], 'slope',[] ...
	);

for itrace = 1:length(traceIndices)
	trace = traceIndices(itrace);
	file = mf.file_ind(trace, 1);
	tracesToPlot{itraces}(itrace).t = mf.t(1:mf.data_len(file,1),mf.t_ind(file,1));
	tracesToPlot{itraces}(itrace).t_sim = mf.t_sim(1:mf.data_sim_len(file,1),mf.t_sim_ind(file,1));
	tracesToPlot{itraces}(itrace).data = mf.data(1:mf.data_len(file,1),trace);
	tracesToPlot{itraces}(itrace).data_sim = mf.data_sim(1:mf.data_sim_len(file,1),trace);
	tracesToPlot{itraces}(itrace).t_event = mf.t_event(trace,1);
	tracesToPlot{itraces}(itrace).noise = mf.noise(file,1) + mf.min_val(trace,1);
	parab_ind = mf.fit_parabola_ind(trace,1);
	if parab_ind > 0
		tracesToPlot{itraces}(itrace).parab = mf.fit_parabola(1:mf.data_sim_len(file,1),parab_ind);
	end
	tracesToPlot{itraces}(itrace).slope = @(y) ...
		mf.rising_edge_scale(trace,1) * (y - mf.rising_edge_xoff(trace,1)) + ...
		mf.rising_edge_yoff(trace,1);
end