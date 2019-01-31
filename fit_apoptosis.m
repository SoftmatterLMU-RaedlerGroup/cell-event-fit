function [] = fit_apoptosis(varargin)
%fit_apoptosis fits fluorescece time traces and extracts event times
%
% Input parameters:
%	varargin(1):	name of model (marker type)
%	varargin(2):	cell array of paths to datafiles
%	varargin(3):	path for results to be saved
%	varargin(4):	print/export options
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

%% Read in data from source files
disp([ get_time 'Started.' ])
Ffile = readInputData(varargin{:});
mf = matfile(Ffile);
disp([get_time 'Read in ' num2str(mf.ntraces) ' ' mf.log_trace ' from ' num2str(mf.ndatafiles) ' ' mf.log_file '.'])

% Close matfile
clear mf;

%% Perform fit, iterating over traces in all files
disp([ get_time 'Starting fitting of all traces ...' ])

Ffile = performParallelTask(Ffile);

%% Export results
disp([ get_time 'Exporting data ...' ])
plotFitResults(Ffile);

%% Get data for final log message
mf = matfile(Ffile);
ntraces = mf.ntraces;
log_trace = mf.log_trace;
ndatafiles = mf.ndatafiles;
log_file = mf.log_file;

%% Remove temporary files
outmode = mf.outmode;
clear mf
if ~outmode.temp
	delete(Ffile);
end

%% Print final log message
disp([get_time 'FINISHED – ' num2str(ntraces) ' ' log_trace ' in ' num2str(ndatafiles) ' ' log_file ' processed.'])

end