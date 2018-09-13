% A wrapper for non-interactive execution of `fit_apoptosis`
%
% It can be used for calling fit_apoptosis with the desired parameters.
% Use this file as a template and adjust it according to your needs.
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

%% (OPTIONAL) Prepare custom parallel pool options
% % Initialize parpool with `min_n_workers` workers
%min_n_workers = 20;
%try
%	pool = parpool(min_n_workers);
%catch
%end
%
% % If parpool initialization failed, fail with exit status 1
%if ~exist('p', 'var') || isempty(pool) || pool.NumWorkers < min_n_workers
%	disp("Parallel pool execution failed, quitting.");
%	exit(1)
%end

%% Specify input parameters for fit_apoptosis

% Name of the marker (model) to be used
markername = 'ros';

% Datafiles to be processed (must be a cell array)
% If the filename is given as an absolute file path, the path is used.
% If the filename is given as a relative path, the path is interpreted
% from the current working directory.
% If no directory but only the filename is given, it will be searched in
% the directory "project/data/", relative to the path of the script.
datafiles = { ...
	'file1', ...
	'file2' ...
};

% (OPTIONAL) Directory to write results to
% If the directory is not given as an absolute path, it will be
% interpreted as the current working directory of the script.
% If this option is not specified, a default directory hierarchy will
% be used.
%outdir = ''

%% Call fit_apoptosis.m
if exist('outdir', 'var')
	fit_apoptosis( markername, datafiles, outdir );
else
	fit_apoptosis( markername, datafiles );
end

%% Close matlab
exit(0)
