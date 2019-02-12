function [ Ffile ] = readInputData(varargin)
%readInputData reads in the given data and writes it in a *.mat-file.
%
% Input parameters
% ================
% varargin{1}:		name of model (marker type)
% varargin{2}:		cell array of paths to datafiles
% varargin{3}:		path for results to be saved
% varargin{4}:		cell array of options
% varargin{5}:		print/export options
%
% varargin{1} and varargin{3} can be a string or a cell array with a length
% equal to the length of varargin{2}.
%
% varargin{4} must be a cell array of the form:
%	{Name1, Value1, ... NameN, ValueN}
% Possible options are:
%	TempDir			Char vector with path to directory for temporary files
%	TimeUnit		Char vector with time unit for plot (default: '[h]')
%	TimeScale		Numeric vector for converting time vector to TimeUnit
%					(one global value or vector with an entry per data file)
%	TimeOffset		Numeric vector to be added to time vector (in TimeUnit)
%					(one global value or vector with an entry per data file)
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

%% Get input parameters
% Get number of input parameters
argc = length(varargin);

% Analyze input parameters
% Model
if argc >= 1
	modelname = varargin{1};

	if ~iscellstr(modelname)
		if ~ischar(modelname)
			error('The model name must be a char array or cell string.')
		end
		modelname = { modelname };
	end
else
	modelname = {};
end

if argc >= 2
	datafilenames = varargin{2};

	if ~iscellstr(datafilenames)
		if ~ischar(datafilenames)
			error('The data file names must be a char array or cell string.')
		end
		datafilenames = { datafilenames };
	end
else
	datafilenames = {};
end

% Directory to write results to
if argc >= 3
	out_dir = varargin{3};

	if isempty(out_dir)
		out_dir = {};
	elseif ischar(out_dir)
		out_dir = { out_dir };
	elseif ~iscellstr(out_dir)
		error('The output directory must be a char vector or cell string.')
	end
else
	out_dir = {};
end

if argc >= 4
	% varargin{4} is the options cell array
	options_in = varargin{4};
	if ~iscell(options_in)
		error('Options must be given as a Key-Value cell array.')
	elseif mod(length(options_in), 2) ~= 0
		error('Options cell array must have an even number of elements.')
	end
else
	options_in = {};
end


if argc >= 5
	outmode_in = varargin{5};
else
	outmode_in = {};
end


%% Get model
if isempty(modelname)
	modelname = {'Lyso', 'TMRM', 'TMRM_interactive', 'CellROX', 'Calcium520', 'Caspase', 'pSIVA', 'PI', 'mRNA', 'Test', 'TMRM_new'};
	modelindex = listdlg('Name', 'Model Missing', ...
						'PromptString', 'Please select a model:', ...
						'SelectionMode', 'single', ...
						'ListSize', [160 115], ...
						'ListString', modelname);
	if isempty(modelindex)
		error('No model name given.');
	else
		modelname = modelname(modelindex);
	end

else
	if length(modelname) > 1
		if argc < 2 || length(datafilenames) ~= length(modelname)
			error(['If there is more than one model given, ' ...
				'the files must also be given as input parameters. ' ...
				'The number of files specified must be equal to the number of models.'])
		end

		if argc > 2 && length(out_dir) ~= length(modelname)
			error(['If there is more than one model given, ' ...
				'the number of output directories must match ' ...
				'the number of models.'])
		end
	end
end

%% Get datafiles
% Get path of this script
[script_path,~,~] = fileparts( mfilename('fullpath') );

if isempty(datafilenames)
	% Ask interactively for filenames
	[name, path] = uigetfile( ...
		{'*.csv;*.txt;*.xls;*.xlsx;*.ods', 'All Tables'; ...
			'*.csv', 'Comma-Separated Values Tables'; ...
			'*.txt', 'Separated Values Tables'; ...
			'*.xls;*.xlsx', 'Microsoft Excel Tables'; ...
			'*.ods', 'OpenDocument Tables'; ...
			'*.*', 'All Files'}, ...
		'Select your data files', ...
		fullfile(script_path, '..', 'project', 'data', filesep), ...
		'MultiSelect','on');

	% Ensure the filenames are in correct format (absolute path, cell array)
	if ~iscellstr(name)
		name = { name };
	end
	datafilenames = fullfile(path, name);
end

% Throw an error if no input files are given
ndatafiles = uint16( length(datafilenames) );

if ndatafiles < 1
	error('No datafiles given.');
end

% Ensure that path numbers given for out_dir and datafilenames match
%	The following cases are allowed:
%	length(out_dir) == length(datafilenames)	% 1:1 match
%	length(out_dir) == 1	% all output to same directory
%	length(out_dir) == 0	% output directory not yet defined, use default
if length(out_dir) ~= ndatafiles ...
		&& length(out_dir) ~= 1 ...
		&& ~isempty(out_dir)
	error(['If there are specified several output directories, ' ...
		' the number of output directories must match the number of input files.'])
end

% Ensure existence of data file names
for i = 1:ndatafiles
	[path, name, extension] = fileparts(datafilenames{i});

	% If a relative file path was given, assume a default directory
	if strcmp(path, '')
	    datafilenames{i} = fullfile(script_path, '..', 'project', 'data', [name extension]);
	end

	% Throw an error if the file does not exist
	if ~exist(datafilenames{i}, 'file')
	    error(['File does not exist: ' datafilenames{i}]);
	end
end

%% Parse options
% Define default values
temp_dir = [];
t_unit = '[h]';
t_scale = [];
t_offset = 0;

for i_opt = 1:2:length(options_in) - 1

	if strcmpi(options_in{i_opt}, 'TempDir')
		% Directory for temporary files given
		temp_dir = options_in{i_opt + 1};

		if ~ischar(temp_dir)
			if ~iscellstr(temp_dir)
				error('The directory for temporary files must be a char array or cell string.')
			elseif length(temp_dir) > 1
				error('Only one directory for temporary files can be specified.')
			else
				temp_dir = temp_dir{1};
			end
		end

	elseif strcmpi(options_in{i_opt}, 'TimeUnit')
		% Time unit to be used for output
		t_unit = options_in{i_opt + 1};

		if ~ischar(t_unit)
			error('TimeUnit must be given as char vector.')
		elseif isempty(t_unit)
			error('TimeUnit is empty.')
		end

	elseif strcmpi(options_in{i_opt}, 'TimeScale')
		% Factor to be multiplied to time values for conversion into hours
		t_scale = options_in{i_opt + 1};

		if ~isnumeric(t_scale) || ~isvector(t_scale)
			error('TimeScale must be a numeric vector.')
		elseif length(t_scale) ~= ndatafiles && length(t_scale) ~= 1
			error(['TimeScale must be a vector of length 1 ' ...
				'or equal to the number of data files.'])
		end

	elseif strcmpi(options_in{i_opt}, 'TimeOffset')
		% Offset to be added to time vector (in units t_unit)
		t_offset = options_in{i_opt + 1};

		if ~isnumeric(t_offset) || ~isvector(t_offset)
			error('TimeOffset must be a numeric vector.')
		elseif length(t_offset) ~= 1 && length(t_offset) ~= ndatafiles
			error(['The length of TimeOffset must be 1 ' ...
				'or equal to the number of data files.'])
		end
	end
end

%% Get output/export options
outmode = struct( ...
	'debug',		true, ...	% Plot debug graph (model-specific)
	'ms',			false, ...	% Plot multistart graph
	'single',		true, ...	% Plot single cell fits
	'total',		true, ...	% Plot all cell data with fits
	'params',		true, ...	% Export parameters table
	'states',		true, ...	% Export states table
	'simulated',	true, ...	% Export simulated values
	'temp',			true ...	% Leave temporary file in directory
	);

if ~isempty(outmode_in)
	outmodes = fieldnames(outmode);

	if ischar(outmode_in)
		% Input argument is char array
		if strcmpi(outmode_in, 'all')
			% Set all values to TRUE
			for m = 1:numel(outmodes)
				outmode.(outmodes{m}) = 1;
			end
		else
			% Set the given value to TRUE, all others to FALSE
			for m = 1:numel(outmodes)
				outmode.(outmodes{m}) = strcmpi(outmodes{m}, outmode_in);
			end
		end
	elseif iscellstr(outmode_in)
		% Input argument is cell string
		for m = 1:numel(outmodes)
			% Set matching values to TRUE, all others to FALSE
			if any(strcmpi(outmodes{m}, outmode_in))
				outmode.(outmode_in{m}) = 1;
			else
				outmode.(outmode_in{m}) = 0;
			end
		end
	elseif isstruct(outmode_in)
		% Input argument is struct
		% Compare with possible outmode values
		outmodes_in = fieldnames(outmode_in);
		for i_om_in = 1:length(outmodes_in)
			for i_om = 1:length(outmodes)
				if strcmpi(outmodes_in{i_om_in}, outmodes{i_om})
					outmode.(outmodes{i_om}) = outmode_in.(outmodes_in{i_om_in});
				end
			end
		end
	else
		% No valid input is given
		error('The mode must be given as char array, cell string, or struct.');
	end
end

% Suppress annoying warnings
warning('off', 'MATLAB:nearlySingularMatrix'); %fmincon
warning('off', 'MATLAB:singularMatrix'); %fmincon
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary'); %parfor

%% Define model specifications
has_interactivity = false;
model_ind(1:ndatafiles,1) = uint8(1);

for imodelname = 1:length(modelname)
	exists_already = false;
	model = model_definitions(modelname{imodelname});

	if exist('models','var')
		% Test if the model exists already in the models struct array
		for imodel = 1:length(models)
			if isequaln(model, models(imodel))
				model_index = imodel;
				exists_already = true;
				break
			end
		end

		if ~exists_already
			model_index = length(models) + 1;
		end
	else
		model_index = 1;
	end

	if ~exists_already
		models(model_index, 1) = model;

		if isa(model.interactive, 'function_handle')
			has_interactivity = true;
		end
	end

	model_ind(imodelname) = model_index;
end

%% Get output directories
% Get current time for unique naming
date_now = datestr(datetime, 'yyyy-mm-dd');
time_now = clock;
time_now = [num2str(time_now(4), '%02d'), num2str(time_now(5), '%02d')];

% If no output directory specified, fallback to default directory
if isempty(out_dir)
	%out_dir = fullfile(script_path, 'project', 'results', filesep);
	out_dir{1} = fullfile(script_path, '..', 'project', 'results', models(1).name, date_now, filesep);
end

% Get directory for temporary files
if isempty(temp_dir)
	if length(out_dir) == 1
		temp_dir = out_dir{1};
	else
		error('No directory for temporary files given.')
	end
end

% Ensure that the temporary file directory exists
if ~exist(temp_dir, 'dir')
	[mkdir_success,mkdir_msg,~] = mkdir(temp_dir);
	if mkdir_success
		clear mkdir_success mkdir_msg;
	else
		error(['Cannot create directory “' temp_dir '”:' ...
			java.lang.System.getProperty('line.separator').char ...
			mkdir_msg])
	end
end

%% Read in the datafiles and write them to temporary matfile

% Define file-wise indexed non-unique arrays
full_path = [];
name_F = [];
extension_F = [];
target_dir = [];
model_name = [];

t = [];
t_sim = [];
format_cell_number = [];

% Define file-wise indexed unique arrays
full_path_ind(ndatafiles,1) = uint16(0);
name_F_ind(ndatafiles,1) = uint16(0);
name_F_len(ndatafiles,1) = uint16(0);
extension_F_ind(ndatafiles,1) = uint8(0);
target_dir_ind(ndatafiles,1) = uint16(0);
target_dir_len(ndatafiles,1) = uint16(0);
model_name_ind(ndatafiles,1) = uint8(0);
model_name_len(ndatafiles,1) = uint8(0);

data_len(ndatafiles,1) = uint16(0);
data_sim_len(ndatafiles,1) = uint16(0);
t_ind(ndatafiles,1) = uint16(0);
t_sim_ind(ndatafiles,1) = uint16(0);
t_sim_res(ndatafiles,1) = 0;
format_cell_number_ind(ndatafiles,1) = uint16(0);
par_num = zeros(ndatafiles, 1, 'uint8');

size_F(ndatafiles,1) = uint32(0);
dataindices(ndatafiles,1) = uint32(0);

% Define trace-wise indexed arrays
data = [];
file_ind(1:ndatafiles,1) = uint16(0);
index_F(1:ndatafiles,1) = uint32(0);

% Define traces counters
ntraces = uint32(0);

% Define global maximum parameter number
max_par_num = uint8(0);

for f = 1:ndatafiles

	% Get data file name information
	[full_path, full_path_ind(f)] = mssina(full_path, datafilenames{f});
	disp([get_time 'Reading file ' num2str(f) ': ' full_path(full_path_ind(f),:)])
	[~,temp_name,temp_ext] = fileparts(full_path(full_path_ind(f),:));

	% Remove trailing NUL bytes from extension (if any)
	%[~,~,temp_ext_ind] = regexpi(temp_ext, '^([^\0]*)\0*$');
	%temp_ext = temp_ext( temp_ext_ind{1}(1) : temp_ext_ind{1}(2) );
	temp_ext = strcat(temp_ext);

	% Insert file name parts into global arrays
	[name_F, name_F_ind(f)] = mssina(name_F, temp_name);
	name_F_len(f,1) = length(temp_name);
	[extension_F, extension_F_ind(f)] = mssina(extension_F, temp_ext);

	% Save model name
	[model_name, model_name_ind(f)] = mssina(model_name, models(model_ind(f)).name);
	model_name_len(f) = length(models(model_ind(f)).name);

	if models(model_ind(f)).par_num > max_par_num
		max_par_num(1) = models(model_ind(f)).par_num;
	end
	par_num(f) = models(model_ind(f)).par_num;

	% Load data from data file
	if strcmpi(temp_ext, '.csv') || strcmpi(temp_ext, '.txt')
		indata = csvread(full_path(full_path_ind(f),:));
	elseif strcmpi(temp_ext, '.xls') || strcmpi(temp_ext, '.xlsx')
		indata = xlsread(full_path(full_path_ind(f),:));
	else
		error(['Datafile has unknown format: ' full_path(full_path_ind(f),:)]);
	end
	clear temp_name temp_ext temp_ext_ind;

	% Get data length
	data_len(f) = size(indata, 1);

	% Get number of traces
	size_F(f) = size(indata, 2) - 1;

	% Ensure that data does not contain negative values
	if models(model_ind(f)).normalize_offset
		data_offset = min(min( indata(:,2:end) ));
		if data_offset < 0
			indata(:,2:end) = indata(:,2:end) - data_offset;
		end
	end

	% Compare lengths of data array and current data
	d_len = size(indata,1) - size(t,1);
	if d_len < 0
		fill_new = -d_len;
		fill_old = 0;
	else
		fill_new = 0;
		fill_old = d_len;
	end
	clear d_len;

	% Convert to desired time unit (required: string t_unit for plots)
	indata(:,1) = scale_t(indata(:,1), t_scale, f, t_offset);

	% Test if current time vector already contained in t
	for i = 1:size(t, 2)
		if size(indata, 1) <= size(t, 1) && isequal(indata(:,1), t(1:data_len(f),i))
			t_ind(f) = i;
		end
	end

	% If current time vector not yet in t, append it
	if t_ind(f) == 0
		% Concatenate t and current time vector
		t = [ [t; zeros(fill_old, size(t,2))] [indata(:,1); zeros(fill_new, 1)] ];
		t_ind(f) = size(t, 2);
	end

	clear fill_old fill_new;

	% Create time vector for simulation with resolution 0.1h
	t_sim_res(f) = .1 * min(diff( t(1:data_len(f), t_ind(f)) ));
	temp_t_sim = ( t(1, t_ind(f)) : t_sim_res(f) : t(data_len(f), t_ind(f)) )';
	data_sim_len(f) = length(temp_t_sim);

	% Test if simulation time vector already contained in t_sim
	for i = 1:size(t_sim, 2)
		if data_sim_len(f) <= size(t_sim, 1) && isequal(temp_t_sim, t_sim(1:data_sim_len(f),i))
			t_sim_ind(f) = i;
		end
	end

	% If time vector for simulation not yet in t_sim, append it
	if t_sim_ind(f) == 0
		if size(t_sim, 2) == 0
			% If there is no t_sim yet, copy t_sim_ind array
			t_sim = temp_t_sim;
		else
			% Build new t_sim array step by step
			% (to avoid warnings caused by concatenating empty arrays)
			fill_t_sim = length(temp_t_sim) - size(t_sim, 1);

			if fill_t_sim > 0
				temp_t_sim_old = [ t_sim; zeros(fill_t_sim, size(t_sim, 2)) ];
			else
				temp_t_sim_old = t_sim;
			end

			if fill_t_sim < 0
				temp_t_sim_new = [ temp_t_sim; zeros(-fill_t_sim, 1) ];
			else
				temp_t_sim_new = temp_t_sim;
			end

			t_sim = [ temp_t_sim_old temp_t_sim_new ];
			clear fill_t_sim temp_t_sim_old temp_t_sim_new;
		end
		t_sim_ind(f) = size(t_sim, 2);
	end
	clear temp_t_sim;

	% Prepare target directory
	if length(out_dir) > 1
		f_out = f;
	else
		f_out = 1;
	end
	[target_dir, target_dir_ind(f)] = mssina(target_dir, out_dir{f_out});
	target_dir_len(f) = length(out_dir{f_out});
	if ~exist(target_dir(target_dir_ind(f),1:target_dir_len(f)), 'dir')
	   [~,~,~] = mkdir(target_dir(target_dir_ind(f),1:target_dir_len(f)));
	end

	% Prepare format string for cell number in files
	[format_cell_number, format_cell_number_ind(f)] = mssina( ...
		format_cell_number, [ '%0' num2str( floor(log10(single(size_F(f)))) + 1 ) 'd' ]);

	% Allocate memory for array of trace indices
	if size(dataindices, 2) < size_F(f)
		dataindices = [ dataindices zeros(ndatafiles, size_F(f)-size(dataindices, 2)) ];
	end

	% Allocate memory for trace in global array of traces
	if size(data, 1) < data_len(f)
		fill_v = data_len(f) - size(data, 1);
	else
		fill_v = 0;
	end

	if size(data, 2) < ( ntraces + size_F(f) )
		fill_h = ntraces + size_F(f) - size(data, 2);
	else
		fill_h = 0;
	end

	data = [ data, zeros(size(data, 1), fill_h); zeros(fill_v, size(data, 2) + fill_h) ];
	clear fill_h fill_v;

	% Populate array of file indices
	file_ind(ntraces+1:ntraces+size_F(f),1) = f;

	% Populate array of file-wide trace indices
	index_F(ntraces+1:ntraces+size_F(f),1) = 0;

	% Populate global array of traces
	for i = 1:size_F(f)
		% Get index of trace in source data file
		index_F(ntraces + i, 1) = i;

		% Get index of trace in global array
		glob_ind = ntraces + index_F(ntraces + i, 1);

		% Write data to global array
		data(1:data_len(f),glob_ind) = indata(:,i+1);

		% Write global trace index to array of dataindices
		dataindices(f,i) = glob_ind;
	end

	% Update size of global array of traces and global number of traces
	ntraces = ntraces + size_F(f);
end

% Get components for correct log messages
if ntraces == 1
	log_trace = 'trace';
else
	log_trace = 'traces';
end

if ndatafiles == 1
	log_file = 'file';
else
	log_file = 'files';
end

%% Create temporary matfile (needed for parallel computing)
% Build name of global matfile
Ffile = fullfile(temp_dir, [ 'TEMP_' time_now '.mat' ]);

% Save input data in temporary matfile
save(Ffile, ...
	'full_path', 'name_F', 'extension_F', 'target_dir', 'model_name', ...
	'full_path_ind', 'name_F_ind', 'name_F_len', 'extension_F_ind', 'target_dir_ind', 'target_dir_len', ...
	'model_name_ind', 'model_name_len', 't', 't_sim', 't_unit', 'format_cell_number', ...
	'data_len', 'data_sim_len', 'size_F', 'dataindices', 'model_ind', 'models', ...
	't_ind', 't_sim_ind', 't_sim_res', 'format_cell_number_ind', 'par_num', ...
	'data', 'file_ind', 'index_F', 'has_interactivity', ...
	'outmode', 'time_now', 'temp_dir', ...
	'ntraces', 'ndatafiles', 'max_par_num', 'log_trace', 'log_file', ...
	'-v7.3');

% Clear file handle to temporary matfile
clear mf;

end