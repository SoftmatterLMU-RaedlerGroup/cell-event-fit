function model = model_definitions(modelname)
%model_definitions provides the model definitions for `fit_apoptosis`
%
% Input:	The name of the model as a char vector
%
% Output:	A struct with the model properties with these fields:
%			name:		The internal name of the model
%			marker:		The human readable marker name
%			normalize_offset: Logical indicator for offset normalization
%			simulate:	The simulation function
%			postproc:	The postprocessing routine
%			preproc:	A preprocessing routine (optional)
%			par_num:	The number of parameters
%			par_min:	Vector of the minimum values the parameters can take
%			par_max:	Vector of the maximum values the parameters can take
%			weights:	Weights of the data points
%			par_log:	Logical or index vector for logarithmic parameters
%			par_names:	Names of the parameters (used on MS plots)
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

%% Allocate model structure
model = struct('name',[], 'marker',[], 'normalize_offset',false, ...
	'simulate',[], 'postproc',false, 'preproc',false, ...
	'par_num',0, 'par_min',[], 'par_max',[], 'weights',[], ...
	'par_log',false, 'par_names',[]);

%% Populate model structure

if regexpi(modelname, '^casp(?:a(?:s(?:e)?)?)?')
	% Caspase
	model.name = 'Caspase_OneSigmoidDecay';
	model.marker = 'Caspase 3/7';
	model.normalize_offset = true;
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-2,-4];
	model.par_max = [+4,+4,+4,+4,+4,0,+4];
	model.par_log = true;
	model.weights = struct('parameter', 7);
	model.par_names = {'A';'B';'t_0';'\alpha';'\gamma';'b';'\sigma'};

elseif strcmpi(modelname, 'pSIVA')
	% pSIVA
	model.name = 'pSIVA';
	model.marker = 'pSIVA';
	model.normalize_offset = true;
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.preproc = @preproc_psiva;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-2,-4];
	model.par_max = [+4,+4,+4,+4,+4,+3,+4];
	model.par_log = true;
	model.weights = struct('parameter', 7);
	model.par_names = {'A';'B';'t_0';'\alpha';'\gamma';'b';'\sigma'};

elseif strcmpi(modelname, 'pi')
	% PI
	model.name = 'PI';
	model.marker = 'PI';
	model.normalize_offset = true;
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-2,-2,-2,-4];
	model.par_max = [+4,+4,+4,+1,+2,+3,+4];
	model.par_log = true;
	model.weights = struct('parameter', 7);
	model.par_names = {'A';'B';'t_0';'\alpha';'\gamma';'b';'\sigma'};

elseif regexpi(modelname, '^(:?cell ?)?ro(?:s|x)$')
	% CellROX
	model.name = 'ROS_Parabola';
	model.marker = 'CellROX';
	model.normalize_offset = true;
	model.simulate = @negative_parabola_sigmoid_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_log = true;
	model.weights = struct('parameter', 7);
	model.par_names = {'a_\text{const}';'a_\text{quad}';'t_\text{vertex}'; ...
		't_\text{step}';'\gamma';'b';'\sigma'};

elseif regexpi(modelname, '^Lyso(:?Tracker)?$')
	% LysoTracker
	model.name = 'Lyso';
	model.marker = 'LysoTracker';
	model.normalize_offset = true;
	model.simulate = @negative_parabola_sigmoid_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_log = true;
	model.weights = struct('parameter', 7);
	model.par_names = {'a_\text{const}';'a_\text{quad}';'t_\text{vertex}'; ...
		't_\text{step}';'\gamma';'b';'\sigma'};

elseif regexpi(modelname, '^TMRM(?:_var(?:iable)?)?$')
	% TMRM (parabola with variable sign)
	model.name = 'TMRM_var';
	model.marker = 'TMRM';
	model.normalize_offset = true;
	model.simulate = @parabola_sigmoid_simulate;
	model.postproc = @postproc_tmrm_kink;
	model.par_num = 8;
	model.par_min = [-6,-4,-1,-4,-4,-2,-4,-4];
	model.par_max = [+4,+4,+1,+4,+4,+2,+4,+4];
	model.par_log = true;
	model.weights = struct('parameter', 8);
	model.par_names = {'a_\text{const}';'|a_\text{quad}|';'sign(a_\text{quad})'; ...
		't_\text{vertex}';'t_\text{step}';'\gamma';'b';'\sigma'};

elseif strcmpi(modelname, 'mRNA')
	% mRNA expression, "trivial model"
	model.name = 'mRNA_trivial';
	model.marker = 'eGFP';
	model.normalize_offset = true;
	model.simulate = @mRNA_trivial_model_simulate;
	model.postproc = @postproc_mRNA_trivial;
	model.par_num = 6;
	model.par_min = [-4,-6,-6,-6,-6,-6];
	model.par_max = [+4,+6,+6,+6,+6,+6];
	model.par_log = true;
	model.weights = struct('parameter', 6);
	model.par_names = {'\text{os}';'\text{sc}';'\delta';'\beta';'t_0';'\sigma'};
	
elseif regexpi(modelname, '^cal(?:cium)?(?:520)?')
	% Calcium marker
	model.name = 'Calcium';
	model.marker = 'Cal520';
	model.normalize_offset = true;
	model.simulate = @LATE_decay_simulate;
	model.postproc = @postproc_LATE_decay;
	model.preproc = @preproc_cal520;
	model.par_num = 10;
	model.par_min = [-3,-2,-3,-2,-1,-3,-3,-2,-2,-3];
	model.par_max = [+3,+4,+4,+4,+2,+3,+3,+2,+2,+3];
	model.par_log = true;
	model.weights = struct('parameter', 10);
	model.par_names = {'C_0';'C_0^\text{amp}';'A';'B';'t_0';'a_1';'a_2';'b_1';'b_2';'\sigma'};

else
	% No model definition found for given model
	error(['There is no model "' modelname '".']);
end
