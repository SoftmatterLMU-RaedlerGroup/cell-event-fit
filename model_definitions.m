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
%			preproc:	A preprocessing routine (optional)
%			postproc:	The postprocessing routine
%			interactive: (optional) function handle for interactivity
%			n_starts:	The number of fitting runs (-1 for default)
%			par_num:	The number of parameters
%			par_min:	Vector of the minimum values the parameters can take
%			par_max:	Vector of the maximum values the parameters can take
%			weights:	Weights of the data points
%			par_log:	Logical or index vector for logarithmic parameters
%			guess:		Set of initial parameter values (optional)
%			par_fun:	Function handle for setting `par_min`, `par_max`, `guess`
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

%% Initialize model structure
model = struct('name',[], 'marker',[], 'normalize_offset',false, ...
	'simulate',[], 'preproc',false, 'postproc',false, 'interactive',[], ...
	'n_starts',-1, 'par_num',0, 'par_min',[], 'par_max',[], 'weights',[], ...
	'par_log',false, 'guess',[], 'par_fun',[], 'par_names',[]);

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

elseif strcmp(modelname, 'Test')
	% Test (linear function)
	model.name = 'Test_linear';
	model.marker = 'line';
	model.simulate = @(t, p) p(1)*t+p(2);
	model.par_num = 2;
	model.par_min = [-100,-100];
	model.par_max = [+100,+100];
	model.par_log = false;
	model.weights = 1;
	model.par_names = {'a','b'};

elseif regexpi(modelname, '^tmrm_?i(?:nter(?:active)?)?$')
	% New TMRM model
	model.name = 'TMRM interactive';
	model.marker = 'TMRM';
	model.simulate = @parabola_TMRM_simulate;
	model.interactive = @parabola_TMRM_interactive;
	model.n_starts = 1;
	model.par_num = 8;
	model.par_min = [-100,-10,-1000, 0,.01,-1000,0,30];
	model.par_max = [ 100, 40, 1000,30, 10,+1000,0,30];
	model.par_log = false;
	model.par_names = {'a_1','x_0','y_0','\tau','a_2','y_\mathrm{end}','t_\mathrm{start}','t_\mathrm{end}'};

elseif regexpi(modelname, '^tmrm_?n(?:ew?)$')
	% New TMRM model
	model.name = 'TMRM new';
	model.marker = 'TMRM';
	model.simulate = @doubleparabola_TMRM_simulate;
	model.par_fun = @doubleparabola_TMRM_guess;
	model.par_num = 11;
	model.par_min = [.1, 0,-1000, 0,-3,.1, 0,-1000, 0,-3,-1000];
	model.par_max = [10,30,+1000,30,+3,10,30,+1000,30,+3,+1000];
	model.par_log = [false,false,false,false,true,false,false,false,false,true,false];
	model.par_names = {'a_1','x_1','y_1','t_1','s_1','a_2','x_2','y_2','t_2','a_3','y_3'};

else
	% No model definition found for given model
	error(['There is no model "' modelname '".']);
end
