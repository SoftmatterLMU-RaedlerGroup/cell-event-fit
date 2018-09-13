function model = model_definitions(modelname)
%model_definitions provides the model definitions for `fit_apoptosis`
%
% Input:	The name of the model as a char vector
%
% Output:	A struct with the model properties with these fields:
%			name:		The internal name of the model
%			marker:		The human readable marker name
%			simulate:	The simulation function
%			postproc:	The postprocessing routine
%			par_num:	The number of parameters
%			par_min:	Vector of the minimum values the parameters can take
%			par_max:	Vector of the maximum values the parameters can take
%			par_names:	Names of the parameters (used on MS plots)
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

%% Allocate model structure
model = struct('name',[], 'marker',[], 'simulate',[], 'postproc',[], ...
	'par_num',0, 'par_min',[], 'par_max',[], 'par_names',[]);

%% Populate model structure

if regexpi(modelname, '^casp(?:a(?:s(?:e)?)?)?')
	% Caspase (recommended model)
	model.name = 'Caspase_OneSigmoidDecay';
	model.marker = 'Caspase 3/7';
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-2,-4];
	model.par_max = [+4,+4,+4,+4,+4,0,+4];
	model.par_names = {'A';'B';'t_0';'\alpha';'\gamma';'b';'\sigma'};

elseif strcmpi(modelname, 'pSIVA')
	% pSIVA (recommended routine for pSIVA)
	model.name = 'pSIVA';
	model.marker = 'pSIVA';
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-2,-4];
	model.par_max = [+4,+4,+4,+4,+4,+3,+4];
	model.par_names = {'A';'B';'t_0';'\alpha';'\gamma';'b';'\sigma'};

elseif strcmpi(modelname, 'pi')
	% PI (recommended routine for PI)
	model.name = 'PI';
	model.marker = 'PI';
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-2,-2,-2,-4];
	model.par_max = [+4,+4,+4,+1,+2,+3,+4];
	model.par_names = {'A';'B';'t_0';'\alpha';'\gamma';'b';'\sigma'};

elseif regexpi(modelname, '^(:?cell ?)?ro(?:s|x)$')
	% CellROX
	model.name = 'ROS_Parabola';
	model.marker = 'CellROX';
	model.simulate = @negative_parabola_sigmoid_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'a_\text{const}';'a_\text{quad}';'t_\text{vertex}'; ...
		't_\text{step}';'\gamma';'b';'\sigma'};

elseif regexpi(modelname, '^Lyso(:?Tracker)?$')
	% LysoTracker (recommended)
	model.name = 'Lyso';
	model.marker = 'LysoTracker';
	model.simulate = @negative_parabola_sigmoid_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'a_\text{const}';'a_\text{quad}';'t_\text{vertex}'; ...
		't_\text{step}';'\gamma';'b';'\sigma'};

elseif regexpi(modelname, '^TMRM(?:_var(?:iable)?)?$')
	% TMRM (parabola with variable sign; recommended routine for TMRM)
	model.name = 'TMRM_var';
	model.marker = 'TMRM';
	model.simulate = @parabola_sigmoid_simulate;
	model.postproc = @postproc_tmrm_kink;
	model.par_num = 8;
	model.par_min = [-6,-4,-1,-4,-4,-2,-4,-4];
	model.par_max = [+4,+4,+1,+4,+4,+2,+4,+4];
	model.par_names = {'a_\text{const}';'|a_\text{quad}|';'sign(a_\text{quad})'; ...
		't_\text{vertex}';'t_\text{step}';'\gamma';'b';'\sigma'};

elseif strcmpi(modelname, 'mRNA')
	% mRNA expression, "trivial model"
	model.name = 'mRNA_trivial';
	model.marker = 'eGFP';
	model.simulate = @mRNA_trivial_model_simulate;
	model.postproc = @postproc_mRNA_trivial;
	model.par_num = 6;
	model.par_min = [-4,-6,-6,-6,-6,-6];
	model.par_max = [+4,+6,+6,+6,+6,+6];
	model.par_names = {'\text{os}';'\text{sc}';'\delta';'\beta';'t_0';'\sigma'};
	
elseif regexpi(modelname, '^cal(?:cium|bryte)?(?:520)?')
	% Calcium marker (recommended model)
	model.name = 'Calcium';
	model.marker = 'Calbryte520';
	model.simulate = @LATE_decay_simulate;
	model.postproc = @postproc_LATE_decay;
	model.par_num = 10;
	model.par_min = [-3,-2,-3,-2,-1,-3,-3,-2,-2,-3];
	model.par_max = [+3,+4,+4,+4,+2,+3,+3,+2,+2,+3];
	model.par_names = {'C_0';'C_0^\text{amp}';'A';'B';'t_0';'a_1';'a_2';'b_1';'b_2';'\sigma'};

%{
elseif strcmpi(modelname, 'pSIVA_old')
	% pSIVA (legacy; made for uncorrected data)
	model.name = 'pSIVA_old';
	model.marker = 'pSIVA';
	model.simulate = @psiva_old_simulate;
	model.par_num = 9;
	model.par_min = [0,1,-1,-2,-2,-1,-2,-2,-4];
	model.par_max = [+3,+3,+2,+2,+3,+2,+2,+3,+4];
	model.par_names = {'A_{offset}';'B_{scale}';'C_{jump1}'; ...
		'D_{steep1}';'E_{decay1}';'F_{shift1}';'G_{start}';'H_{ramp}';'\sigma'};

elseif regexpi(modelname, '^casp(?:a(?:s(?:e)?)?)?_old')
	% Caspase (legacy; does not implement step in asymptotic behaviour)
	model.name = 'Caspase_VarSigmoidDecay';
	model.marker = 'Caspase 3/7';
	model.simulate = @casp_varsigmoid_decay_simulate;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'A_{offset}';'B_{difference}';'C_{ascent}'; ...
		'D_{steepness}';'E_{height}';'F_{decay}';'\sigma'};

elseif strcmpi(modelname, 'pi_old')
	% PI (legacy; threshold-based post processing)
	model.name = 'PI_Sigmoid2Decay';
	model.marker = 'PI';
	model.simulate = @ALL_sigmoid_decay_simulate;
	model.postproc = @postproc_ALL;
	model.par_num = 7;
	model.par_min = [-4,-3,-3,-2,-3,-3,-4];
	model.par_max = [+4,+4,+4,+2,+4,+4,+4];
	model.par_names = {'A_{offset}';'B_{scale}';'C_{jump}';'D_{steep}'; ...
		'E_{decay}';'F_{shift}';'\sigma'};

elseif strcmpi(modelname, 'pi_deriv2')
	% PI (experimental; not for productive use)
	model.name = 'PI_deriv2';
	model.marker = 'PI';
	model.simulate = @pi_sigmoid2_decay_simulate;
	model.par_num = 7;
	model.par_min = [-4,-3,-3,-2,-3,-3,-4];
	model.par_max = [+4,+4,+4,+2,+4,+4,+4];
	model.par_names = {'A_{offset}';'B_{scale}';'C_{jump}';'D_{steep}'; ...
		'E_{decay}';'F_{shift}';'\sigma'};
%}

elseif regexpi(modelname, '^lyso(?:tracker)?1$')
	% LysoTracker (experimental)
	model.name = 'Lyso1';
	model.marker = 'LysoTracker';
	model.simulate = @lyso_1_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 10;
	model.par_min = [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'t_0';'a_0';'t_1';'b_0';'b_1';'-b_2';'t_2';'c';'d';'\sigma'};

elseif regexpi(modelname, '^lyso(?:tracker)?2$')
	% LysoTracker (experimental)
	model.name = 'Lyso2';
	model.marker = 'LysoTracker';
	model.simulate = @lyso_2_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 8;
	model.par_min = [-4,-4,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'t_0';'b_0';'b_1';'-b_2';'t_1';'c';'d';'\sigma'};

elseif regexpi(modelname, '^lyso(?:tracker)?3$')
	% LysoTracker (experimental)
	model.name = 'Lyso3';
	model.marker = 'LysoTracker';
	model.simulate = @lyso_3_simulate;
	model.postproc = @postproc_kink;
	model.par_num = 7;
	model.par_min = [-4,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'t_0';'b_0';'-b_2';'t_1';'c';'d';'\sigma'};

elseif regexpi(modelname, '^TMRM_neg(?:ative)?$')
	% TMRM (negative parabola)
	model.name = 'TMRM_neg';
	model.marker = 'TMRM';
	model.simulate = @negative_parabola_sigmoid_simulate;
	model.postproc = @postproc_tmrm_kink;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'a_0';'-a_2';'t_\text{vertex}';'t_\text{step}'; ...
		'\alpha_\text{steep}';'d_\text{end}';'\sigma'};

elseif regexpi(modelname, '^TMRM_pos(?:itive)?$')
	% TMRM (positive parabola)
	model.name = 'TMRM_pos';
	model.marker = 'TMRM';
	model.simulate = @positive_parabola_sigmoid_simulate;
	model.postproc = @postproc_tmrm_kink;
	model.par_num = 7;
	model.par_min = [-6,-4,-4,-4,-4,-4,-4];
	model.par_max = [+4,+4,+4,+4,+4,+4,+4];
	model.par_names = {'a_0';'a_2';'t_\text{vertex}';'t_\text{step}'; ...
		'\alpha_\text{steep}';'d_\text{end}';'\sigma'};

elseif strcmpi(modelname, 'pSIVA_ALL')
	% pSIVA (deprecated)
	model.name = 'pSIVA';
	model.marker = 'pSIVA';
	model.simulate = @ALL_sigmoid_decay_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-4,-3,-3,-2,-3,-3,-4];
	model.par_max = [+4,+4,+4,+2,+4,+4,+4];
	model.par_names = {'A';'B';'t_\text{step}'; '\alpha';'\gamma'; ...
		't_\text{decay}';'\sigma'};

elseif regexpi(modelname, 'pi_?extrapol')
	% PI (deprecated)
	model.name = 'PI_extrapol';
	model.marker = 'PI';
	model.simulate = @LATE_simulate;
	model.postproc = @postproc_LATE_extrapol;
	model.par_num = 7;
	model.par_min = [-4,-3,-3,-2,-3,-3,-4];
	model.par_max = [+4,+4,+4,+2,+4,+4,+4];
	model.par_names = {'A';'B';'t_\text{step}';'\alpha';'\gamma'; ...
		't_\text{decay}';'\sigma'};

else
	% No model definition found for given model
	error(['There is no model "' modelname '".']);
end
