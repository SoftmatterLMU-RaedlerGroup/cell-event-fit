mf = matfile('debug_inter.mat');

if ~exist('private', 'var')
	private = 0;
end

[is_to_fit, S_inter, data_indices, par_fun, private] = ...
	parabola_TMRM_interactive(mf.idx, mf.t, mf.data, mf.model, mf.S, mf.mf, private);