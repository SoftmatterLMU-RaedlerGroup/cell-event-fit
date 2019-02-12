% Get and read file
[filename,path] = uigetfile({'*.csv;*.txt';'*.*'}, 'Select file', fullfile('project', 'data'));
if filename == 0
	disp('Stop execution.')
	return
end
filepath = fullfile(path,filename);

D = csvread(filepath);
t = (D(:,1) - 1) / 6;
D = D(:,2:end);
n_traces = size(D, 2);

% For each trace, guess initial parameters
model = model_definitions('tmrm_new');
params = zeros(n_traces, model.par_num);

for iT = 1:n_traces
	S = doubleparabola_TMRM_guess(D(:,iT), t, model);
	params(iT,:) = S.guess;
end

% Plot each trace and its guess
for iT = 1:n_traces
	ax = gca(figure);
	hold(ax, 'on');

	d = D(:,iT);
	m = doubleparabola_TMRM_simulate(t, params(iT,:));

	plot(ax, t, d);
	plot(ax, t, m);
	plot(ax, t, parabola(t, params(iT, 1), params(iT, 2), params(iT, 3)), '--');
	plot(ax, t, parabola(t, params(iT, 6), params(iT, 7), params(iT, 8)), ':');

	ylim(ax, [min([d;m]), max([d;m])])
	xlabel(ax, 'Time [h]')
	ylabel(ax, 'Fluorescence [a.u.]')
	title(ax, sprintf('Trace %02d', iT));
	legend(ax, ["raw", "model", "parabola 1", "parabola 2"]);

end
clear ax d m


function y = parabola(t, a, x0, y0)
	y = a * (t - x0).^2 + y0;
end