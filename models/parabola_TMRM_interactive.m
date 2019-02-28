function [is_to_fit, S_inter, data_indices, par_fun, private] = ...
	parabola_TMRM_interactive(idx, t, data, model, S, mf, private)
%PARABOLA_TMRM_INTERACTIVE interactive function for TMRM fitting
%
% Arguments:
%	idx			index of current trace (in matfile)
%	t			time vector of raw data
%	data		raw data to be fitted
%	model		model struct
%	S			struct with results from fitting
%	mf			handle to global matfile
%	private		private data of this function specific for this trace
%
% Returns:
%	is_to_fit		boolean whether to fit or to use return values of this
%					function
%	S_inter			struct of return results (like `S`)
%	data_indices	indices of points in `data` and `t` for fitting
%	par_fun			parameter function (see `model.par_fun`)
%	private			private trace-specific data of this function to be saved
%
% Copyright © 2019 Daniel Woschée <daniel.woschee@physik.lmu.de>
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
data_indices = [];
par_fun = [];

if isempty(S)
	is_to_fit = true;
	S_inter = [];
	return
else
	is_to_fit = false;
	S_inter = struct;
end

ACTION_ACCEPT = 'accept';
ACTION_DISCARD = 'discard';
ACTION_REFIT = 'refit';
selected_action = ACTION_ACCEPT;
have_results_changed = false;

% DEBUG
% debug_mat = 'debug_inter.mat';
% if ~exist(debug_mat, 'file')
% 	M = matfile(debug_mat);
% 	M.idx = idx;
% 	M.t = t;
% 	M.data = data;
% 	M.model = model;
% 	M.S = S;
% 	M.mf = struct( ...
% 		'index_F', mf.index_F, ...
% 		'file_ind', mf.file_ind, ...
% 		'name_F', mf.name_F, ...
% 		'name_F_ind', mf.name_F_ind, ...
% 		'name_F_len', mf.name_F_len, ...
% 		'size_F', mf.size_F ...
% 		);
% end

index_F = mf.index_F(idx, 1);
file_ind = mf.file_ind(idx, 1);
filename = mf.name_F(mf.name_F_ind(file_ind,1), 1:mf.name_F_len(file_ind,1));
size_F = mf.size_F(file_ind, 1);

a_par = S.params(1);
x0 = S.params(2);
y0 = S.params(3);
t_break = S.params(4);
a_break = S.params(5);
y_end = S.params(6);
t_start = S.params(7);
t_end = S.params(8);

if isfield(S, 't_event')
	t_event = S.t_event;
else
	t_event = t_break;
end

par_min = model.par_min;
par_max = model.par_max;
if isa(model.par_fun, 'function_handle')
	par_props = model.par_fun(data, t, model);
	if isfield(par_props, 'min')
		par_min = par_props.min;
	end
	if isfield(par_props, 'max')
		par_max = par_props.max;
	end
end


BUTTON_WIDTH = 2.5;
BUTTON_HEIGHT = 1;

STATE_OUTSIDE = 0;
STATE_INSIDE = 1;
STATE_PRE_EVENT = 2;
STATE_EVENT = 3;
STATE_PRE_START = 4;
STATE_START = 5;
STATE_PRE_END = 6;
STATE_END = 7;
STATE_PRE_PAR_A = 8;
STATE_PAR_A = 9;
STATE_PRE_PAR_VTX = 10;
STATE_PAR_VTX = 11;
STATE_PRE_Y_END = 12;
STATE_Y_END = 13;
STATE_PRE_BREAK = 14;
STATE_BREAK = 15;
current_state = STATE_OUTSIDE;
mouse_down_info = struct();

fig = figure('MenuBar', 'none', ...
	'NumberTitle', 'off', ...
	'Name', 'Fitting result inspection', ...
	...'WindowStyle', 'modal', ...
	'CloseRequestFcn', @window_close_btn, ...
	'Units', 'centimeters' ...
	);
fig.InnerPosition(3:4) = [20 + 1.1 * BUTTON_WIDTH, 15];

ax = axes(fig, ...
	'Units', 'centimeters', ...
	'NextPlot', 'add', ...
	'XLimMode', 'manual', ...
	'YLimMode', 'manual');
btn_accept = uicontrol(fig, ...
	'Style', 'pushbutton', ...
	'Units', 'centimeters', ...
	'String', 'Accept', ...
	'Tooltip', 'Accept these values for the trace', ...
	'UserData', ACTION_ACCEPT, ...
	'Callback', @select_action);
btn_discard = uicontrol(fig, ...
	'Style', 'pushbutton', ...
	'Units', 'centimeters', ...
	'String', 'Discard', ...
	'Tooltip', 'Discard the trace; set event time to NaN', ...
	'UserData', ACTION_DISCARD, ...
	'Callback', @select_action);
btn_refit = uicontrol(fig, ...
	'Style', 'pushbutton', ...
	'Units', 'centimeters', ...
	'String', 'Fit again', ...
	'Tooltip', 'Fit the trace again', ...
	'UserData', ACTION_REFIT, ...
	'Enable', 'off', ...
	'Callback', @select_action);

fig.SizeChangedFcn = @resize_figure;
resize_figure();

ax.XLim = [min(t) max(t)];
data_amp = max(data) - min(data);
ax.YLim = [min(data) - .05 * data_amp, max(data) + .05 * data_amp];

plot(ax, t, data);
p_fit = plot(ax, t, NaN(size(t))); %plot(ax, t, model.simulate(t, S.params));
p_evt = plot(ax, [t_event t_event], ax.YLim, '-', 'Color', [0 .5 0], 'LineWidth', 1);
p_before = rectangle('Parent', ax, 'FaceColor', [.2 .2 .2 .3], 'EdgeColor', 'none', 'Visible', 'off');
p_after = rectangle('Parent', ax, 'FaceColor', [.2 .2 .2 .3], 'EdgeColor', 'none', 'Visible', 'off');

xlabel(ax, 'Time [h]');
ylabel(ax, 'Flurescence [a.u.]');
title(ax, sprintf('Model: %s\nFile: %s [Trace %0*d]', ...
	model.name, filename, floor(log10(double(size_F)))+1, index_F), ...
	'Interpreter', 'none');

% DEBUG: display derivative change at breakdown
p_info = text(ax, ax.XLim(2), ax.YLim(2), '', ...
	'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
	'FontSize', 10);

p_handle = plot(ax, NaN, NaN, 'Marker', 'o', 'MarkerSize', 6, ...
	'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');

current_t_begin = NaN;
current_y_break = NaN;
current_t_stop = NaN;
current_y_stop = NaN;
current_break_deriv = NaN;
current_deriv_change = NaN;
current_amplitude = NaN;
update_params()
update_info_msg();

fig.WindowButtonMotionFcn = @move_mouse;
fig.WindowButtonDownFcn = @mouse_down;
fig.WindowButtonUpFcn = @mouse_release;

% Display initial limits, if defined
if isfinite(t_start) && (t_start > ax.XLim(1))
	p_before.Position = [ax.XLim(1), ax.YLim(1), t_start - ax.XLim(1), diff(ax.YLim)];
	p_before.Visible = 'on';
end
if isfinite(t_end) && (t_end < ax.XLim(2))
	p_after.Position = [t_end, ax.YLim(1), ax.XLim(2) - t_end, diff(ax.YLim)];
	p_after.Visible = 'on';
end

% When Figure is closed, return results
waitfor(fig);
switch selected_action
	case ACTION_ACCEPT
		is_to_fit = false;
		S_inter.t_event = single(t_event);
		if have_results_changed
			S_inter.params = single(current_params());
			S_inter.logPost = NaN('single');
			S_inter.fit_type = -5;
			if isfinite(t_event)
				S_inter.event_deriv = current_break_deriv;
			else
				S_inter.event_deriv = NaN;
			end
			S_inter.custom_data_values(1) = current_deriv_change;
			S_inter.amplitude = current_amplitude;
		else
			S_inter.fit_type = 5;
		end

	case ACTION_REFIT
		is_to_fit = true;
		params = current_params();
		par_min(4) = max(t_start, min(t));
		par_min(7) = params(7);
		par_min(8) = params(8);
		par_max(4) = min(t_end, max(t));
		par_max(7) = params(7);
		par_max(8) = params(8);
		par_fun = @(~,~,~) struct( ...
			'guess', double(params), ...
			'min', par_min, ...
			'max', par_max, ...
			'n_starts', 5);
		data_indices = find((t >= t_start) & (t <= t_end));

	case ACTION_DISCARD
		is_to_fit = false;
		S_inter.t_event = NaN('single');
		if have_results_changed
			S_inter.params = single(current_params());
			S_inter.logPost = NaN('single');
			S_inter.fit_type = -5;
		else
			S_inter.fit_type = 5;
		end
end


	function resize_figure(~,~,~)
		% Callback for figure resize events
		fig_width = fig.InnerPosition(3);
		fig_height = fig.InnerPosition(4);

		ax.OuterPosition = [0, 0, fig_width - 1.1 * BUTTON_WIDTH, fig_height];
		btn_accept.Position = [fig_width - 1.05 * BUTTON_WIDTH, ...
			fig_height - 1.05 * BUTTON_HEIGHT, ...
			BUTTON_WIDTH, BUTTON_HEIGHT];
		btn_refit.Position = [fig_width - 1.05 * BUTTON_WIDTH, ...
			fig_height - 2.15 * BUTTON_HEIGHT, ...
			BUTTON_WIDTH, BUTTON_HEIGHT];
		btn_discard.Position = [fig_width - 1.05 * BUTTON_WIDTH, ...
			fig_height - 3.25 * BUTTON_HEIGHT, ...
			BUTTON_WIDTH, BUTTON_HEIGHT];
	end


	function move_mouse(~,~,~)
		% Callback for mouse motion events
		mx = ax.CurrentPoint(1, 1);
		my = ax.CurrentPoint(1, 2);
		deltax = 0.02 * diff(ax.XLim);
		deltay = 0.02 * diff(ax.YLim);

		% Check mouse position relative to limits
		if mx >= ax.XLim(1) && mx <= ax.XLim(2) && my >= ax.YLim(1) && my <= ax.YLim(2)
			is_in_ax = true;
		else
			is_in_ax = false;
		end
		if is_in_ax && (mx > current_t_begin) && (mx < current_t_stop)
			is_between_limits = true;
		else
			is_between_limits = false;
		end

		% State dependent actions
		if (current_state == STATE_EVENT) || (current_state == STATE_BREAK)
			if mx <= current_t_begin + deltax
				t_break = current_t_begin + deltax;
			elseif mx >= current_t_stop - deltax
				t_break = current_t_stop - deltax;
			else
				t_break = mx;
			end

			switch_event(is_in_ax);
			update_params();

		elseif current_state == STATE_START
			if mx >= t_break - deltax
				t_start = t_break - deltax;
			elseif mx >= current_t_stop - deltax
				t_start = current_t_stop - deltax;
			else
				t_start = mx;
			end
			if t_start <= ax.XLim(1)
				t_start = ax.XLim(1);
				p_before.Visible = 'off';
			else
				p_before.Position = [ax.XLim(1), ax.YLim(1), t_start - ax.XLim(1), diff(ax.YLim)];
				p_before.Visible = 'on';
			end
			update_params();

		elseif current_state == STATE_END
			if mx <= t_break + deltax
				t_end = t_break + deltax;
			elseif isfinite(t_start) && (mx <= t_start + deltax)
				t_end = t_start + deltax;
			else
				t_end = mx;
			end
			if t_end >= ax.XLim(2)
				t_end = ax.XLim(2);
				p_after.Visible = 'off';
			else
				p_after.Position = [t_end, ax.YLim(1), ax.XLim(2) - t_end, diff(ax.YLim)];
				p_after.Visible = 'on';
			end
			update_params();

		elseif current_state == STATE_Y_END
			dy = (my - mouse_down_info.anchor(2));
			y_end = mouse_down_info.y_end + dy;
			if y_end < par_min(6)
				y_end = par_min(6);
			elseif y_end > par_max(6)
				y_end = par_max(6);
			end
			ys = mouse_down_info.y_stop;
			denom = current_t_stop - t_break + (mx - mouse_down_info.anchor(1))^3;
			if denom < 0
				denom = 0;
			end
			a_break = - log(abs( (ys + dy - y_end) / (current_y_break - y_end) )) / denom;
			if a_break < par_min(5)
				a_break = par_min(5);
			elseif a_break > par_max(5)
				a_break = par_max(5);
			end
			update_params();

		elseif current_state == STATE_PAR_A
			a_par = (my - y0) / (mx - x0)^2;
			fit_at_begin = model.simulate(current_t_begin, current_params());
			if fit_at_begin > ax.YLim(2)
				a_par = (ax.YLim(2) - y0) / (current_t_begin - x0)^2;
			elseif fit_at_begin < ax.YLim(1)
				a_par = (ax.YLim(1) - y0) / (current_t_begin - x0)^2;
			end
			update_params();

		elseif current_state == STATE_PAR_VTX
			x0 = mouse_down_info.x0 + mx - mouse_down_info.anchor(1);
			if x0 < par_min(2)
				x0 = par_min(2);
			elseif x0 > par_max(2)
				x0 = par_max(2);
			end

			y0 = mouse_down_info.y0 + my - mouse_down_info.anchor(2);
			if y0 < par_min(3)
				y0 = par_min(3);
			elseif y0 > par_max(3)
				y0 = par_max(3);
			end
			update_params;

		elseif is_in_ax && (abs(current_y_stop - my) < deltay) && ...
				(mx >= current_t_stop - deltax)
			current_state = STATE_PRE_Y_END;

		elseif is_between_limits && (mx >= current_t_begin) && ...
				(mx <= t_break) && (abs(my - model.simulate(mx, current_params())) < deltay)
			if (abs(mx - x0) < deltax) || ...
					((x0 <= current_t_begin) && (mx < current_t_begin + deltax)) || ...
					((x0 >= t_break) && (mx > t_break - 3 * deltax))
				current_state = STATE_PRE_PAR_VTX;
			else
				current_state = STATE_PRE_PAR_A;
			end

		elseif is_between_limits && (abs(mx - t_break) < deltax) && (abs(my - current_y_break) < deltay)
				current_state = STATE_PRE_BREAK;

		elseif is_between_limits && (~isfinite(t_event) || abs(mx - t_event) < deltax)
			current_state = STATE_PRE_EVENT;

		elseif is_in_ax && (abs(mx - current_t_begin) < deltax)
			current_state = STATE_PRE_START;

		elseif is_in_ax && (abs(mx - current_t_stop) < deltax)
			current_state = STATE_PRE_END;

		elseif is_between_limits
			current_state = STATE_INSIDE;

		else
			current_state = STATE_OUTSIDE;
		end
		update_pointer();
		update_handle(mx, my);
	end


	function mouse_down(~,~,~)
		% Callback for mouse down events
		mouse_down_info.anchor = ax.CurrentPoint(1, 1:2);
		switch current_state
			case STATE_PRE_EVENT
				current_state = STATE_EVENT;
				switch_event();
				update_params();
			case STATE_PRE_START
				current_state = STATE_START;
			case STATE_PRE_END
				current_state = STATE_END;
			case STATE_PRE_Y_END
				current_state = STATE_Y_END;
				mouse_down_info.y_end = y_end;
				mouse_down_info.y_stop = current_y_stop;
			case STATE_PRE_BREAK
				current_state = STATE_BREAK;
				switch_event();
				update_params();
			case STATE_PRE_PAR_A
				current_state = STATE_PAR_A;
			case STATE_PRE_PAR_VTX
				current_state = STATE_PAR_VTX;
				mouse_down_info.x0 = x0;
				mouse_down_info.y0 = y0;
		end
		update_pointer()
	end


	function mouse_release(~,~,~)
		% Callback for mouse release events
		switch current_state
			case STATE_EVENT
				current_state = STATE_PRE_EVENT;
			case STATE_START
				current_state = STATE_PRE_START;
			case STATE_END
				current_state = STATE_PRE_END;
			case STATE_Y_END
				current_state = STATE_PRE_Y_END;
			case STATE_BREAK
				current_state = STATE_PRE_BREAK;
			case STATE_PAR_A
				current_state = STATE_PRE_PAR_A;
			case STATE_PAR_VTX
				current_state = STATE_PRE_PAR_VTX;
		end
		mouse_down_info = struct();
		move_mouse();
	end


	function update_pointer()
		% Set mouse pointer appearance based on current state
		switch current_state
			case {STATE_EVENT, STATE_BREAK, STATE_PAR_A, STATE_PAR_VTX}
				fig.Pointer = 'fleur';
			case {STATE_PRE_EVENT, STATE_PRE_BREAK, STATE_PRE_PAR_A, STATE_PRE_PAR_VTX}
				fig.Pointer = 'hand';
			case STATE_INSIDE
				fig.Pointer = 'arrow';
			case STATE_PRE_START
				fig.Pointer = 'left';
			case STATE_START
				fig.Pointer = 'left';
			case STATE_PRE_END
				fig.Pointer = 'right';
			case STATE_END
				fig.Pointer = 'right';
			case STATE_PRE_Y_END
				fig.Pointer = 'hand';
			case STATE_Y_END
				fig.Pointer = 'fleur';
			otherwise
				fig.Pointer = 'arrow';
		end
	end


	function update_handle(mx, my)
		% Show (or hide) mouse handle for parameter adjustment
		switch current_state
			case {STATE_PRE_EVENT, STATE_EVENT, STATE_PRE_BREAK, STATE_BREAK}
				p_handle.XData = t_break;
				p_handle.YData = current_y_break;
				p_handle.Marker = 'o';

			case {STATE_PRE_Y_END, STATE_Y_END}
				if mx < ax.XLim(2) && current_t_stop < ax.XLim(2)
					if mx < current_t_stop
						p_handle.XData = current_t_stop;
					else
						p_handle.XData = mx;
					end
				else
					p_handle.XData = ax.XLim(2);
				end
				p_handle.YData = current_y_stop;
				p_handle.Marker = 'o';

			case {STATE_PRE_PAR_A, STATE_PAR_A}
				p_handle.XData = mx;
				p_handle.YData = model.simulate(mx, current_params(), 'parabola');
				p_handle.Marker = 'o';

			case {STATE_PRE_PAR_VTX, STATE_PAR_VTX}
				p_handle.XData = mx;
				p_handle.YData = model.simulate(mx, current_params(), 'parabola');
				if x0 < current_t_begin
					p_handle.Marker = '<';
				elseif x0 > current_t_stop
					p_handle.Marker = '>';
				else
					p_handle.Marker = 'v';
				end

			otherwise
				p_handle.XData = NaN;
				p_handle.YData = NaN;
		end
	end


	function select_action(src, ~, ~)
		% Callback for action buttons
		selected_action = src.UserData;
		delete(fig);
	end


	function P = current_params()
		% Return array of current parameter values
		P = zeros(1, 8);
		P(1) = a_par;
		P(2) = x0;
		P(3) = y0;
		P(4) = t_break;
		P(5) = a_break;
		P(6) = y_end;
		P(7) = t_start;
		P(8) = t_end;
	end


	function update_params()
		% Refresh intermediate values and plot
		cur_par = current_params();
		current_t_begin = min(t);
		if isfinite(t_start) && t_start > current_t_begin
			current_t_begin = t_start;
		end
		current_y_break = model.simulate(t_break, cur_par);
		current_t_stop =  max(t);
		if isfinite(t_end) && t_end < current_t_stop
			current_t_stop = t_end;
		end
		current_y_stop = model.simulate(current_t_stop, cur_par);

		p_fit.YData = model.simulate(t, cur_par);

		current_amplitude = max(p_fit.YData) - min(p_fit.XData);
		current_break_deriv = -a_break * (current_y_break - y_end);
 		parab_deriv = 2 * a_par * (t_break - x0);
 		current_deriv_change = (parab_deriv - current_break_deriv) / current_amplitude;
		current_amplitude = max(p_fit.YData) - min(p_fit.YData);

		update_info_msg();
		check_param_change();
	end


	function check_param_change()
		% Check if parameters have been changed
		changed_params = S.params ~= current_params();
		has_changed = any(changed_params([ones(1, 6, 'logical') ~isnan(t_start) ~isnan(t_end)]));
		if ~has_changed && ( (t_event ~= S.t_event) && ~(isnan(t_event) && isnan(S.t_event)) )
			has_changed = true;
		end
		if ~has_changed && (isnan(t_start) && isnan(S.params(7)))
			has_changed = true;
		end
		if ~has_changed && (isnan(t_end) && isnan(S.params(8)))
			has_changed = true;
		end

		if has_changed ~= have_results_changed
			have_results_changed = has_changed;
			if has_changed
				btn_refit.Enable = 'on';
			else
				btn_refit.Enable = 'off';
			end
		end
	end

	function switch_event(new_state)
		% Switch event on or off
		if nargin < 1
			if isfinite(t_event)
				t_event = NaN;
			else
				t_event = t_break;
			end
		elseif new_state
			t_event = t_break;
		else
			t_event = NaN;
		end
		p_evt.XData = [t_event, t_event];
	end

	function update_info_msg()
		% Update info message in upper right corner
		p_info.String = sprintf('t_{event} = %.2f\n\\delta = %.2f\n\\Delta = %.2f', ...
			t_event, current_break_deriv/current_amplitude, current_deriv_change);
	end

	function window_close_btn(~,~,~)
		% Window close callback
		msgbox('Please choose an action to close this window.', ...
			'Interactive TMRM fit', 'warn', 'modal');
	end

end
