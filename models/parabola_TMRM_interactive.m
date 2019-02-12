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

if isempty(private)
	private = 0;
else
	private = private + 1;
end

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
% M = matfile('debug_inter.mat');
% M.idx = idx;
% M.t = t;
% M.data = data;
% M.model = model;
% M.S = S;
% M.mf = struct('index_F', mf.index_F);

index_F = mf.index_F(idx, 1);

t_event = S.params(4);
t_start = S.params(7);
t_end = S.params(8);


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
current_state = STATE_OUTSIDE;

fig = figure('MenuBar', 'none', ...
	'NumberTitle', 'off', ...
	'Name', 'Fitting result inspection', ...
	...'WindowStyle', 'modal', ...
	'Units', 'centimeters' ...
	);
fig.InnerPosition(3:4) = [20 + 1.1 * BUTTON_WIDTH, 15];

ax = axes(fig, 'Units', 'centimeters');
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
fig.WindowButtonMotionFcn = @move_mouse;
fig.WindowButtonDownFcn = @mouse_down;
fig.WindowButtonUpFcn = @mouse_release;

resize_figure();

plot(ax, t, data, t, model.simulate(t, S.params));
hold(ax, 'on');
p_evt = plot(ax, [t_event t_event], ax.YLim, '-', 'Color', [0 .5 0], 'LineWidth', 1);
p_before = rectangle('Parent', ax, 'FaceColor', [.2 .2 .2 .3], 'EdgeColor', 'none', 'Visible', 'off');
p_after = rectangle('Parent', ax, 'FaceColor', [.2 .2 .2 .3], 'EdgeColor', 'none', 'Visible', 'off');

xlabel(ax, 'Time [h]');
ylabel(ax, 'Flurescence [a.u.]');
title(ax, sprintf('Trace %d [Model: %s]', index_F, model.name));

text(ax, mean(ax.XLim), ax.YLim(2), sprintf('private=%d', private), ...
	'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

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
			S_inter.params = NaN(1, model.par_num, 'single');
			S_inter.logPost = NaN('single');
			S_inter.fit_type = -5;
		else
			S_inter.fit_type = 5;
		end

	case ACTION_REFIT
		is_to_fit = true;
		params = S.params;
		if isfinite(t_event)
			params(2) = t_event;
			params(4) = t_event;
		end
		params(7) = t_start;
		params(8) = t_end;

		par_min = model.par_min;
		par_min(7) = t_start;
		par_min(8) = t_end;

		par_max = model.par_max;
		par_max(7) = t_start;
		par_max(8) = t_end;
		par_fun = @(~,~,~) struct( ...
			'guess', double(params), ...
			'min', par_min, ...
			'max', par_max);
		data_indices = find((t >= t_start) & (t <= t_end));
		disp(data_indices)

	case ACTION_DISCARD
		is_to_fit = false;
		S_inter.t_event = NaN('single');
		if have_results_changed
			S_inter.params = NaN(1, model.par_num, 'single');
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
		delta = 0.01 * diff(ax.XLim);

		% Check mouse position relative to limits
		if mx >= ax.XLim(1) && mx <= ax.XLim(2) && my >= ax.YLim(1) && my <= ax.YLim(2)
			%fprintf(1, 'x=%f, y=%f\n', mx, my); % DEBUG
			is_in_ax = true;
		else
			is_in_ax = false;
		end
		if is_in_ax && (~isfinite(t_start) || mx > t_start) && (~isfinite(t_end) || mx < t_end)
			is_between_limits = true;
		else
			is_between_limits = false;
		end

		% State dependent actions
		if current_state == STATE_EVENT
			if is_between_limits
				t_event = mx;
			else
				t_event = NaN;
			end
			p_evt.XData = [t_event, t_event];
			check_param_change();

		elseif current_state == STATE_START
			if isfinite(t_event) && (mx >= t_event - delta)
				t_start = t_event - delta;
			elseif isfinite(t_end) && (mx >= t_end - delta)
				t_start = t_end - delta;
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
			check_param_change();

		elseif current_state == STATE_END
			if isfinite(t_event) && (mx <= t_event + delta)
				t_end = t_event + delta;
			elseif isfinite(t_start) && (mx <= t_start + delta)
				t_end = t_start + delta;
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
			check_param_change();

		elseif is_between_limits && (~isfinite(t_event) || abs(mx - t_event) < delta)
			current_state = STATE_PRE_EVENT;

		elseif is_in_ax && ( (isfinite(t_start) && abs(mx - t_start) < delta) || (abs(mx - ax.XLim(1)) < delta) )
			current_state = STATE_PRE_START;

		elseif is_in_ax && ( (isfinite(t_end) && abs(mx - t_end) < delta) || (abs(mx - ax.XLim(2)) < delta) )
			current_state = STATE_PRE_END;

		elseif is_between_limits
			current_state = STATE_INSIDE;

		else
			current_state = STATE_OUTSIDE;
		end
		update_pointer()
	end


	function mouse_down(~,~,~)
		% Callback for mouse down events
		if current_state == STATE_PRE_EVENT
			current_state = STATE_EVENT;
			move_mouse();
			return
		elseif current_state == STATE_PRE_START
			current_state = STATE_START;
		elseif current_state == STATE_PRE_END
			current_state = STATE_END;
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
		end
		move_mouse();
	end


	function update_pointer()
		% Set mouse pointer appearance based on current state
		switch current_state
			case STATE_EVENT
				fig.Pointer = 'fleur';
			case STATE_PRE_EVENT
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
			otherwise
				fig.Pointer = 'arrow';
		end
	end


	function select_action(src, ~, ~)
		% Callback for action buttons
		selected_action = src.UserData;
		close(fig);
	end


	function check_param_change()
		% Check if parameters have been changed
		new_change = t_event ~= S.params(4) && ~(isnan(t_event) && isnan(S.params(4)));
		new_change = new_change || ...
			(t_start ~= S.params(7) && ~(isnan(t_start) && isnan(S.params(7))));
		new_change = new_change || ...
			(t_end ~= S.params(8) && ~(isnan(t_end) && isnan(S.params(8))));

		if new_change ~= have_results_changed
			have_results_changed = new_change;
			if new_change
				btn_refit.Enable = 'on';
			else
				btn_refit.Enable = 'off';
			end
		end
	end

end
