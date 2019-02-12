function S = doubleparabola_TMRM_guess(data, t, model, ~)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initialize parameter properties
% if ~isempty(model.par_min)
% 	par_min = model.par_min;
% else
% 	par_min = -inf(1, model.par_num);
% end
% if ~isempty(model.par_max)
% 	par_max = model.par_max;
% else
% 	par_max = inf(1, model.par_num);
% end
if ~isempty(model.guess)
	guess = model.guess;
else
	guess = zeros(1, model.par_num);
end

% Define parameter indices
i_a1 = 1;
i_x1 = 2;
i_y1 = 3;
i_t1 = 4;
i_s1 = 5;
i_a2 = 6;
i_x2 = 7;
i_y2 = 8;
i_t2 = 9;
i_a3 = 10;
i_y3 = 11;

%% Calculate general features
[max_val, max_ind] = max(data);
[min_val, ~] = min(data);
amp = max_val - min_val;

%% Guess parameters for first parabola

% Search for local extrema before global maximum
thresh_change = 0.08 * amp;

FIRST_IS_MAX = 1;
FIRST_IS_MIN = -1;
first_extremum = 0;

current_max = data(1);
current_max_start = current_max;

current_min = data(1);
current_min_start = current_min;

total_max = -inf;
i_total_max = 0;

total_min = inf;
i_total_min = 0;

prev_val = data(1);

for iP = 2:max_ind
	val = data(iP);

	if val > prev_val
		if ~isfinite(current_max)
			current_max = val;
			current_max_start = prev_val;
		elseif val > current_max
			current_max = val;
			if current_max_start + thresh_change <= current_max
				if first_extremum == FIRST_IS_MAX && isfinite(total_min)
					break
				elseif ~isfinite(total_min)
					first_extremum = FIRST_IS_MAX;
				end
				current_min = NaN;
				current_min_start = NaN;
				if current_max > total_max
					total_max = current_max;
					i_total_max = iP;
				end
			end
		end

	elseif val < prev_val
		if ~isfinite(current_min)
			current_min = val;
			current_min_start = prev_val;
		elseif val < current_min
			current_min = val;
			if current_min_start - thresh_change >= current_min
				if first_extremum == FIRST_IS_MIN && isfinite(total_max)
					break
				elseif ~isfinite(total_max)
					first_extremum = FIRST_IS_MIN;
				end
				current_max = NaN;
				current_max_start = NaN;
				if current_min < total_min
					total_min = current_min;
					i_total_min = iP;
				end
			end
		end
	end

	prev_val = val;
end

% Guess parameters
if first_extremum == FIRST_IS_MAX && isfinite(total_min)
	x1 = t(i_total_max);
	y1 = total_max;
	a1 = (data(1) - y1) / (t(1) - x1)^2;

elseif first_extremum == FIRST_IS_MIN && isfinite(total_max)
	x1 = t(i_total_min);
	y1 = total_min;
	a1 = (data(1) - y1) / (t(1) - x1)^2;

else
	x1 = t(max_ind);
	y1 = max_val;
	if max_ind == 1
		a1 = 0;
	else
		a1 = (data(1) - y1) / (t(1) - x1)^2;
	end
end

guess(i_a1) = a1;
guess(i_x1) = x1;
guess(i_y1) = y1;

%% Guess parameters for step
y3 = data(end);

% Search last decay over at least 3% of amplitude
i_decay_max = 0;
h_decay_current = 0;
h_decay_max = 0;
d_old = data(end);
for i_decay_current = length(data)-1:-1:max_ind
	d = data(i_decay_current);
	if d > d_old
		h_decay_current = h_decay_current + d - d_old;
	elseif d < d_old
		if h_decay_current > h_decay_max
			h_decay_max = h_decay_current;
			i_decay_max = i_decay_current + 1;
		end
		if h_decay_max >= 0.03 * amp
			break
		end
		h_decay_current = 0;
	end
	d_old = d;
end
if h_decay_max < 0.05 * amp && h_decay_current > h_decay_max
% 	h_decay_max = h_decay_current;
	i_decay_max = i_decay_current;
	h_decay_max = data(i_decay_max) - data(end);
% 	for i_next = i_decay_max+1:length(data)
% 		if data(i_next) < data(i_decay_max)
% 			h_decay_max = data(i_decay_max) - data(i_next);
% 			break
% 		end
% 	end
end
% At this point, `i_decay_max` is the start of the decay where we want to
% settle the step.

% Calculate decay steepness based on half-life time
y_half = data(i_decay_max) - h_decay_max / 2;
[~,i_half] = min(abs(data(i_decay_max:end) - y_half));
i_half = i_half + i_decay_max - 1;
t2 = t(i_decay_max);
t_half = t(i_half) - t(i_decay_max);
dy_half = (data(i_decay_max) - data(i_half)) / h_decay_max;

if dy_half <= 0
	a3 = 1;
else
	a3 = - log(dy_half) / t_half;
end

guess(i_t2) = t2;
guess(i_y3) = y3;
guess(i_a3) = a3;

%% Guess parameters for second parabola
% Check height of trace between maximum and t2
idx_t1 = 3 * max_ind;
for m = 3:-1:1
	idx_t1 = m * max_ind;
	if idx_t1 >= i_decay_max
		continue
	else
		break
	end
end
t1 = t(idx_t1);

idx_middle = idx_t1 + floor((idx_t1 + i_decay_max) / 2) - 1;
if idx_middle < 1 || idx_middle > length(data)
	idx_middle = idx_t1;
end
y_middle = data(idx_middle);

% Check if vertex of second parabola is at middle or at t2
if abs(y_middle - data(idx_t1)) < abs(y_middle - data(i_decay_max))
	x2 = t(idx_middle);
	y2 = y_middle;
	x_ref = t2;
	y_ref = data(i_decay_max);
else
	x2 = t2;
	y2 = data(i_decay_max);
	x_ref = t(idx_middle);
	y_ref = y_middle;
end

a2 = (y_ref - y2) / (x_ref - x2)^2;

guess(i_t1) = t1;
guess(i_s1) = 1;
guess(i_a2) = a2;
guess(i_x2) = x2;
guess(i_y2) = y2;


%S = struct('min', par_min, 'max', par_max', 'guess', guess);
S = struct('guess', guess);
end

