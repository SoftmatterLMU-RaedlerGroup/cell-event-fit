function Y = doubleparabola_TMRM_simulate(t, params)
%DOUBLEPARABOLA_TMRM_SIMULATE Summary of this function goes here
%   Detailed explanation goes here

%% Get parameters
% First parabola
a1 = params(1);
x1 = params(2);
y1 = params(3);

% Transition to second parabola (time t1, smoothness s1)
t1 = params(4);
s1 = params(5);

% Second parabola
a2 = params(6);
x2 = params(7);
y2 = params(8);

% Transition to exponential decay
t2 = params(9);

% Exponential decay (steepness a3, final level y3)
a3 = params(10);
y3 = params(11);

%% Calculation
Y = zeros(size(t));

idx = t < t2;
Y(idx) = combine(t(idx), a1, x1, y1, a2, x2, y2, t1, s1);

idx = t >= t2;
Y(idx) = y3 + (combine(t2, a1, x1, y1, a2, x2, y2, t1, s1) - y3) .* exp(-a3 * (t(idx) - t2));

end


function y = combine(t, a1, x1, y1, a2, x2, y2, t1, s1)
	%sigmoid = 1 ./ (1 + exp(- s1 * (t - t1)));
% 	s1 = abs(s1);
% 	sigmoid_ = [zeros(size(t(t <= t1 - s1 / 2))); ...
% 		1/s1 * (t((t > t1 - s1 / 2) & (t < t1 + s1 / 2)) - (t1 - s1 / 2)); ...
% 		ones(size(t(t >= t1 + s1 / 2)))];
	sigmoid_ = sigmoid(t, t1, s1);

	y = (1 - sigmoid_) .* parabola(t, a1, x1, y1);
	y = y + sigmoid_ .* parabola(t, a2, x2, y2);
end

function y = sigmoid(t, t0, s)
	s = abs(s);
	t_start = t0 - s / 2;
	t_end = t0 + s / 2;

	y = zeros(size(t));

	idx = (t > t_start) & (t < t_end);
	y(idx) = (sin((t(idx) - t0) ./ s * pi) + 1) * .5;

	idx = t >= t_end;
	y(idx) = 1;
end


function y = parabola(t, a, x0, y0)
	y = a * (t - x0).^2 + y0;
end