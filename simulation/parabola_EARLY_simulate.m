function Y = parabola_EARLY_simulate(t, params, component)
%DOUBLEPARABOLA_EARLY_SIMULATE model function for interactive early marker model
%
% The model function consists of a parabolic part. At time `t_break`, the
% function breaks down with an exponential decay to level `y_end`.
% For times outside of the interval [`t_start`, `t_end`], a constant value
% corresponding to the nearest interval limit is returned.
%
% Parameter meaning:
%	params(1)	a_par		Quadratic coefficient of parabola
%	params(2)	x0			Time coordinate of parabola vertex
%	params(3)	y0			Vertical coordinate of parabola vertex
%	params(4)	t_break		Breakdown time (in hours)
%	params(5)	a_break		Decay coefficient
%	params(6)	y_end		Final level of decay
%	params(7)	t_start		Start of relevant data interval
%	params(8)	t_end		End of relevant data interval
%
% Arguments:
%	t			Time vector (in hours)
%	params		Parameter vector
%	component	(optional) string indicating single component to plot
%
% Returns:
%	Vector of calculated values
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

%% Get parameters
% First parabola
a_par = params(1);
x0 = params(2);
y0 = params(3);

% Transition to exponential decay
t_break = params(4);

% Exponential decay (steepness a3, final level y3)
a_break = params(5);
y_end = params(6);

% Specify valid time range
t_start = params(7);
t_end = params(8);

% If requested, compute single component
if nargin >= 3
	switch component
		case 'parabola'
			Y = parabola(t);
		case 'decay'
			Y = decay(t);
		otherwise
			error('Unknown component: %s', component);
	end
	return
end

%% Calculation
Y = zeros(size(t));

idx = t <= t_start;
Y(idx) = parabola(t_start);

idx = (t > t_start) & (t <= t_break);
Y(idx) = parabola(t(idx));

idx = (t > t_break) & (t < t_end);
Y(idx) = decay(t(idx));

idx = t >= t_end;
Y(idx) = decay(t_end);


	function y = decay(t)
		y = y_end + (parabola(t_break) - y_end) .* exp(-a_break * (t - t_break));
	end


	function y = parabola(t)
		y = a_par * (t - x0).^2 + y0;
	end
end