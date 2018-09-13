function Y = LATE_decay_simulate(t, params)
% Provides a peak function with bleaching decay
% 
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%	params(1): global bleaching decay rate
%	params(2): global bleaching decay amplitude
% 	params(3): offset of the function
% 	params(4): peak height
% 	params(5): temporal position of peak
% 	params(6): steepness of left side of peak
% 	params(7): steepness of right side of peak
% 	params(8): width of left side of peak
% 	params(9): width of right side of peak
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

	C0 = params(1);
	C0_amp = params(2);
	A = params(3);
	B = params(4);
	t_0 = params(5);
	a1 = params(6);
	a2 = params(7);
	b1 = params(8);
	b2 = params(9);

	idx_left = t < t_0;

	Y = zeros(length(t), 1);

	Y( idx_left) = exp( -((t_0 - t( idx_left)) ./ a1).^b1 );
	Y(~idx_left) = exp( -((t(~idx_left) - t_0) ./ a2).^b2 );
	Y = Y * B + A + C0_amp .* decay(t, 0, C0);
end


function D = decay(t, t0, gamma)
D = exp( (t0 - t) .* gamma );
end
