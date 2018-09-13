function Y = LATE_simulate(t, params)
% Provides a sigmoid function with parameters
% 
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
% 	params(1): offset of the function
% 	params(2): peak height
% 	params(3): temporal position of onset
% 	params(4): steepness of onset
% 	params(5): steepness of decay
% 	params(6): constant end value
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

	A = params(1);
	B = params(2);
	t_0 = params(3);
	alpha = params(4);
	gamma = params(5);
	b = params(6);

	Y = zeros(length(t), 1);

	for i = 1:length(t)
		Y(i) = A + B * sigmoid(t(i)) ...
			* ( b + decay(t(i)) );
	end

	function S = sigmoid(t)
		S = 1 / ( 1 + exp( (t_0 - t) * alpha ) );
	end

	function D = decay(t)
		D = exp( (t_0 - t) * gamma );
		
		if D > 100
			D = Inf;
		end
	end
end
