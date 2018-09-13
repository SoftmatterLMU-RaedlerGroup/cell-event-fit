function Y = pi_heaviside_simulate(t, params)
% Provides a heaviside function with parameters
% 
% Arguments:
% t: time
% params:  coefficients of the heaviside functions:
%    params(1): offset of the function
%    params(2): scaling factor of the function
%    params(3): jump point of the fuction
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

Y = zeros(length(t), 1);

% Y is theta-shaped: Y = params(1) + params(2) * theta(x-params(3))
for i = 1:length(t)
    if t(i) < params(3)
        Y(i) = params(1);
    elseif t(i) > params(3)
        Y(i) = params(2) + params(1);
    elseif t(i) == params(3)
        Y(i) = params(2) / 2 + params(1);
    end
end

