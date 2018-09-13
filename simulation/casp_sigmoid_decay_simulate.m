function Y = casp_sigmoid_decay_simulate(t, params)
% Provides a sigmoid function with parameters
% 
% Arguments:
% t: time
% params:  coefficients of the sigmoid functions:
%    params(1): offset of the function
%    params(2): height difference
%    params(3): ascent position
%    params(4): peak height
%    params(5): decay constant
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


for i = 1:length(t)
    Y(i) = params(1) ...
        + params(2) * ( 1 - 1 / (1 + (t(i) / params(3))^20 ) ) ...
        * params(4) * exp( (params(3)-t(i)) * params(5) );
end
