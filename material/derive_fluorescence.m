%% derive fluorescence.m
%
% This script calculates the first, second, and third derivative of a given
% curve and plots them.
%
% Required:
% t		Time vector
% d		Data matrix
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

%% Allocate memory
dn = zeros(size(d));
d1 = zeros(size(d));
d2 = zeros(size(d));
d3 = zeros(size(d));

%% Normalize
for j = 1:size(d,2)
	d_max = max(d(:,j));
	d_min = min(d(:,j));
	amp = d_max - d_min;
	
	for i = 1:size(d,1)
		dn(i,j) = ( d(i,j) - d_min ) / amp;
	end
end

%% Derive
for j = 1:size(d,2)			% iterate over cells
	
	% first derivative
	for i = 2:size(d,1)
		d1(i,j) = dn(i,j) - dn(i-1,j);
	end
	d1(1,:) = NaN;
	
	% second derivative
	for i = 3:size(d,1)
		d2(i,j) = d1(i,j) - d1(i-1,j);
	end
	d2(1:2,:) = NaN;
	
	% third derivative
	for i = 4:size(d,1)
		d3(i,j) = d2(i,j) - d2(i-1,j);
	end
	d3(1:3,:) = NaN;
end

%% Plot
subplot(4,1,1)
plot(t,dn)
title('Simulation')
xlabel('Time [h]')
ylabel('Fluorescence [a.u.]')

subplot(4,1,2)
plot(t,d1)
title('First derivative')
xlabel('Time [h]')
ylabel('Slope [a.u.]')

subplot(4,1,3)
plot(t,d2)
title('Second derivative')
xlabel('Time [h]')
ylabel('Curvature [a.u.]')

subplot(4,1,4)
plot(t,d3)
title('Third derivative')
xlabel('Time [h]')
ylabel('Change of Curvature [a.u.]')