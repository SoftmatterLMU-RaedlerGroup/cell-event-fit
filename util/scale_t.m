function [ t_scaled ] = scale_t( t, t_scale, f_ind, t_offset )
%SCALE_T properly scales time vector to fit a given unit
%
% Arguments:
% ==========
%	t			numeric time vector to be scaled
%	t_scale		scale factor (empty, scalar or vector with length equal to
%				number of data files) [optional]
%	f_ind		index of current data file [optional]
%	t_offset	offset value to be added to scaled vector [optional]
%
% Returns:
% ========
%	t_scaled	scaled time vector
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

%% Get default values
if nargin < 2
	t_scale = [];
end
if nargin < 3
	f_ind = 0;
end
if nargin < 4
	is_offset_defined = false;
	t_offset = 0;
else
	is_offset_defined = true;
end

%% Get correct value for t_scale
if length(t_scale) > 1 && f_ind > 0 && f_ind <= length(t_scale)
	t_scale = t_scale(f_ind);
elseif length(t_scale) < 1
	t_scale = get_t_scale(t);
end

%% Get correct value for t_offset
if is_offset_defined
	if length(t_offset) > 1
		try
			t_offset = t_offset(f_ind);
		catch
			error('f_ind musst be an index to t_offset.')
		end
	elseif isempty(t_offset)
		t_offset = 0;
	end
end

%% Calculate t_scaled
t_scaled = t .* t_scale + t_offset;


	function [ t_scale ] = get_t_scale( t )
	%GET_T_SCALE Heuristically finds a scale value for a time vector
	%   t_scale is a heuristically found factor for converting the time
	%	vector t to unit hours:
	%		t .* t_scale = t[h]

	t_max = max(t);

	if t_max < 50
		% Assume that t is already in unit 'hours'
		t_scale = 1;
	elseif t_max < 1000
		% Assume that t is in unit 'frames', assume frame rate 1/10min
		t_scale = 1 / 6;
		if ~is_offset_defined && t(1) == 1
			t_offset = -1;
		end
	else
		% Large values; assume that t is in unit 'seconds'
		t_scale = 1 / 3600;
	end
	end % end of `get_t_scale`
end % end of `scale_t`