%  varargout = notchFilterDesign(fn,varargin)
%
%  DESCRIPTION
%  Designs a bank of notch filters at the normalised frequencies FN with a
%  normalilsed bandwidth BWN and filter order FILTORDER. Depending on the
%  numer of output arguments, NOTCHFILTERDESIGN will return: 1. A 'cascade'
%  class filter H consisting in a sequence of order FILTORDER d2fsos filters
%  in cascade; 2. The B and A coefficients of a sequence of order-2 filters
%  put in cascade by convolving their time impulse responses. For this
%  second method FILTORDER is ignored. Both methods give the same result
%  for FILTORDER = 2.
% 
%  INPUT VARIABLES
%  - fn: vector of normalised frequencies at which the notch will be applied.
%  - bwn (varargin{1}): vector of normalised bandwidths for the specified
%    frequencies FN. BWN can be a one-element vector, in which case the
%    same value will be applied to all frequencies.
%  - filtOrder (varargin{2}): order for the notch filters.
%   
%  OUTPUT VARIABLES
%  - h (varargout{1}): 'cascade' class filter containing a sequence of
%    df2sos notch filters of order FILTORDER.
%  - [b,a] (varargout{1},varargout{2}): B and A coefficients of a sequence
%    of order-2 notch filters put in cascade.
%
%  FUNCTION CALL
%   ...  = NOTCHFILTERDESIGN(fn)
%   ...  = NOTCHFILTERDESIGN(fn,bwn)
%   ...  = NOTCHFILTERDESIGN(fn,bwn,filtOrder)
%   h    = NOTCHFILTERDESIGN(fn,...)
%  [b,a] = NOTCHFILTERDESIGN(fn,...)
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)

%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  23 May 2021

function varargout = notchFilterDesign(fn,varargin)
    
% Error Control
narginchk(0,3)
if ~isempty(fn) && (~isnumeric(fn) || ~isvector(fn) ...
        || any(fn < 0 | fn > 1))
    error(['The normalised input frequency FN must be numeric vector '...
        'with values between 0 and 1'])
end

% Retrieve Input Arguments
switch nargin
    case 1
        bwn = 0.02*max(fn)*ones(length(fn),1); % quality factor Q = MAX(FN)/BWN = 50
        filtOrder = 2;
    case 2
        bwn = varargin{1};
        filtOrder = 2;
    case 3
        bwn = varargin{1};
        filtOrder = varargin{2};
end

% Error Control
if ~isnumeric(bwn) || ~isvector(bwn) || any(bwn < 0)
    bwn = 0.02*max(fn); % quality factor Q = MAX(FN)/BWN = 50
    warning(['The normalised bandwidth BWN must be a numeric value '...
        'or vector of values higher than 0. BWN = %0.1f will be used'],bwn)
end
if length(bwn) == 1
    bwn = bwn*ones(length(fn),1);
    if length(bwn) ~= length(bwn)
        error('LENGTH(BWN) must be 1 or equal to LENGTH(FN)')
    end
end
if ~isnumeric(filtOrder) || rem(filtOrder,1) || length(filtOrder) > 1 ...
        || filtOrder < 2
    filtOrder = 2;
    warning('The filter order must be an integer higher than 2')
end

switch nargout
    case 1
        % Create Cascade of Notch Filters
        nFreq = length(fn);
        h0 = dfilt.df2sos(1,1); % initialise filter
        for k = 1:nFreq
            h0(k) = design(fdesign.notch('N,F0,BW',filtOrder,fn(k),bwn(k)));
            if ~isstable(h0(k))
                error(['Notch filter at FN = %0.3e (and below) '...
                    'is not stable'],fn(k))
            end
        end
        varargout{1} = dfilt.cascade(h0); % df2sos notch filters in cascade
    case 2
        % Create Cascade of Notch Filters (Order 2 only)
        nFreq = length(fn);
        a0_prev = 1;
        b0_prev = 1;
        for k = 1:nFreq
            [b0,a0] = iirnotch(fn(k),bwn(k));
            b0_prev = conv(b0_prev,b0);
            a0_prev = conv(a0_prev,a0);
        end
        varargout{1} = b0_prev; % b coefficients of notch filters in cascade
        varargout{2} = a0_prev; % a coefficients of notch filters in cascade
end


