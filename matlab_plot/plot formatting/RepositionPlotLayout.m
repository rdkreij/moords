function [ pos ] = RepositionPlotLayout( varargin )

if nargin == 3
    P = varargin{1};
    rows = varargin{2};
    cols = varargin{3}; 
elseif nargin == 1
    P = varargin{1};
    rows = 1:size(P, 1);
    cols = 1:size(P, 2);
else
    error('Unrecognised inputs')
end

dc=diff(cols)==1;
if ~all(dc) && ~isempty(dc)
    error
end
dr=diff(rows)==1;
if ~all(dc) && ~isempty(dc)
	error
end

LHX=P{rows(end),cols(1)}(1);
LHY=P{rows(end),cols(1)}(2);

RHX=P{rows(1),cols(end)}(1) + P{rows(1),cols(end)}(3);
RHY=P{rows(1),cols(end)}(2)   + P{rows(1),cols(end)}(4);

pos = [LHX LHY RHX-LHX RHY-LHY];

end

