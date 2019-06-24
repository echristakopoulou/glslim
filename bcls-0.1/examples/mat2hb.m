function mat2hb( file, title, key, A, b, bl, bu, c, x )

% ----------------------------------------------------------------------
% mat2hb.m: Convert data for a BCLS problem to Harwell-Boeing format.
% ----------------------------------------------------------------------
% mat2hb( file, title, key, A, b, bl, bu, c, x )
%
% Write out the data for the problem
%
%   minimize   || ( A    )     (b) ||
%              || (      ) x - ( ) ||  +  c'x  s.t.  bl <= x <= bu,
%      x       || (damp I)     (0) ||
%
% to a format usable by the BCLS commandline utility bcsol.
%
% file        (required) is the root of the two filenames written:
%
%             file.hbf:  The Harwell-Boeing file that describes
%                        A  and  b.
%             file.opt:  The file that describes
%                        bl, bu, c, x.
%                        These 4 n-vectors are simply listed
%                        consecutively in the file.  No missing
%                        entries are allowed.
%
% title       (optional) is the "title" field of the Harwell-Boeing file.
%
% key         (optional) is the "key" field of the Harwell-Boeing file.
%
% A           (required) is the sparse m-by-n matrix.
%
% b           (required) is the right-hand side.
%
% bl          (optional) is the vector of lower bounds.  May have
%             elements that are -inf.  If bl is missing or empty, then
%               bl(1:n) = -inf  is used.
%
% bu          (optional) is the vector of upper bounds.  May have
%             elements that are +inf.  If bu is missing or empty, then
%               bu(1:n) = +inf  is used.
%
% c           (optional) is the linear term.  If  c  is missing or
%             empty, then
%               c(1:n) = 0  is used.
%
% x           (optional) is an estimate of the solution.  If x is
%             missing or empty, then
%               x(1:n) = 0  is used.
%
% $Revision: 290 $ $Date: 2007-03-04 22:01:54 -0800 (Sun, 04 Mar 2007) $
% ----------------------------------------------------------------------

if nargin < 5
   error('The first 5 arguments are required');
end

[m, n] = size(A);

% Check strings.
if ~ischar(file) || length(file)==0
   error('"file" must be a nonempy string')
end
if ~ischar(title) || length(title)==0
   error('"title" must be a nonempty string')
end
if ~ischar(key) || length(key)==0
   error('"key" must be a nonempty string')
end
   
% Check that required data isn't empty.
if isempty(A)
   error('A  cannot be empty.');
elseif ~issparse(A)
   error('A  must be sparse.');
end

% Check b.
if isempty(b)
   error('b  cannot be empty');
elseif ~isequal( size(b), [m 1] )
   error('b must be m-by-1.');
end

% Check bl.
if ~exist('bl','var') || isempty(bl)
   bl = repmat(-1e20, n, 1 );
elseif ~isequal( size(bl), [n 1] )
   error('bl must be n-by-1.');
else
   bl = max( bl, -1e20 );
end

% Check bu.
if ~exist('bu','var') || isempty(bu)
   bu = repmat( 1e20, n, 1 );
elseif ~isequal( size(bu), [n 1] )
   error('bu must be n-by-1');
else
   bu = min( bu, 1e20 );
end

% Check c.
if ~exist('c','var') || isempty(c)
   c = zeros(n,1);
elseif ~isequal( size(c), [n 1] )
   error('c must be n-by-1.')
end

% Check x.
if ~exist('x','var') || isempty(x)
   x = zeros(n,1);
elseif ~isequal( size(x), [n 1] )
   error('x must be n-by-1.');
end

% Don't allow sparse vectors.
if issparse(b ), b  = full(b ); end
if issparse(bl), bl = full(bl); end
if issparse(bu), bu = full(bu); end
if issparse(c ), c  = full(c ); end
if issparse(x ), x  = full(x ); end

% Call mat2hbmex.
mat2hbmex(file, title, key, A, b, bl, bu, c, x );