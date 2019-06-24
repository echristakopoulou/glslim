function [x, g, info] = bcls( A, b, bl, bu, x, c, damp, options );

% ----------------------------------------------------------------------
% bcls.m: Bound constrained linear least squares.
% ----------------------------------------------------------------------
% [x, g, info] = bcls( A, b, bl, bu, x, c, damp, options )
%
% or
%
% [x, g, info] = bcls( Aprod, b, bl, bu, x, c, damp, options )
%
% Solves the bound-constrained linear least-squares problem
%
%   minimize   || ( A    )     (b) ||
%              || (      ) x - ( ) ||  +  c'x  s.t.  bl <= x <= bu,
%      x       || (damp I)     (0) ||
% 
% where
%
% A        is an  m x n  dense or sparse matrix.
%
% Aprod    is a function with the prototype
%
%          [r] = Aprod( mode, ix, z )
%
%          where  r = y + A(:,ix) * z(ix)  if  mode == 1
%          and    r = x + A' * z           if  mode == 2.
%
%          Here, ix is a vector that holds a subset of indices of
%          1:n.  In the case where mode == 2, only the components of
%          r(ix) will be used; the other components will be ignored.
%
% b        is the RHS vector (length m).
%
% bl, bu   (optional) are n-vectors of lower and upper bounds on
%          x .  If these are empty or missing, then they will be
%          substitued with -inf and +inf, respectively.
%
% x        (optional) is an n-vector that gives an estimate of the
%          solution.  If x is empty or missing, then x = 0 is the
%          default.
%
% c        (optional) is an n-vector that defines the linear term.
%          If c is empty or missing, then c = 0 is the default.
%
% damp     (optional) is a scalar regularization term.  If damp is
%          empty or missing, then damp = 0 is used.
%
% options  (optional) is a structure with optional parameters that
%           control the bcls algorithm.  The special call
%
%          [options] = bcls();
%
%          returns an empty structure with valid option fields:
%
%          .n                number of columns in A
%          .print_level      0 = (none) ... = 6 (most)
%          .proj_search      0 = first exact min along projected path
%                                (default)
%                            1 = projected backtrack
%          .newton_step      0 = LSQR (default), 1 = CGLS
%          .opt_tol          default is 1e-6
%          .major_itns       default is  5 * n
%          .minor_itns       default is 10 * n
%          .Usolve           function handle to preconditioning routine.
%
%  The (optional) preconditioning routine Usolve must have the prototype
%
%          z = Usolve( mode, ix, q )
%
%  where   z  solves  U z = q  if  mode == 1
%  and     z  solves  U'z = q  if  mode == 2.
%
% As is the case with Aprod, ix is a vector that holds a subset of
% "participating" indices in 1:n.
%
%
% OUTPUTS
%
% x        is a solution (if all went well!).
%
% g        is the gradient of the objective at x.
%
% info     is a structure that describes the status of the BCLS solve,
%          and has fields
%
%          .dInf             the largest gradient among the free vars
%          .jInf             the variable index with the largest gradient
%          .major_itns       number of major iterations
%          .minor_itns       total number of minor (CG) iterations
%          .exit             exit status
%                            0 = converged
%                            1 = too many major iterations
%                            2 = too many minor iterations
%                            3 = undefined exit (bug?)
%                            4 = exit with direction of inifinite descent
%                            5 = bounds are inconsistent
%                            6 = linesearch failure
%                            These are exits requested by the user's routines:
%                            100 = Aprod
%                            110 = Usolve
%          .message          string with BCLS's exit status
%          .time             total solution time (seconds).
%
% Copyright 2006, Michael P. Friedlander, Univ of British Columbia
% 
% $Revision: 290 $ $Date: 2007-03-04 22:01:54 -0800 (Sun, 04 Mar 2007) $

% Create a structure to hold valid options.
defaultOpts.print_level = [];
defaultOpts.proj_search = [];
defaultOpts.opt_tol     = [];
defaultOpts.major_itns  = [];
defaultOpts.minor_itns  = [];
defaultOpts.usolve      = [];
defaultOpts.minor_file  = '';
defaultOpts.newton_step = [];

if nargin == 0
%  Return an empty options structure.
   x = defaultOpts;
   return
elseif nargin < 2
   error('The first 2 arguments are required');
end

% Problem size.  If A is an operator, we'll need to discover n.
m = length(b);

if     isnumeric(A),                n = size(A, 2);
elseif nargin > 2 && ~isempty(bl),  n = length(bl);
elseif nargin > 3 && ~isempty(bu),  n = length(bu);
elseif nargin > 4 && ~isempty(x ),  n = length(x );
elseif nargin > 5 && ~isempty(c ),  n = length(c );
elseif nargin > 7 && ~isempty(options) && isfield(options,'n')
                                    n = options.n;
else
   error(sprintf('No. of vars must be specified in "options" structure\n'));
end

% Make sure that A is either a function or sparse matrix.
if isa(A, 'function_handle') || isa(A, 'char')
   % Relax.  This is an implicit matrix.
elseif ~isnumeric(A)
   error('Explicit matrix  A  must be sparse.')
end
   
% Set default starting point and linear term if none given.
if nargin <  3 || isempty(bl)  , bl   = repmat(-inf, n, 1); end
if nargin <  4 || isempty(bu)  , bu   = repmat( inf, n, 1); end
if nargin <  5 || isempty(x )  , x    = zeros(n,1);         end
if nargin <  6,                  c    = [];                 end
if nargin <  7 || isempty(damp), damp = 0;                  end

% Don't allow sparse vectors.
if issparse(b ), b  = full(b ); end
if issparse(c ), c  = full(c ); end
if issparse(x ), x  = full(x ); end
if issparse(bl), bl = full(bl); end
if issparse(bu), bu = full(bu); end

% Replace all +inf and -inf with "large" numbers.
bl = max( bl, -1e20 );
bu = min( bu,  1e20 );

% Grab all user-defined options that are valid.
if nargin >= 7 && exist('options','var') && ~isempty(options)
   validFields = fieldnames( defaultOpts );
   usersFields = fieldnames( options     );
   % Check if all fields in the user-provided struc are valid.
   for i = 1 : length(usersFields)
       iField = usersFields{i};
       if ~ismember( iField, validFields )
          error('%s.%s is not a valid options field.', ...
                inputname(8),iField);
       end
   end
   for i = 1 : length(validFields)
       iField = validFields{i};
       if isfield( options, iField )
          defaultOpts.(iField) = options.(iField);
       end
   end
end

% Check that the user-defined options are valid.
if ~isa(defaultOpts.print_level,'numeric')
   error('Incorrect, print_level option');
end
if ~isa(defaultOpts.proj_search,'numeric')
   error('Incorrect proj_search option');
end
if ~isa(defaultOpts.opt_tol,'numeric')
   error('Incorrect proj_search option');
end
if ~isa(defaultOpts.major_itns,'numeric')
   error('Incorrect major_itns option');
end
if ~isa(defaultOpts.minor_itns,'numeric')
   error('Incorrect minor_itns option');
end
if ~isempty(defaultOpts.usolve) && ~isa(defaultOpts.usolve,'function_handle')
   error('Incorrect usolve option');
end
if ~isa(defaultOpts.minor_file,'char')
   error('Incorrect minor_file option');
end

% Call the MEX interface to BCLS.
[x, g, dInf, jInf, major_itns, minor_itns, ...
 status_num, status_msg, tot_time] = ...
    bclsmex( m, n, A, b, bl, bu, x, c, damp, defaultOpts );

% Assemble output data.
info.dInf       = dInf;
info.jInf       = jInf;
info.major_itns = major_itns;
info.minor_itns = minor_itns;
info.exit       = status_num;
info.message    = status_msg;
info.time       = tot_time;