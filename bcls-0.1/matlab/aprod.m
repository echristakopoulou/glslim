function [r] = aprod(mode, ix, z);
% ----------------------------------------------------------------------
% aprod.m: Matrix-vector routine for testing Matlab interface to BCLS.
% ----------------------------------------------------------------------
% Called by test.m. 
%
% Copyright 2006, Michael P. Friedlander, Univ of British Columbia
%
% $Revision: 268 $ $Date: 2006-05-24 13:26:44 -0700 (Wed, 24 May 2006) $

% Note that the length of z will change depending on the mode:
% If mode = 1, then length(z) = n, and
% if mode = 2, then length(z) = m.
% In this particular example, m = n.
n = length(z);

if mode <= 0
%  Aprod is being called in one of two modes:
%     "initialization" (mode = -1)
%  or "termination"    (mode = -2).
%  In either case, return some dummy value for r in order to avoid
%  a warning message from Matlab.
   r = [];

elseif mode == 1 % r = A(:,ix) * x(ix)

%  The vector of indices  ix  indicates which columns of A need to
%  be involved in the mat-vec product.  One easy way of dealing
%  with this is just to zero out all components of x that are not
%  in ix.
  
   x     = zeros(n,1);
   x(ix) = z(ix);
  
   r = 4 * x;
   r(2:n) = r(2:n) - 2 * x(1:n-1);
   r(1:n-1) = r(1:n-1) - x(2:n);

elseif mode == 2 % r = A'y

% BCLS will ignore all components of r that are not in ix.  So we
% are free to ignore the index vector ix in this case.

  r = 4 * z;
  r(1:n-1) = r(1:n-1) - 2 * z(2:n);
  r(2:n) = r(2:n) - z(1:n-1);
  
end
