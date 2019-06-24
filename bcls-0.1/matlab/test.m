function test
% ----------------------------------------------------------------------
% aprod.m: Test driver for Matlab interface to BCLS.
% ----------------------------------------------------------------------
% Calls matrix-vector routine aprod.m
%
% Copyright 2006, Michael P. Friedlander, Univ of British Columbia
%
% $Revision: 273 $ $Date: 2006-09-04 15:59:04 -0700 (Mon, 04 Sep 2006) $

rand('state',0);            % Make experiments reproducible.
options.print_level = 2;    % Increase print level.

% ----------------------------------------------------------------------
% Explicit sparse matrix.
% ----------------------------------------------------------------------

fprintf('\nBCLS Test 1: Explicit random matrix\n');
fprintf(  '-----------------------------------\n');

m    =  1000;
n    =   500;
A    =  sprand(m, n, 0.01);
b    =  randn(m, 1);
bl   =  zeros(n, 1);
bu   =  inf*ones(n,1);
c    =  randn(n, 1);
x0   =  zeros(n, 1);
y0   =  zeros(m, 1);
z0   =  zeros(n, 1);
damp =  0;

options.major_itns  = 40;
options.usolve      = @(mode, ix, z)usolve(mode, ix, z, A);

[Xsol, g, info] = bcls(A, b, bl, bu, x0, c, damp, options);

% ----------------------------------------------------------------------
% Implicit matrix.
% ----------------------------------------------------------------------

fprintf('\nBCLS Test 2: Implicit matrix\n');
fprintf(  '----------------------------\n');

options.usolve = [];      % No preconditioning for this problem.

n    = 100;               % No. of columns (and rows) of A.
b    = aprod( 1, 1:n, ones(n,1) );
bl   = zeros(n,1);        % Solution must be nonnegative.
bu   = [];                % No upper bounds on variables. 
c    = [];                % No linear term.
x    = rand(n,1);         % Use random starting vector.

[xsol, g, info] = bcls(@aprod, b, bl, bu, x, c, 0, options);
