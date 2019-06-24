function savehb( out_file, A, b, bl, bu )
% savehb  Create a Harwell-Boeing file from A and b.
% [] = savehb( out_file, A, b )  saves the sparse matrix A and the
% right-hand-side vector b to "out_file.ab" using the Harwell-Boeing
% format.
%
% [] = savehb( out_file, A, b, bl, bu )  also saves the bounds bl
% and bu into a two-column plain ascii file named outfile.lu.
%
% If A is dense, it's first converted to sparse format.

if nargin < 3
   error('Not enough arguments required.');
elseif nargin < 5
   error('Must include lower *and* upper bounds.');
end

if nargin == 5
   bounds =  true;
end

if ~issparse(A)
   A = sparse(A);
end

% Write out A and b.
dm2hb( [out_file '.ab'], ...
       A, b, 'Created by dm2hb', 'Matlab', [], [], 3 );

% Write out the bounds.
if bounds
   blbu  = [bl bu];
   fname = [out_file '.lu'];
   save( fname, 'blbu', '-ascii', '-double')
end