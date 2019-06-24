function r = usolve(mode, ix, z, A)

persistent U

if mode == -1     % Initialize the preconditioner.
  U = qr( A(:,ix), 0 );
  r = [];
elseif mode == -2 % Preconditioner no longer needed.
  clear U;
elseif mode == 1  % Solve  U r = z.
  r = U \ z;
elseif mode == 2  % Solve  U'r = z.
  r = U'\ z;
end
