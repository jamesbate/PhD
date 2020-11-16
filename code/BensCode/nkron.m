function ans = nkron(first, varargin)
% nkron.m : Kronecker tensor product for 2 or more arguments
% originally written by Julio Barreiro

ans = kron(first,varargin{1});
if nargin > 2
  ans = nkron(ans,varargin{2:end});
end
