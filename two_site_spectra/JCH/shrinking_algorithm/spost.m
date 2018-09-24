function S = spost(A)
% Function to find the superoperator associated with premultiplying density
% matrix by A

% Get dimension of operator:
d = size(A,2);

S = kron(A.', speye(d));

end