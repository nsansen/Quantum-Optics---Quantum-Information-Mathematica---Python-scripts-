function S = spre(A)
% Function to find the superoperator associated with premultiplying density
% matrix by A

% Get dimension of operator:
d = size(A,2);

S = kron(speye(d), A);

end