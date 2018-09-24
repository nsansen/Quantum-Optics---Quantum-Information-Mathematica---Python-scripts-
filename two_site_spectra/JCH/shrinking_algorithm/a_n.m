function operator = a_n(n_max)

% Function to return the destruction operator in a Fock state basis with
% number states from 0 to n-1

operator = spdiags(sqrt([0:n_max+1].'), 1, n_max+1, n_max+1);

end