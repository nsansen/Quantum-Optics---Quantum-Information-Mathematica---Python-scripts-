function a = a_n_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>

a = kron(speye(2), a_n(n_max));
a(n_max+1,:) = [];
a(:,n_max+1) = [];

end