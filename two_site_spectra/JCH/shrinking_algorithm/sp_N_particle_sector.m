function sp = sp_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>

sp = kron([0 1;0 0], speye(n_max+1));
sp(n_max+1,:) = [];
sp(:,n_max+1) = [];

end