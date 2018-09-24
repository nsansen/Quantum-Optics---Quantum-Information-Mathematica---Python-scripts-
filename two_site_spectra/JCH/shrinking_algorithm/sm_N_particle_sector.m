function sm = sm_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>

sm = kron([0 0;1 0], speye(n_max+1));
sm(n_max+1,:) = [];
sm(:,n_max+1) = [];

end