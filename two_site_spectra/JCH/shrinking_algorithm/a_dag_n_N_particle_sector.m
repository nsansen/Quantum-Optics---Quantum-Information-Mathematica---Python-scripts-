function ad = a_dag_n_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>
ad = kron(speye(2), a_n_dag(n_max));
ad(n_max+1,:) = [];
ad(:,n_max+1) = [];

end