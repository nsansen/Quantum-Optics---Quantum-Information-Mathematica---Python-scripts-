function delta_omega_c = return_driving_laser_detuning_lowest_two_particle_resonance(g, delta, J)

M = 2;  % Start with just two sites
n_max = 2;  % We'll only be interested in a maximum of 2 particle in the whole system, so this should also be 2

omega_d = 10000;
omega_0 = omega_d - delta;

H = jch_hamiltonian_multi_site_pbc(M, n_max, omega_d, omega_0, g, J, 0);

[v,d] = eig(full(H));

evs = diag(d);

two_particle_state_indices = find((evs > 0.95*omega_d*M) & (evs < 1.05*omega_d*M));

state_frequency = evs(two_particle_state_indices(1));   % This is the frequency of the state in the bare frame

delta_omega_c = omega_d - state_frequency/M;    % This is the driving laser detuning necessary to drive this state

end