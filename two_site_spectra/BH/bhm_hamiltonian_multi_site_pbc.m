function H = bhm_hamiltonian_multi_site_pbc(M, n_max, delta_omega_d, J, U, Omega)
% Function to return the BH Hamiltonian for periodic boundary conditions
% M - number of lattice sites
% n_max - boson truncation
% delta_omega_d - difference between boson frequency and driving frequency:
% delta_omega_d = omega_boson - omega_laser
% delta_omega_0 - difference between atomic frequency and driving frequency:
% n_max (maximum photon number at which we truncate the
% Fock space basis at each site), A (hopping energy)

b = a_n(n_max + 1); % Destruction operator for one site
b_dag = a_n_dag(n_max + 1); % Creation operator

num = b_dag*b;    % Number operator

ket_g = [0;1];
ket_e = [1;0];

H = sparse((n_max+1)^M,(n_max+1)^M);   % Initialise Hamiltonian

for m = 1:M-1     % Loop over the lattice sites
    
    b_m = tensor_matrix(b, M, m);
    b_dag_m = tensor_matrix(b_dag, M, m);
    b_m_plus_1 = tensor_matrix(b, M, m+1);
    b_dag_m_plus_1 = tensor_matrix(b_dag, M, m+1);
    num_m = tensor_matrix(num, M, m);
    
    hopping_term = b_dag_m*b_m_plus_1 + b_dag_m_plus_1*b_m;
    
    H = H + delta_omega_d*num_m - J*hopping_term + U*b_dag_m*b_dag_m*b_m*b_m + Omega*(b_m + b_dag_m);   % Assemble the Hamiltonian for each lattice site contribution
    
end

% Need to also include the last on-site term:
b_M = tensor_matrix(b, M, M);
b_dag_M = tensor_matrix(b_dag, M, M);
num_M = tensor_matrix(num, M, M);

H = H + delta_omega_d*num_M + U*b_dag_M*b_dag_M*b_M*b_M + Omega*(b_M + b_dag_M);   % Add the last site's on-site interaction contribution

% Now need the terms which make periodic boundary conditions (i.e. hopping
% between sites 1 and M)
b_1 = tensor_matrix(b, M, 1);
b_dag_1 = tensor_matrix(b_dag, M, 1);
hopping_term_pbc_pbc = b_dag_M*b_1 + b_dag_1*b_M;

H = H - J*hopping_term_pbc_pbc;

end