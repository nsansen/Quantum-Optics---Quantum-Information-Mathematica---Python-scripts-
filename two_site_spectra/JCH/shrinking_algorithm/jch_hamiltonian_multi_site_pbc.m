function H = jch_hamiltonian_multi_site_pbc(M, n_max, delta_omega_d, delta_omega_0, g, A, Omega)
% Function to return the BH Hamiltonian for OBC
% M - number of lattice sites
% n_max - boson truncation
% delta_omega_d - difference between boson frequency and driving frequency:
% delta_omega_d = omega_boson - omega_laser
% delta_omega_0 - difference between atomic frequency and driving frequency:
% n_max (maximum photon number at which we truncate the
% Fock space basis at each site), A (hopping energy)

b = a_n(n_max); % Destruction operator for one site
b_dag = a_n_dag(n_max); % Creation operator

% Set up identity operators for the atomic and photonic parts

eye_a = eye(2);
eye_p = eye(n_max+1);

b = kron(eye_a, b);
b_dag = kron(eye_a, b_dag);
num = b_dag*b;    % Number operator

ket_g = [0;1];
ket_e = [1;0];

% Atomic raising, lowering, and number operators
sigma_p = kron([0 1;0 0], eye_p);
sigma_m = kron([0 0;1 0], eye_p);
num_a = sigma_p*sigma_m;

dim_H = (2*(n_max+1))^M;
H = sparse(dim_H,dim_H);   % Initialise Hamiltonian

for m = 1:M-1     % Loop over the lattice sites
    
    b_m = tensor_matrix(b, M, m);
    b_dag_m = tensor_matrix(b_dag, M, m);
    b_m_plus_1 = tensor_matrix(b, M, m+1);
    b_dag_m_plus_1 = tensor_matrix(b_dag, M, m+1);
    num_m = tensor_matrix(num, M, m);
    num_a_m = tensor_matrix(num_a, M, m);
    sigma_p_m = tensor_matrix(sigma_p, M, m);
    sigma_m_m = tensor_matrix(sigma_m, M, m);
    
    hopping_term = b_dag_m*b_m_plus_1 + b_dag_m_plus_1*b_m;
    
    H = H + delta_omega_d*num_m + delta_omega_0*num_a_m - A*hopping_term + g*(sigma_p_m*b_m + sigma_m_m*b_dag_m) + Omega*(b_m + b_dag_m);   % Assemble the Hamiltonian for each lattice site contribution
    
end

% Need to also include the last on-site term:
b_M = tensor_matrix(b, M, M);
b_dag_M = tensor_matrix(b_dag, M, M);
num_M = tensor_matrix(num, M, M);
num_a_M = tensor_matrix(num_a, M, M);
sigma_p_M = tensor_matrix(sigma_p, M, M);
sigma_m_M = tensor_matrix(sigma_m, M, M);

H = H + delta_omega_d*num_M + delta_omega_0*num_a_M + g*(sigma_p_M*b_M + sigma_m_M*b_dag_M) + Omega*(b_M + b_dag_M);   % Add the last site's on-site interaction contribution

% Hopping term:
b_1 = tensor_matrix(b, M, 1);
b_dag_1 = tensor_matrix(b_dag, M, 1);
hopping_term = b_1*b_dag_M + b_dag_1*b_M;

H = H - A*hopping_term;

end