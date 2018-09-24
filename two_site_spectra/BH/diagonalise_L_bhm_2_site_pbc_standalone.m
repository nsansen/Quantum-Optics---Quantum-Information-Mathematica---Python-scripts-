function diagonalise_L_bhm_2_site_pbc_standalone(job_number)
% Script to build the flattened Liovillian and the diagonalise it:
% clear

if ischar(job_number), job_number = str2num(job_number); end;

save_prefix = 'BH_ED_drive_1_part_res_increase_J_compare_iTEBD_';
save_name = [save_prefix num2str(job_number) '.mat'];

% Parameters:

n_max = 3  % Photon truncation parameter
M = 2  % Number of sites in ringa

dim_H = (n_max+1)^M;

gamma = 0.2;

Omega = 0.03;

U = 1;
J_list = logspace(-2, 1, 51);

% J = J_list(job_number);

var_list = J_list;

rho_store = zeros(dim_H, dim_H, length(var_list));

for loop = 1:length(var_list)
    
    disp(num2str(length(var_list) - loop))
    
    % Set what needs to be set:
    J = J_list(loop);
    delta_omega_d = 2*J;
    
    rho_ss = get_steady_state_rho(M, n_max, delta_omega_d, J, U, Omega, gamma);
    rho_store(:,:,loop) = rho_ss;
    
end

% save(save_name, 'n_max', 'M', 'delta_omega_d_list', 'gamma', 'Omega', 'J', 'U', 'rho_store')
save(save_name, 'n_max', 'M', 'gamma', 'Omega', 'J_list', 'U', 'rho_store')
disp(['Saved data to file ' save_name])

end

%%%%%%%%%%%%%%%%5 Auxillary functions %%%%%%%%%%%%%%%%%%%%%%%
function rho_ss = get_steady_state_rho(M, n_max, delta_omega_d, J, U, Omega, gamma)
% Function to return the steady state given the system parameters

% Choose the boundary conditions:
H = bhm_hamiltonian_multi_site_pbc(M, n_max, delta_omega_d, J, U, Omega);
% H = bhm_hamiltonian_multi_site_obc(M, n_max, delta_omega_d, J, U, Omega);

L = flattened_liouvillian(H, n_max, gamma, M);

[min_evector, min_evalue] = eigs(L,1,'SM');
dim = sqrt(length(min_evector));
rho_ss = reshape(min_evector, dim, dim);
rho_ss = rho_ss/trace(rho_ss);

end
%%%%%%%%%%%%%%%
function L = flattened_liouvillian(H_unitary, n_max, gamma, M)
% Function to return the `flattened' Liouvillian operator, size N^2 times
% N^2, to premultiply a flattened density vector (size N^2 times 1)
% H_unitary is the Hermitian part Hamiltonian
% n_max is photon truncation parameter
% gamma is the loss rate from each cavity
% M is the system size

% Construct the creation and annhilation operators for the field:
c = a_n(n_max+1);
c_dag = a_n_dag(n_max+1);

L = 1/i*(spre(H_unitary) - spost(H_unitary));    % The Hermitian commutator

for m = 1:M     % Loop over the lattice sites
    
    c_m = tensor_matrix(c, M, m);
    c_dag_m = tensor_matrix(c_dag, M, m);
    
    L = L + gamma/2*(2*spre(c_m)*spost(c_dag_m) - spre(c_dag_m*c_m) - spost(c_dag_m*c_m));    % The Liouvillian for the field
    
end

end
%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%5
function H = bhm_hamiltonian_multi_site_obc(M, n_max, delta_omega_d, J, U, Omega)
% Function to return the BH Hamiltonian for OPEN boundary conditions
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

end
%%%%%%%%%%%%%%%%
function operator = a_n(n_max)

% Function to return the destruction operator in a Fock state basis with
% number states from 0 to n-1

operator = spdiags(sqrt([0:n_max].'), 1, n_max, n_max);

end
%%%%%%%%%%%%%%%%
function operator = a_n_dag(n_max)

% Function to return the creation operator in a Fock state basis with
% number states from 0 to n-1

operator = spdiags(sqrt([1:n_max-1].'), -1, n_max, n_max);

end
%%%%%%%%%%%%%%
function matrix_big = tensor_matrix(matrix, M, index)
% Function to return the tensor product: I * I * ... * matrix * I * ... *
% I, where matrix is in the index'th position, and '*' denotes a tensor
% multiplication, and there are M matrices all together in the product

id = speye(size(matrix)); % The identity matrix to be tensored successively
matrix_big = 1; % Initialise the output matrix (NB kron(1, matrix) = matrix

for k = 1:M
    
    if (k == index)
        next_matrix = matrix;
    else
        next_matrix = id;
    end
    
    matrix_big = kron(matrix_big, next_matrix);
    
end

end
%%%%%%%%%%%%%%%%%%%%%%
function S = spost(A)
% Function to find the superoperator associated with premultiplying density
% matrix by A

% Get dimension of operator:
d = size(A,2);

S = kron(A.', speye(d));

end
%%%%%%%%%%%%%%%%%%%%%%
function S = spre(A)
% Function to find the superoperator associated with premultiplying density
% matrix by A

% Get dimension of operator:
d = size(A,2);

S = kron(speye(d), A);

end