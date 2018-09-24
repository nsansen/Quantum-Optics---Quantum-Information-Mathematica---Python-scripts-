function two_site_shrinking_algorithm_slice_const_delta_standalone(job_number)
% To take slices through a phase diagram
% Script to build the flattened Liovillian and the diagonalise it:
% clear

if ischar(job_number), job_number = str2num(job_number); end;

save_prefix = 'two_site_JCH_phase_diagram_slice_const_delta_';
save_name = [save_prefix num2str(job_number) '.mat'];

global spre_a spost_a spre_ad spost_ad spre_ada spost_ada spre_a_spost_ad spre_sm_spost_sp spre_spsm spost_spsm spre_adsm spost_adsm spre_asp spost_asp spre_photon_hopping spost_photon_hopping spre_trans_photon_hopping spost_trans_photon_hopping

% Size of the system before shrinking:
M = 2;  
n_max = 4;

% Parameters:

J = 0.01;

omega_d = 10000;    % A large cavity frequency to find the two-particle modes

g = 1;

delta_sym = 0;  % For each different J we'll set delta = delta_sym + 2J

J_list = logspace(-3,1,43);

gamma_p = 0.01; % Cavity loss rate
gamma_a = 0;

Omega = 0.3*gamma_p;   % Cavity driving strength

% Set the variable which is changed internally in this file:
var_list = J_list;

% Load in the shrinking data...
shrinking_data = load('shrinking_data_for_M_2_n_max_4_N_excs_8');

M_shrink = shrinking_data.M_shrink;
M_shrink_perm = shrinking_data.M_shrink_perm;
add_matrix = shrinking_data.add_matrix;
connected_elements_store = shrinking_data.connected_elements_store;
non_repeated_elements = shrinking_data.non_repeated_elements;

dim_rho_after_first_shrink = size(M_shrink,2);

% Define the shrunken system operators using the shrink data:
[spre_a, spost_a, spre_ad, spost_ad, spre_ada, spost_ada, spre_a_spost_ad, spre_sm_spost_sp, spre_spsm, spost_spsm, spre_adsm, spost_adsm, spre_asp, spost_asp, spre_photon_hopping, spost_photon_hopping, spre_trans_photon_hopping, spost_trans_photon_hopping] = build_shrunken_system_operators(M_shrink, add_matrix, M_shrink_perm, M, n_max);

% Initialise store matrices:
num_1_store = zeros(1,length(var_list));
numsq_1_store = zeros(1,length(var_list));
ee_1_store = zeros(1,length(var_list));
corr_store = zeros(1,length(var_list));

np_na_1_store = zeros(1,length(var_list));

fock_1_store = zeros(n_max + 1,length(var_list));

% Define operators:
a = a_n_N_particle_sector(n_max);
ad = a_dag_n_N_particle_sector(n_max);
num = ad*a;
num_1 = tensor_matrix(num, M, 1);
num_2 = tensor_matrix(num, M, 2);

num_1_shrink = M_shrink.'*num_1*M_shrink;
num_2_shrink = M_shrink.'*num_2*M_shrink;

numsq = ad*a*ad*a;
numsq_1 = tensor_matrix(numsq, M, 1);
numsq_1_shrink = M_shrink.'*numsq_1*M_shrink;

for loop_n = 0:n_max
    
    kn = [zeros(loop_n,1); 1; zeros(n_max - loop_n,1)];
    proj_n = kron(eye(2),kn*kn');
    proj_n(n_max+1,:) = [];
    proj_n(:,n_max+1) = [];
    proj_n = kron(proj_n, eye(2*n_max+1));
    proj_n_1_shrink{loop_n+1} = M_shrink.'*proj_n*M_shrink;
    
end

sp = sp_N_particle_sector(n_max);
sm = sm_N_particle_sector(n_max);
N_ee = sp*sm;
N_ee_1 = tensor_matrix(N_ee, M, 1);

np_na = num*N_ee;
np_na_1 = tensor_matrix(np_na, M, 1);
np_na_1_shrink = M_shrink.'*np_na_1*M_shrink;

N_ee_1_shrink = M_shrink.'*N_ee_1*M_shrink;

corr_op = tensor_matrix(ad*a, M, 1)*tensor_matrix(ad*a, M, 2);
corr_op_shrink = M_shrink.'*corr_op*M_shrink;

warning off

for loop_var = 1:length(var_list)
    
    % Countdown:
    length(var_list) - loop_var
    
    % Set whatever needs to be set:
    J = J_list(loop_var);
    delta = delta_sym + 2*J;
    
    % Find the two-particle mode:
    H = jch_hamiltonian_multi_site_pbc(M, M, omega_d, omega_d-delta, g, J, 0); % Use n_max = M here...
    
    [v,d] = eig(full(H));
    
    evs = diag(d);
    
    % One particle states:
    one_particle_state_indices = find((evs > 0.95*omega_d) & (evs < 1.05*omega_d));
    one_particle_state_frequencies = evs(one_particle_state_indices);
    
    % Two particle states:
    two_particle_state_indices = find((evs > 0.95*omega_d*2) & (evs < 1.05*omega_d*2));
    two_particle_state_frequencies = evs(two_particle_state_indices);
    delta_omega_d = omega_d - two_particle_state_frequencies(1)/2;   % Choose the lowest lying two particle mode to drive
    delta_omega_0 = delta_omega_d - delta;
    
    L_shrink = build_flattened_liouvillian_in_reduced_space(delta_omega_d, delta_omega_0, g, J, Omega, gamma_p, gamma_a, M, n_max);
    
    disp('Diagonalising shrunken L...')
    [min_evector, min_evalue] = eigs(L_shrink,1,'SM');
    
    filled_rho = fill_out_shrunken_rho(dim_rho_after_first_shrink, min_evector, connected_elements_store, non_repeated_elements);
    
    num_1_store(loop_var) = trace(filled_rho*num_1_shrink);
    ee_1_store(loop_var) = trace(filled_rho*N_ee_1_shrink);
    numsq_1_store(loop_var) = trace(filled_rho*numsq_1_shrink);
    
    np_na_1_store(loop_var) = trace(filled_rho*np_na_1_shrink);
    
    corr_store(loop_var) = trace(filled_rho*corr_op_shrink);
   
    for loop_n = 0:n_max
        
        fock_1_store(loop_n + 1, loop_var) = trace(filled_rho*proj_n_1_shrink{loop_n+1});
    
    end
    
end

% Post processing quantities:
total_num_ss_store = num_1_store + ee_1_store
g2_ss_store = (numsq_1_store - num_1_store)./num_1_store.^2
var_total_num_ss_store = numsq_1_store + 2*np_na_1_store + ee_1_store

save(save_name, 'num_1_store', 'numsq_1_store', 'g2_ss_store', 'var_total_num_ss_store', 'var_list')
disp(['Saved data to ' save_name '... Quitting!'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spre_a_matrices, spost_a_matrices, spre_ad_matrices, spost_ad_matrices, spre_ada_matrices, spost_ada_matrices, spre_a_spost_ad_matrices, spre_sm_spost_sp_matrices, spre_spsm_matrices, spost_spsm_matrices, spre_adsm_matrices, spost_adsm_matrices, spre_asp_matrices, spost_asp_matrices, spre_photon_hopping_matrices, spost_photon_hopping_matrices, spre_trans_photon_hopping_matrices, spost_trans_photon_hopping_matrices] = build_shrunken_system_operators(M_shrink, add_matrix, M_shrink_perm, M, n_max)
% Generate all the relevant shrunken super-operators in our problem --
% these are shrunk via restriction to a maximum particle number subspace,
% and also via symmetry considerations

spre_a_matrices = [];
spost_a_matrices = [];
spre_ad_matrices = [];
spost_ad_matrices = [];
spre_ada_matrices = [];
spost_ada_matrices = [];
spre_a_spost_ad_matrices = [];
spre_sm_spost_sp_matrices = [];
spre_spsm_matrices = [];
spost_spsm_matrices = [];
spre_adsm_matrices = [];
spost_adsm_matrices = [];
spre_asp_matrices = [];
spost_asp_matrices = [];
spre_photon_hopping_matrices = [];
spost_photon_hopping_matrices = [];
spre_trans_photon_hopping_matrices = [];
spost_trans_photon_hopping_matrices = [];

% Single cavity operators:
a = kron(speye(2), spdiags(sqrt([0:n_max+1].'), 1, n_max+1, n_max+1));
ad = kron(speye(2), spdiags(sqrt([1:n_max+1].'), -1, n_max+1, n_max+1));

sp = kron([0 1; 0 0], speye(n_max+1));
sm = kron([0 0; 1 0], speye(n_max+1));

asp = a*sp;
adsm = ad*sm;

a(n_max+1,:) = [];
a(:,n_max+1) = [];

ad(n_max+1,:) = [];
ad(:,n_max+1) = [];

sp(n_max+1,:) = [];
sp(:,n_max+1) = [];

sm(n_max+1,:) = [];
sm(:,n_max+1) = [];

asp(n_max+1,:) = [];
asp(:,n_max+1) = [];

adsm(n_max+1,:) = [];
adsm(:,n_max+1) = [];

% First build the matrices whole:
for site = 1:M
   
    spre_a_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*(tensor_matrix(a, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spost_a_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*(tensor_matrix(a, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spre_ad_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*(tensor_matrix(ad, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spost_ad_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*(tensor_matrix(ad, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    
    spre_ada_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*(tensor_matrix(ad*a, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spost_ada_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*(tensor_matrix(ad*a, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    
    spre_a_spost_ad_matrices{site} = M_shrink_perm.'*(spre(M_shrink.'*(tensor_matrix(a, M, site))*M_shrink)*spost(M_shrink.'*(tensor_matrix(ad, M, site))*M_shrink))*add_matrix*M_shrink_perm;
    spre_sm_spost_sp_matrices{site} = M_shrink_perm.'*(spre(M_shrink.'*(tensor_matrix(sm, M, site))*M_shrink)*spost(M_shrink.'*(tensor_matrix(sp, M, site))*M_shrink))*add_matrix*M_shrink_perm;
    
    spre_spsm_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*(tensor_matrix(sp*sm, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spost_spsm_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*(tensor_matrix(sp*sm, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    
    spre_adsm_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*(tensor_matrix(adsm, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spost_adsm_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*(tensor_matrix(adsm, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    
    spre_asp_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*(tensor_matrix(asp, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    spost_asp_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*(tensor_matrix(asp, M, site))*M_shrink)*add_matrix*M_shrink_perm;
    
    spre_photon_hopping_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*tensor_matrix(ad, M, site)*tensor_matrix(a, M, mod(site,M)+1)*M_shrink)*add_matrix*M_shrink_perm; 
    spost_photon_hopping_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*tensor_matrix(ad, M, site)*tensor_matrix(a, M, mod(site,M)+1)*M_shrink)*add_matrix*M_shrink_perm; 
    
    spre_trans_photon_hopping_matrices{site} = M_shrink_perm.'*spre(M_shrink.'*tensor_matrix(a, M, site)*tensor_matrix(ad, M, mod(site,M)+1)*M_shrink)*add_matrix*M_shrink_perm; 
    spost_trans_photon_hopping_matrices{site} = M_shrink_perm.'*spost(M_shrink.'*tensor_matrix(a, M, site)*tensor_matrix(ad, M, mod(site,M)+1)*M_shrink)*add_matrix*M_shrink_perm; 

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = spre(A)
% Function to find the superoperator associated with premultiplying density
% matrix by A

% Get dimension of operator:
d = size(A,2);

S = kron(speye(d), A);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = spost(A)
% Function to find the superoperator associated with premultiplying density
% matrix by A

% Get dimension of operator:
d = size(A,2);

S = kron(A.', speye(d));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ad = a_dag_n_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>
ad = kron(speye(2), a_n_dag(n_max));
ad(n_max+1,:) = [];
ad(:,n_max+1) = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = a_n_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>

a = kron(speye(2), a_n(n_max));
a(n_max+1,:) = [];
a(:,n_max+1) = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function operator = a_n(n_max)

% Function to return the destruction operator in a Fock state basis with
% number states from 0 to n-1

operator = spdiags(sqrt([0:n_max+1].'), 1, n_max+1, n_max+1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function operator = a_n_dag(n_max)

% Function to return the creation operator in a Fock state basis with
% number states from 0 to n-1

operator = spdiags(sqrt([1:n_max+1].'), -1, n_max+1, n_max+1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp = sp_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>

sp = kron([0 1;0 0], speye(n_max+1));
sp(n_max+1,:) = [];
sp(:,n_max+1) = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sm = sm_N_particle_sector(n_max)
% returns the full photonic raising operator in a basis that has at most
% n_max excitations, i.e. we exclude the ket |e,n_max>

sm = kron([0 0;1 0], speye(n_max+1));
sm(n_max+1,:) = [];
sm(:,n_max+1) = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L_shrink = build_flattened_liouvillian_in_reduced_space(delta_omega_d, delta_omega_0, g, A, Omega, gamma_p, gamma_a, M, n_max)
% most of the parameters are

global spre_a spost_a spre_ad spost_ad spre_ada spost_ada spre_a_spost_ad spre_sm_spost_sp spre_spsm spost_spsm spre_adsm spost_adsm spre_asp spost_asp spre_photon_hopping spost_photon_hopping spre_trans_photon_hopping spost_trans_photon_hopping

dim_shrunk_L = size(spre_a{1},1);

L_shrink = sparse(dim_shrunk_L,dim_shrunk_L);   % Initialise shrunken Liouvillian

for m = 1:M     % Loop over the lattice sites
    
    % Hamiltonian part:
    L_shrink = L_shrink + 1/1i*(delta_omega_d*(spre_ada{m} - spost_ada{m}) + ...
        delta_omega_0*(spre_spsm{m} - spost_spsm{m}) + ...
        g*(spre_adsm{m} - spost_adsm{m} + spre_asp{m} - spost_asp{m}) + ...
         Omega*(spre_a{m} - spost_a{m} + spre_ad{m} - spost_ad{m}) - ...
         A*(spre_photon_hopping{m} - spost_photon_hopping{m} + spre_trans_photon_hopping{m} - spost_trans_photon_hopping{m}));
    
     L_shrink = L_shrink + gamma_p/2*(2*spre_a_spost_ad{m} - spre_ada{m} - spost_ada{m});
     L_shrink = L_shrink + gamma_a/2*(2*spre_sm_spost_sp{m} - spre_spsm{m} - spost_spsm{m});
     
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filled_rho = fill_out_shrunken_rho(full_dim, min_evector, connected_elements_store, non_repeated_elements)

filled_rho = sparse(full_dim, full_dim);

% Fill in the subset of values that were explicitly solved for:
filled_rho(non_repeated_elements) = min_evector;
num_connected_paths = length(connected_elements_store);

for loop = 1:num_connected_paths
    
    if ne(length(connected_elements_store{loop}),1) % If this path is not just an element connected to itself
       
        % Add the columns of L corresponding to the redundant elements:
        elements_to_fill_this_connection = connected_elements_store{loop}(2:end);
        filled_rho(elements_to_fill_this_connection) = filled_rho(connected_elements_store{loop}(1));
        
    else
        
        continue
        
    end
    
end

filled_rho = filled_rho/trace(filled_rho);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%5