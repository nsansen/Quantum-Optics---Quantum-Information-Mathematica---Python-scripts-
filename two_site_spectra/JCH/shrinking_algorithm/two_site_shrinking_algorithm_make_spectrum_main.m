clear
clc

global spre_a spost_a spre_ad spost_ad spre_ada spost_ada spre_a_spost_ad spre_sm_spost_sp spre_spsm spost_spsm spre_adsm spost_adsm spre_asp spost_asp spre_photon_hopping spost_photon_hopping spre_trans_photon_hopping spost_trans_photon_hopping

% Size of the system before shrinking:
M = 2;  
n_max = 4;

% Parameters:

J = 0.001;

g = 1;

delta = 2*J;

gamma_p = 0.05; % Cavity loss rate
gamma_a = 0;

Omega = 0.3*gamma_p;   % Cavity driving strength

delta_omega_d_list = linspace(-1.5*g + 2*J, 1.5*g + 2*J,101);

% Set the variable which is changed internally in this file:
var_list = delta_omega_d_list;

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

num_deloc_store = zeros(1,length(var_list));
numsq_deloc_store = zeros(1,length(var_list));
ee_deloc_store = zeros(1,length(var_list));

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

% Number operator for the delocalised mode:
a_deloc = 1/sqrt(3)*(tensor_matrix(a,M,1) + tensor_matrix(a,M,2) + tensor_matrix(a,M,3));
ad_deloc = 1/sqrt(3)*(tensor_matrix(ad,M,1) + tensor_matrix(ad,M,2) + tensor_matrix(ad,M,3));
num_deloc = ad_deloc*a_deloc;
num_deloc_shrink = M_shrink.'*num_deloc*M_shrink;
numsq_deloc = ad_deloc*a_deloc*ad_deloc*a_deloc;
numsq_deloc_shrink = M_shrink.'*numsq_deloc*M_shrink;

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

% Delocalised atomic modes:
sm_deloc = 1/sqrt(3)*(tensor_matrix(sm,M,1) + tensor_matrix(sm,M,2) + tensor_matrix(sm,M,3));
sp_deloc = 1/sqrt(3)*(tensor_matrix(sp,M,1) + tensor_matrix(sp,M,2) + tensor_matrix(sp,M,3));

ee_deloc = sp_deloc*sm_deloc;
ee_deloc_shrink = M_shrink.'*ee_deloc*M_shrink;

corr_op = tensor_matrix(ad*a, M, 1)*tensor_matrix(ad*a, M, 2);
corr_op_shrink = M_shrink.'*corr_op*M_shrink;

warning off

for loop_var = 1:length(var_list)
    
    % Countdown:
    length(var_list) - loop_var
    
    % Set whatever needs to be set:
    delta_omega_d = var_list(loop_var);
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
    
    num_deloc_store(loop_var) = trace(filled_rho*num_deloc_shrink);
    numsq_deloc_store(loop_var) = trace(filled_rho*numsq_deloc_shrink);
    ee_deloc_store(loop_var) = trace(filled_rho*ee_deloc_shrink);
   
    for loop_n = 0:n_max
        
        fock_1_store(loop_n + 1, loop_var) = trace(filled_rho*proj_n_1_shrink{loop_n+1});
    
    end
    
end

%% Post processing quantities:
total_num_ss_store = num_1_store + ee_1_store
g2_ss_store = (numsq_1_store - num_1_store)./num_1_store.^2
var_total_num_ss_store = numsq_1_store + 2*np_na_1_store + ee_1_store