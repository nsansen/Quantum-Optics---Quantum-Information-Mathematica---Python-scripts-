function nikos_time_evolve_2_sites(job_number)
% Script to build the flattened Liovillian and the diagonalise it:
% clear

global L

if ischar(job_number), job_number = str2num(job_number); end;

save_prefix = 'nikos_2_site_lossless_increasing_N0_';
save_name = [save_prefix num2str(job_number) '.mat'];

% Parameters:

M = 2  % Number of sites in ring

J = 1;  % Coherent coupling betweenr resonators

N0_list = [1 2 3 4 7 10];

N0 = N0_list(job_number); % Initial number of photons in the left resonator

n_max = N0;   % Photon truncation parameter

D_H = 2*(n_max+1) - 1;  % Local dimension of Hilbert space

dim_H = (2*(n_max+1)-1)^M;

g_c = 2.8*sqrt(N0)*J;   % The semi-classical critical atom-resonator coupling

delta = 0;

delta_omega_d = 0;
delta_omega_0 = delta_omega_d - delta;

gamma_p = 0.05; % Cavity loss rate
gamma_a = 0.05;    % Atom loss rate

Omega = 0;   % Cavity driving strength

g_list = linspace(0, 2.5, 11)*g_c;

% Time range for integration:
t = linspace(0, 20*J, 251);

% Set the variable which is changed internally in this file:
var_list = g_list;

% Initial condition
b = a_n(n_max+1); % Destruction operator for one site
b_dag = a_n_dag(n_max+1); % Creation operator

% Set up identity operators for the atomic and photonic parts
eye_a = eye(2);
eye_p = eye(n_max+1);

b = kron(eye_a, b);
b_dag = kron(eye_a, b_dag);
num = b_dag*b;    % Number operator

num(n_max+1,:) = [];
num(:,n_max+1) = [];

b(n_max+1,:) = [];
b(:,n_max+1) = [];

b_dag(n_max+1,:) = [];
b_dag(:,n_max+1) = [];

num_L = kron(num, eye(D_H));
num_R = kron(eye(D_H), num);

sigma_p = kron([0 1;0 0], eye_p);
sigma_m = kron([0 0;1 0], eye_p);
num_a = sigma_p*sigma_m;

sigma_p(n_max+1,:) = [];
sigma_p(:,n_max+1) = [];

sigma_m(n_max+1,:) = [];
sigma_m(:,n_max+1) = [];

D_H = 2*(n_max+1) - 1;  % Local dimension of Hilbert space
y00_blank = [1;zeros(D_H - 1,1)];
y00_L = sigma_m*y00_blank;
for n = 1:N0    
    y00_L = b_dag*y00_L;
end
y00_L = y00_L/sqrt(y00_L'*y00_L);

y00_R = sigma_m*y00_blank;

y00 = kron(y00_L, y00_R);

rho00 = y00*y00';

rho_store = zeros(dim_H, dim_H, length(var_list));

num_t_store = zeros(length(var_list), length(t));
numsq_t_store = zeros(length(var_list), length(t));
num_R_t_store = zeros(length(var_list), length(t));
densdens_cross_t_store = zeros(length(var_list), length(t));

for loop = 1:length(var_list)
    
    disp(num2str(length(var_list) - loop))
    
    % Set what needs to be set:
    g = g_list(loop);
    
%     H = jch_hamiltonian_multi_site_pbc(M, n_max, delta_omega_d, delta_omega_0, g, J, Omega);
    H = jch_hamiltonian_multi_site_obc(M, n_max, delta_omega_d, delta_omega_0, g, J, Omega);

    L = flattened_liouvillian(H, n_max, gamma_p, gamma_a, M);
    
    rho0_flat = reshape(rho00, [D_H^4, 1]);
    
    for loop_t = 1:length(t)-1
        
        length(t)-1 - loop_t
        
        t_int = [t(loop_t) t(loop_t+1)];
                
        [TS, YS] = ode45(@masterkernel, t_int, rho0_flat);
        
        rho0_flat = YS(end,:).';
        
        rho_t = reshape(rho0_flat, [D_H^2, D_H^2]);
        
        num_t_store(loop, loop_t) = trace(rho_t*num_L);
        numsq_t_store(loop, loop_t) = trace(rho_t*num_L*num_L);
        num_R_t_store(loop, loop_t) = trace(rho_t*num_R);
        densdens_cross_t_store(loop, loop_t) = trace(rho_t*num_R*num_L);
        
    end
    
%     rho_ss = get_steady_state_rho(M, n_max, delta_omega_d, delta_omega_0, g, J, Omega, gamma_p, gamma_a);
    
end

save(save_name, 'n_max', 'M', 'gamma_p', 'gamma_a', 'Omega', 'J', 'g', 'num_t_store', 'numsq_t_store', 'num_R_t_store', 'densdens_cross_t_store')
disp(['Saved data to file ' save_name])

end

%%%%%%%%%%%%%%%%5 Auxillary functions %%%%%%%%%%%%%%%%%%%%%%%
function rho_ss = get_steady_state_rho(M, n_max, delta_omega_d, delta_omega_0, g, J, Omega, gamma_p, gamma_a)
% Function to return the steady state given the system parameters

H = jch_hamiltonian_multi_site_pbc(M, n_max, delta_omega_d, delta_omega_0, g, J, Omega);

L = flattened_liouvillian(H, n_max, gamma_p, gamma_a, M);

[min_evector, min_evalue] = eigs(L,1,'SM');
dim = sqrt(length(min_evector));
rho_ss = reshape(min_evector, dim, dim);
rho_ss = rho_ss/trace(rho_ss);

end
%%%%%%%%%%%%%%%
function L = flattened_liouvillian(H_unitary, n_max, gamma_p, gamma_a, M)
% Function to return the `flattened' Liouvillian operator, size N^2 times
% N^2, to premultiply a flattened density vector (size N^2 times 1)
% H_unitary is the Hermitian part Hamiltonian
% n_max is photon truncation parameter
% gamma is the loss rate from each cavity
% M is the system size

% Construct the creation and annhilation operators for the field:
c = kron(speye(2), a_n(n_max+1));
c(n_max+1,:) = [];
c(:,n_max+1) = [];
c_dag = kron(speye(2), a_n_dag(n_max+1));
c_dag(n_max+1,:) = [];
c_dag(:,n_max+1) = [];

sp = kron([0 1;0 0], speye(n_max+1));
sm = kron([0 0;1 0], speye(n_max+1));

sp(n_max+1,:) = [];
sp(:,n_max+1) = [];
sm(n_max+1,:) = [];
sm(:,n_max+1) = [];

L = 1/i*(spre(H_unitary) - spost(H_unitary));    % The Hermitian commutator

for m = 1:M     % Loop over the lattice sites
    
    c_m = tensor_matrix(c, M, m);
    c_dag_m = tensor_matrix(c_dag, M, m);
    
    sp_m = tensor_matrix(sp, M, m);
    sm_m = tensor_matrix(sm, M, m);
    
    L = L + gamma_p/2*(2*spre(c_m)*spost(c_dag_m) - spre(c_dag_m*c_m) - spost(c_dag_m*c_m));    % The Liouvillian for the field
    L = L + gamma_a/2*(2*spre(sm_m)*spost(sp_m) - spre(sp_m*sm_m) - spost(sp_m*sm_m));    % The Liouvillian for the atom
    
end

end
%%%%%%%%%%%%%%%%%%%%%
function H = jch_hamiltonian_multi_site_pbc(M, n_max, delta_omega_d, delta_omega_0, g, A, Omega)
% Function to return the BH Hamiltonian for OBC
% M - number of lattice sites
% n_max - boson truncation
% delta_omega_d - difference between boson frequency and driving frequency:
% delta_omega_d = omega_boson - omega_laser
% delta_omega_0 - difference between atomic frequency and driving frequency:
% n_max (maximum photon number at which we truncate the
% Fock space basis at each site), A (hopping energy)

D_H = 2*(n_max+1) - 1;

b = a_n(n_max+1); % Destruction operator for one site
b_dag = a_n_dag(n_max+1); % Creation operator

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

bd_sm = b_dag*sigma_m;
b_sp = b*sigma_p;

% Now take out the |e,n_max> state adter we've made all our definitions
b(n_max+1,:) = [];
b(:,n_max+1) = [];

b_dag(n_max+1,:) = [];
b_dag(:,n_max+1) = [];

num(n_max+1,:) = [];
num(:,n_max+1) = [];

sigma_p(n_max+1,:) = [];
sigma_p(:,n_max+1) = [];

sigma_m(n_max+1,:) = [];
sigma_m(:,n_max+1) = [];

num_a(n_max+1,:) = [];
num_a(:,n_max+1) = [];

bd_sm(n_max+1,:) = [];
bd_sm(:,n_max+1) = [];

b_sp(n_max+1,:) = [];
b_sp(:,n_max+1) = [];

dim_H = (2*(n_max+1)-1)^M;
H = sparse(dim_H,dim_H);   % Initialise Hamiltonian

for m = 1:M-1     % Loop over the lattice sites
    
    b_m = tensor_matrix(b, M, m);
    b_dag_m = tensor_matrix(b_dag, M, m);
    b_m_plus_1 = tensor_matrix(b, M, m+1);
    b_dag_m_plus_1 = tensor_matrix(b_dag, M, m+1);
    bd_sm_m = tensor_matrix(bd_sm, M, m);
    b_sp_m = tensor_matrix(b_sp, M, m);
    num_m = tensor_matrix(num, M, m);
    num_a_m = tensor_matrix(num_a, M, m);
    sigma_p_m = tensor_matrix(sigma_p, M, m);
    sigma_m_m = tensor_matrix(sigma_m, M, m);
    
    hopping_term = b_dag_m*b_m_plus_1 + b_dag_m_plus_1*b_m;
    
    H = H + delta_omega_d*num_m + delta_omega_0*num_a_m - A*hopping_term + g*(bd_sm_m + b_sp_m) + Omega*(b_m + b_dag_m);   % Assemble the Hamiltonian for each lattice site contribution
    
end

% Need to also include the last on-site term:
b_M = tensor_matrix(b, M, M);
b_dag_M = tensor_matrix(b_dag, M, M);
num_M = tensor_matrix(num, M, M);
num_a_M = tensor_matrix(num_a, M, M);
sigma_p_M = tensor_matrix(sigma_p, M, M);
sigma_m_M = tensor_matrix(sigma_m, M, M);
bd_sm_M = tensor_matrix(bd_sm, M, M);
b_sp_M = tensor_matrix(b_sp, M, M);

H = H + delta_omega_d*num_M + delta_omega_0*num_a_M + g*(bd_sm_M + b_sp_M) + Omega*(b_M + b_dag_M);   % Add the last site's on-site interaction contribution

% Hopping term:
b_1 = tensor_matrix(b, M, 1);
b_dag_1 = tensor_matrix(b_dag, M, 1);
hopping_term = b_1*b_dag_M + b_dag_1*b_M;

H = H - A*hopping_term;

end
%%%%%%%%%%%%%%%%%%%%%
function H = jch_hamiltonian_multi_site_obc(M, n_max, delta_omega_d, delta_omega_0, g, A, Omega)
% Function to return the BH Hamiltonian for OBC
% M - number of lattice sites
% n_max - boson truncation
% delta_omega_d - difference between boson frequency and driving frequency:
% delta_omega_d = omega_boson - omega_laser
% delta_omega_0 - difference between atomic frequency and driving frequency:
% n_max (maximum photon number at which we truncate the
% Fock space basis at each site), A (hopping energy)

D_H = 2*(n_max+1) - 1;

b = a_n(n_max+1); % Destruction operator for one site
b_dag = a_n_dag(n_max+1); % Creation operator

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

bd_sm = b_dag*sigma_m;
b_sp = b*sigma_p;

% Now take out the |e,n_max> state adter we've made all our definitions
b(n_max+1,:) = [];
b(:,n_max+1) = [];

b_dag(n_max+1,:) = [];
b_dag(:,n_max+1) = [];

num(n_max+1,:) = [];
num(:,n_max+1) = [];

sigma_p(n_max+1,:) = [];
sigma_p(:,n_max+1) = [];

sigma_m(n_max+1,:) = [];
sigma_m(:,n_max+1) = [];

num_a(n_max+1,:) = [];
num_a(:,n_max+1) = [];

bd_sm(n_max+1,:) = [];
bd_sm(:,n_max+1) = [];

b_sp(n_max+1,:) = [];
b_sp(:,n_max+1) = [];

dim_H = (2*(n_max+1)-1)^M;
H = sparse(dim_H,dim_H);   % Initialise Hamiltonian

for m = 1:M-1     % Loop over the lattice sites
    
    b_m = tensor_matrix(b, M, m);
    b_dag_m = tensor_matrix(b_dag, M, m);
    b_m_plus_1 = tensor_matrix(b, M, m+1);
    b_dag_m_plus_1 = tensor_matrix(b_dag, M, m+1);
    bd_sm_m = tensor_matrix(bd_sm, M, m);
    b_sp_m = tensor_matrix(b_sp, M, m);
    num_m = tensor_matrix(num, M, m);
    num_a_m = tensor_matrix(num_a, M, m);
    sigma_p_m = tensor_matrix(sigma_p, M, m);
    sigma_m_m = tensor_matrix(sigma_m, M, m);
    
    hopping_term = b_dag_m*b_m_plus_1 + b_dag_m_plus_1*b_m;
    
    H = H + delta_omega_d*num_m + delta_omega_0*num_a_m - A*hopping_term + g*(bd_sm_m + b_sp_m) + Omega*(b_m + b_dag_m);   % Assemble the Hamiltonian for each lattice site contribution
    
end

% Need to also include the last on-site term:
b_M = tensor_matrix(b, M, M);
b_dag_M = tensor_matrix(b_dag, M, M);
num_M = tensor_matrix(num, M, M);
num_a_M = tensor_matrix(num_a, M, M);
sigma_p_M = tensor_matrix(sigma_p, M, M);
sigma_m_M = tensor_matrix(sigma_m, M, M);
bd_sm_M = tensor_matrix(bd_sm, M, M);
b_sp_M = tensor_matrix(b_sp, M, M);

H = H + delta_omega_d*num_M + delta_omega_0*num_a_M + g*(bd_sm_M + b_sp_M) + Omega*(b_M + b_dag_M);   % Add the last site's on-site interaction contribution

end
%%%%%%%%%%%%%%%%5
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drho_flat_shrink = masterkernel(t,rho_flat_shrink)
% ODE kernel of the spins MF master equation

global L

drho_flat_shrink = L*rho_flat_shrink;

end