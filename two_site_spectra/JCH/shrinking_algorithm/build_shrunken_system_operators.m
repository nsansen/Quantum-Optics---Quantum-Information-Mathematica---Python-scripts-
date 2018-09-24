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