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