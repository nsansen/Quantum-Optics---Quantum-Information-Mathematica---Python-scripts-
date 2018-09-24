% Script to load data for PBC BHM

tot_num_jobs = 15;

file_prefix = 'two_site_spectrum_BH_increasing_U_';
% file_prefix = 'two_site_spectrum_BH_increasing_J_';

start_job = 8;
end_job = 8;

numjobs = end_job - start_job + 1;

close all
figure_to_use = 20;
figure(figure_to_use)

numjobs_loaded = 0;

track_feature_freq_store = zeros(1, numjobs);

for loop = 1:numjobs
    
    job_to_load = start_job + loop - 1;
    
    file_name = [file_prefix num2str(job_to_load) '.mat'];
    data_loaded = 0;
    
    try
        
        disp(['Loading file ' file_name])
        data = load(file_name);
        data_loaded = 1;
        numjobs_loaded = numjobs_loaded + 1;
        
    catch me
        disp(['Error! File ' file_name ' not loaded!'])
    end
    
    if data_loaded == 1
        
        if numjobs_loaded == 1
            
            M = data.M;
            n_max = data.n_max;
            a = spdiags(sqrt([0:n_max+1].'), 1, n_max+1, n_max+1);
            ad = spdiags(sqrt([1:n_max].'), -1, n_max+1, n_max+1);
            a_1 = tensor_matrix(a,M,1);
            ad_1 = tensor_matrix(ad,M,1);
            a_2 = tensor_matrix(a,M,2);
            ad_2 = tensor_matrix(ad,M,2);
            
%             [a_j_matrices, ad_j_matrices]  = define_momentum_space_operators(M, n_max);
            
            var_list = data.delta_omega_d_list;
%             var_list = data.var_list;
            
            num_ss_store = zeros(numjobs,length(var_list));
            numsq_ss_store = zeros(numjobs,length(var_list));
            corr_ss_store = zeros(numjobs,length(var_list));
            g2_ss_store = zeros(numjobs,length(var_list));
            
%             mom_mode_occupations_store = zeros(numjobs, length(var_list), M);
            
        end
        
        for loop_var = 1:length(var_list)
            
            rho = data.rho_store(:,:,loop_var);
            num_ss_store(loop, loop_var) = trace(rho*ad_1*a_1);
            numsq_ss_store(loop, loop_var) = trace(rho*ad_1*a_1*ad_1*a_1);
            corr_ss_store(loop, loop_var) = trace(rho*ad_1*a_2);
            
            if loop_var == 151
                
                pause_here = 1
                
            end
            
%             for mom_mode = 1:M
%                 
%                 mom_mode_occupations_store(loop, loop_var, mom_mode) = trace(rho*ad_j_matrices(:,:,mom_mode)*a_j_matrices(:,:,mom_mode));
%             
%             end
            
        end
        
    end
    
    % Find the lowest-lying two-particle mode for this job:
    omega_d = 10000;    % A large cavity frequency to find the two-particle modes
    H = bhm_hamiltonian_multi_site_pbc(data.M, data.M, omega_d, data.J, data.U, 0); % Use n_max = M here...
    
    [v,d] = eig(full(H));
    
    evs = diag(d);
    
    % One particle states:
    one_particle_state_indices = find((evs > 0.95*omega_d) & (evs < 1.05*omega_d));
    one_particle_state_frequencies = evs(one_particle_state_indices);
    
    % Two particle states:
    two_particle_state_indices = find((evs > 0.95*omega_d*2) & (evs < 1.05*omega_d*2));
    two_particle_state_frequencies = evs(two_particle_state_indices);
    delta_omega_d_track = omega_d - two_particle_state_frequencies(1)/2;   % Choose the lowest lying two particle mode to drive
    track_feature_freq_store(loop) = delta_omega_d_track;
    
    g2_ss_store = (numsq_ss_store - num_ss_store)./num_ss_store.^2;
    
    % Plotting:
    figure(figure_to_use)
    hold on
    g2_plot = g2_ss_store(loop,:);
    
        plot(-data.delta_omega_d_list/(2*data.U), log10(num_ss_store(loop,:)),'k','LineWidth',2)
    
%     figure(10)
%     hold on
    plot(-data.delta_omega_d_list/(2*data.U), log10(real(g2_plot)),'g','LineWidth',2)
    
    % xlim([-.75 1.25])
%     plot([-.75 1.25], [0 0], 'k:')
    xlabel('\Delta_c / U', 'FontSize', 14)
    ylabel('log_{10} NESS exp. vals.', 'FontSize', 14)
    
    
    y_lim = ylim;
        
    plot(-delta_omega_d_track/(2*data.U)*[1 1], y_lim)
        
    
end