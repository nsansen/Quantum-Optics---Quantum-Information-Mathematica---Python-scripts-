% Script to load data for PBC BHM

tot_num_jobs = 21;

file_prefix = 'two_site_spectrum_JCH_increasing_delta_';

start_job = 1;
end_job = tot_num_jobs;

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
        
        %     % Find the lowest-lying two-particle mode for this job:
        %     omega_d = 10000;    % A large cavity frequency to find the two-particle modes
        %     H = bhm_hamiltonian_multi_site_pbc(data.M, data.M, omega_d, data.J, data.U, 0); % Use n_max = M here...
        %
        %     [v,d] = eig(full(H));
        %
        %     evs = diag(d);
        %
        %     % One particle states:
        %     one_particle_state_indices = find((evs > 0.95*omega_d) & (evs < 1.05*omega_d));
        %     one_particle_state_frequencies = evs(one_particle_state_indices);
        %
        %     % Two particle states:
        %     two_particle_state_indices = find((evs > 0.95*omega_d*2) & (evs < 1.05*omega_d*2));
        %     two_particle_state_frequencies = evs(two_particle_state_indices);
        %     delta_omega_d_track = omega_d - two_particle_state_frequencies(1)/2;   % Choose the lowest lying two particle mode to drive
        %     track_feature_freq_store(loop) = delta_omega_d_track;
        
        % Plotting:
        %     figure(figure_to_use)
        hold on
        g2_plot = g2_ss_store(loop,:);
        plot(-data.var_list, log10(real(data.g2_ss_store)),'g','LineWidth',2)
        plot(-data.var_list, log10(data.num_1_store),'k','LineWidth',2)
        
        
        % xlim([-.75 1.25])
        %     plot([-.75 1.25], [0 0], 'k:')
        xlabel('\Delta_c / U', 'FontSize', 14)
        ylabel('log_{10} NESS exp. vals.', 'FontSize', 14)
        
        
        y_lim = ylim;
        
        %     plot(-delta_omega_d_track/(2*data.U)*[1 1], y_lim)
         
    end
    
end