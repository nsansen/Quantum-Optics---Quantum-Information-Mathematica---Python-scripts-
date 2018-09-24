% Script to load data for PBC BHM

file_prefix = 'BH_ED_drive_1_part_res_increase_J_compare_iTEBD_';
% file_prefix = 'BH_ED_linear_spectrum_OBC_';

start_job = 1;
end_job = 1;

numjobs = end_job - start_job + 1;

% close all
figure_to_use = 20;
figure(figure_to_use)

numjobs_loaded = 0;

track_feature_freq_store = zeros(1, numjobs);

J_list = logspace(-2, 1, 51);
var_list = J_list;

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
            
            num_ss_store = zeros(numjobs,length(var_list));
            numsq_ss_store = zeros(numjobs,length(var_list));
            corr_ss_store = zeros(numjobs,length(var_list));
            g2_ss_store = zeros(numjobs,length(var_list));
            
        end
        
        for loop_var = 1:length(var_list)
            
            rho = data.rho_store(:,:,loop_var);
            num_ss_store(loop, loop_var) = trace(rho*ad_1*a_1);
            numsq_ss_store(loop, loop_var) = trace(rho*ad_1*a_1*ad_1*a_1);
            corr_ss_store(loop, loop_var) = trace(rho*ad_1*a_2);
            
        end
        
    end 
    
end

% Post-processing:
g2 = (numsq_ss_store - num_ss_store)./num_ss_store.^2;

figure(30)
% plot(log10(J_list/0.1), num_ss_store,'m')
plot(linspace(-3,3,51), num_ss_store,'m.')