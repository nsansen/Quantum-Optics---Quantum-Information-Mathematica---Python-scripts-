clear
clc

N0_list_lossless = [1 2 3 4 7 10 20 40 50];
N0_list_lossy = [1 2 3 4 7 10];

g2_time_averaged_lossless = NaN*ones(length(N0_list_lossless), 101);
g2_cross_time_averaged_lossless = NaN*ones(length(N0_list_lossless), 101);
g2_time_averaged_lossy = NaN*ones(length(N0_list_lossy), 101);
g2_cross_time_averaged_lossy = NaN*ones(length(N0_list_lossy), 101);


Z_time_averaged_store_lossless = NaN*ones(length(N0_list_lossless), 101);
Z_time_averaged_store_lossy = NaN*ones(length(N0_list_lossy), 101);

for loop_N0 = 1:length(N0_list_lossless)
    
    N0 = N0_list_lossless(loop_N0);
    
%     try
        data_lossless = load(['./data_from_cluster_nikos_time_evolve_lossless/nikos_2_site_lossless_increasing_N0_' num2str(loop_N0) '.mat']);
%     catch me
%         continue
%     end
    
    g2_t_store = (data_lossless.numsq_t_store - data_lossless.num_t_store)./data_lossless.num_t_store.^2;
    
    g2_time_averaged_lossless(loop_N0, :) = (mean(g2_t_store(:,1:end-1),2));
    
%     try
        
        Z_t_store_lossless = (data_lossless.num_t_store - data_lossless.num_R_t_store)./(data_lossless.num_t_store + data_lossless.num_R_t_store);
        
        Z_time_averaged_store_lossless(loop_N0, :) = (mean(Z_t_store_lossless(:,1:end-1),2));
        
%         g2_cross_t_store = data_lossless.densdens_cross_t_store./(data_lossless.num_t_store.*data_lossless.num_R_t_store);
        g2_cross_t_store = data_lossless.densdens_cross_t_store./(data_lossless.num_t_store + data_lossless.num_R_t_store).^2;
        
        g2_cross_time_averaged_lossless(loop_N0, :) = (mean(g2_cross_t_store(:,1:end-1),2));
        
%     catch me
%         continue
%     end
    
end

% Z_time_averaged_store_lossless = mean(Z_t_store_lossless, 2);

for loop_N0 = 1:length(N0_list_lossy)
    
    N0 = N0_list_lossy(loop_N0);
    
    try
        data_lossy = load(['./data_from_cluster_nikos_time_evolve/nikos_2_site_lossless_increasing_N0_' num2str(loop_N0) '.mat'])
%         data_lossy = load(['./data_from_cluster_nikos_time_evolve/nikos_2_site_lossless_increasing_N0_lower_loss' num2str(loop_N0) '.mat'])

        g2_cross_t_store_lossy = data_lossy.densdens_cross_t_store./(data_lossy.num_t_store.*data_lossy.num_R_t_store);
        
        g2_cross_time_averaged_lossless(loop_N0, :) = (mean(g2_cross_t_store(:,1:end-1),2));
        
    catch me
        continue
    end
    
    g2_t_store = (data_lossy.numsq_t_store - data_lossy.num_t_store)./data_lossy.num_t_store.^2;
    
    size(g2_t_store)
    
    g2_time_averaged_lossy(loop_N0, :) = (mean(g2_t_store(:,1:end-1),2));
    
    Z_t_store_lossy = (data_lossy.num_t_store - data_lossy.num_R_t_store)./(data_lossy.num_t_store + data_lossy.num_R_t_store);    
    
    Z_time_averaged_store_lossy(loop_N0, :) = (mean(Z_t_store_lossy(:,1:end-1),2));
    
end

%% Plotting
figure(101)
clf
plot(linspace(0, 2.5, 101), Z_time_averaged_store_lossy.')
legend('N_0 = 1', '2', '3', '4', '7', '10')
set(gca, 'FontSize', 16)
xlabel('g / g_c')
ylabel('\bar{Z}')
title('Photon imbalance for lossy case for inreasing initial photon number')

figure(102)
clf
plot(linspace(0, 2.5, 101), g2_time_averaged_lossy.')
legend('N_0 = 1', '2', '3', '4', '7', '10')
set(gca, 'FontSize', 16)
xlabel('g / g_c')
ylabel('N0')
title('Time averaged g2 for lossy case for inreasing initial photon number')

figure(103)
clf
plot(linspace(0, 2.5, 101), Z_time_averaged_store_lossless.')
legend('N_0 = 1', '2', '3', '4', '7', '10', '20', '40', '50')
set(gca, 'FontSize', 16)
xlabel('g / g_c')
ylabel('\bar{Z}')
title('Photon imbalance for lossless case for inreasing initial photon number')

figure(104)
clf
plot(linspace(0, 2.5, 101), g2_time_averaged_lossless.')
legend('N_0 = 1', '2', '3', '4', '7', '10', '20', '40', '50')
set(gca, 'FontSize', 16)
xlabel('g / g_c')
ylabel('N0')
title('Time averaged g2 for lossless case for inreasing initial photon number')