clc
clear
close all

data_low = load('two_site_spectrum_JCH_small_hopping_1');
data_high = load('two_site_spectrum_JCH_high_hopping_1');

figure(1)
clf
hold on
plot(-data_low.var_list, log10(real(data_low.num_1_store)),'b', 'LineWidth', 2)
% plot(-data.var_list, log10(real(data_low.)),'r', 'LineWidth', 2)
plot(-data_low.var_list, log10(real(data_low.g2_ss_store)),'g', 'LineWidth', 2)
xlabel('\Delta_c / g', 'FontSize', 14)
ylabel('log_{10} NESS exp vals', 'FontSize', 14)
xlim([-1.5 1.5])
ylim([-4 2])
x_lim = xlim;
y_lim = ylim;
plot(x_lim, 0*[1 1], 'k:')
set(gcf, 'Color', 'w')
export_fig 'two_site_spectrum_JCH_small_loss_small_hopping' '-pdf'

figure(2)
clf
hold on
plot(-(data_high.var_list - 2*10), log10(real(data_high.num_1_store)),'b', 'LineWidth', 2)
% plot(-data.var_list, log10(real(data_high.)),'r', 'LineWidth', 2)
plot(-(data_high.var_list - 2*10), log10(real(data_high.g2_ss_store)),'g', 'LineWidth', 2)
xlabel('\Delta_c / g', 'FontSize', 14)
ylabel('log_{10} NESS exp vals', 'FontSize', 14)
xlim([-1.5 1.5])
ylim([-4 2])
x_lim = xlim;
y_lim = ylim;
plot(x_lim, 0*[1 1], 'k:')
set(gcf, 'Color', 'w')
export_fig 'two_site_spectrum_JCH_small_loss_high_hopping' '-pdf'