%% grid, box, background color
grid on; box on
set(gcf,'color','w');

%% figure size

%% set labels
ylabel('Reactive Power [kVAr]');
xlabel('Distance [kft]');
%%
ylabel('Active Power [kW]');
xlabel('Distance [kft]');

%% fontsize
fs = 14;
set(gca,'fontsize',fs);
set(get(gca,'xlabel'),'fontsize',14);
set(get(gca,'ylabel'),'fontsize',14);

%% limits
ylim([-2000 1000]);