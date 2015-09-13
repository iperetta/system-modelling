PID = input('Enter with the desired PID: ');

A = dlmread(['GP_' num2str(PID) 'bestfit_vec.txt']);
B = dlmread(['GP_' num2str(PID) 'meanfit_vec.txt']);
close all
plot(A(5:end),'-k','LineWidth',2)
hold on
plot(B(5:end),'--k','LineWidth',2)
legend('Best fitness', 'Mean fitness (population)')
xlabel('# GENERATIONS','FontSize',14)
ylabel('FITNESS','FontSize',14)
hold off
set(gca,'FontSize',14)
grid on