% Main program to execute disaggregation of load curve

% Aggregated daily load curve

clear all
clc

AggreLoadCurve = [0.310, 0.300, 0.297, 0.296, 0.300, 0.351, 0.433, 0.507, 0.477, 0.523, 0.543, 0.516, 0.479, 0.480, 0.463, 0.489, 0.536, 0.526, 0.507, 0.500, 0.520, 0.470, 0.407, 0.333];
Hrs = 1:24;
Ratings = [2,1,7];
NumLoad = 2;
IndiLoad = DisaggregLoad(AggreLoadCurve,NumLoad,0.5);

subplot(2,1,1)
for i = 1:1:NumLoad
    plot(Hrs,IndiLoad(:,i));hold on;
end
grid on
set(gca,'FontSize',20)
xlabel('Hours')
xlim([1,24])
ylabel('Individual Load Curves','FontSize',20)
% legend('Individual Load-1','Individual Load-2','Individual Load-3')


subplot(2,1,2)
plot(Hrs,AggreLoadCurve,'r-*','LineWidth',4);hold on;
plot(Hrs,sum(IndiLoad,2),'b-^','LineWidth',2)
legend('Aggregate Load','Sum of Individual Loads')
grid on
set(gca,'FontSize',20)
xlabel('Hours')
xlim([1,24])
ylabel('Comparision of Load Curves')




