yu=linspace(1,11,100);
figure()
plot(yu,(12+2*yu)/3,'k--')
hold on
plot(yu,(14-yu),'k--')

plot(yu,(12+1.9*yu)/3,'m')
plot(yu,(12+2.1*yu)/3,'c')

plot(yu,(14-0.9*yu),'m')
plot(yu,(14-1.1*yu),'c')

plot(yu,(-3+3*yu),'k')
plot(yu,30-3*yu,'k')

xlabel('yu')
ylabel('yl')
ylim([0 inf])

xticks(0:1:11)
yticks(0:1:30)
grid on
hold on

saveas(gcf,'ToyExample2RobustUpperLevel.jpg')

%% Robust Lower Level

yu=linspace(0,12,100);
figure()

plot(yu,(12+2*yu)/3,'k')
hold on

plot(yu,(14-yu),'k')

plot(yu,(-3+3*yu),'k--')
plot(yu,30-3*yu,'k--')

plot(yu,(-3+3.1*yu),'c')
plot(yu,30-2.9*yu,'c')

plot(yu,(-3+2.9*yu),'m')
plot(yu,30-3.1*yu,'m')

xlabel('yu')
ylabel('yl')
xticks(0:1:11)
yticks(0:1:35)
ylim([0 inf])
grid on
saveas(gcf,'ToyExample2RobustLowerLevel.jpg')