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




plot(yu,(-3+3*yu)/0.9,'c')
plot(yu,(30-3*yu)/0.9,'c')

plot(yu,(-3+3*yu)/1.1,'m')
plot(yu,(30-3*yu)/1.1,'m')

xlabel('yu')
ylabel('yl')
xticks(0:1:11)
yticks(0:1:35)
ylim([0 inf])
grid on
saveas(gcf,'ToyExample2RobustLowerLevel.jpg')


%% Two dimensional toy problem (not bilevel) 
x=linspace(0,9,100);
figure(5)
plot(x,(30+25*x)/(20),'k')
hold on
xlim([0 9])
xticks(0:1:9)
ylim([0 9])
yticks(0:1:9)
grid on

plot(x,(10-x)/(2),'k')
plot(x,(15-2*x)/(-1),'k')
plot(x,(-15+2*x)/(-10),'k')

plot(x,(30+24.9*x)/(20.1),'m')
plot(x,(10-1.1*x)/(2.1),'m')
plot(x,(15-2.1*x)/(-0.9),'m')
plot(x,(-15+1.9*x)/(-9.9),'m')

plot(x,(30+25.1*x)/(19.9),'c')
plot(x,(10-0.9*x)/(1.9),'c')
plot(x,(15-1.9*x)/(-1.1),'c')
plot(x,(-15+2.1*x)/(-10.1),'c')

ylabel('y')
xlabel('x')

saveas(gcf,'ToyExample1Robust1D.jpg')

%% Toy Problem 1 (not bilevel) with 2d uncertainty

