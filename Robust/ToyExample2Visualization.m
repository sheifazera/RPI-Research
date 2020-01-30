A=[-2,3;1,1];
B=[-3,1;3,1]; %Lower level constraints
b=[-3;30]; %Lower level constraint RHS
d=[12;14];

P=0.1*[1,0;0,1]; %WHERE DO WE GET P? (a+Pxi)^Tx<=d, norm(xi,2)<=1
xi1=[1;0];
xi2=[0;1];
xi1n=[-1;0];
xi2n=[0;-1];

yu=linspace(1,10,100);
%yl=linspace(0,10,100);

figure()
plot(yu,-A(1,1)/A(1,2)*yu+d(1)/A(1,2),'k')
hold on
grid on
xlabel('yu')
ylabel('yl')
plot(yu,-A(2,1)/A(2,2)*yu+d(2)/A(2,2),'k')
plot(yu,-B(1,1)/B(1,2)*yu+b(1)/B(1,2),'k') % Lower level
plot(yu,-B(2,1)/B(2,2)*yu+b(2)/B(2,2),'k') % Lower level
%%
plot(yu,-(A(1,1)+P(1,:)*xi1)/(A(1,2)+P(2,:)*xi1)*yu+d(1)/(A(1,2)+P(2,:)*xi1),'-m')
plot(yu,-(A(2,1)+P(1,:)*xi1)/(A(2,2)+P(2,:)*xi1)*yu+d(2)/(A(2,2)+P(2,:)*xi1),'-c')

plot(yu,-(A(1,1)+P(1,:)*xi2)/(A(1,2)+P(2,:)*xi2)*yu+d(1)/(A(1,2)+P(2,:)*xi2),'--m')
plot(yu,-(A(2,1)+P(1,:)*xi2)/(A(2,2)+P(2,:)*xi2)*yu+d(2)/(A(2,2)+P(2,:)*xi2),'--c')

plot(yu,-(A(1,1)+P(1,:)*xi1n)/(A(1,2)+P(2,:)*xi1n)*yu+d(1)/(A(1,2)+P(2,:)*xi1n),'-mo')
plot(yu,-(A(2,1)+P(1,:)*xi1n)/(A(2,2)+P(2,:)*xi1n)*yu+d(2)/(A(2,2)+P(2,:)*xi1n),'-co')

plot(yu,-(A(1,1)+P(1,:)*xi2n)/(A(1,2)+P(2,:)*xi2n)*yu+d(1)/(A(1,2)+P(2,:)*xi2n),'-mx')
plot(yu,-(A(2,1)+P(1,:)*xi2n)/(A(2,2)+P(2,:)*xi2n)*yu+d(2)/(A(2,2)+P(2,:)*xi2n),'-cx')