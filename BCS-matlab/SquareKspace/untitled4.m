clear all
clc

count=100;
a=rand(count,3);
figure
hold on
for x=1:count
quiver3(0,0,0,a(x,1),a(x,2),a(x,3),0,'Color',[a(x,1),a(x,2),a(x,3)])
end
hold off
axis equal