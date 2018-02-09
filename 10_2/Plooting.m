load('Travel.mat')

hold on
plot(e_measured(1,:),e_measured(2,:),'mo')
plot(t,x2.*180/pi,'m',t,x2'.*180/pi,'mo')
ylabel('lambda')
hold off