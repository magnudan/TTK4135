load('Data.mat')

Time = data(1,:);
Travel = data(2,:);
Travel_rate = data(3,:);
Pitch = data(4,:);
Pitch_rate = data(5,:);
Elevation = data(6,:);
Elevation_rate = data(7,:);

f = figure(2);
subplot(611)
plot(Time, Travel', 'r'),grid
ylabel('lambda (rad)')
subplot(612)
plot(Time, Travel_rate, 'r'),grid
ylabel('r (rad)')
subplot(613)
plot(Time, Pitch, 'r'),grid
ylabel('p (rad)')
subplot(614)
plot(Time, Pitch_rate, 'r'),grid
ylabel('p\_dot (rad)')
subplot(615)
plot(Time, Elevation, 'r'),grid
ylabel('e (rad)')
subplot(616)
plot(Time, Elevation_rate, 'r'),grid
ylabel('e\_dot (rad)')
xlabel('Time (s)')

set(f, 'Units', 'Centimeters');
pos = get(f, 'Position');
set(f, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(f, 'States10_1', '-dpdf', '-r0')