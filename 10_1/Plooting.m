load('Data.mat')

Time = data(1,:);
Travel = data(2,:);
Travel_rate = data(3,:);
Pitch = data(4,:);
Pitch_rate = data(5,:);
Elevation = data(6,:);
Elevation_rate = data(7,:);

bb = figure(2);
subplot(211)
plot(Time, Travel', 'r'),grid
ylabel('lambda (rad)')
subplot(212)
plot(Time, Travel_rate, 'r'),grid
ylabel('r (rad)')
xlabel('Time (s)')

cc = figure(3);
subplot(211)
plot(Time, Pitch, 'r'),grid
ylabel('p (rad)')
subplot(212)
plot(Time, Pitch_rate, 'r'),grid
ylabel('p\_dot (rad)')
xlabel('Time (s)')

dd = figure(4);
subplot(211)
plot(Time, Elevation, 'r'),grid
ylabel('e (rad)')
subplot(212)
plot(Time, Elevation_rate, 'r'),grid
ylabel('e\_dot (rad)')
xlabel('Time (s)')

set(bb, 'Units', 'Centimeters');
pos = get(bb, 'Position');
set(bb, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(bb, '10_1_Travel', '-dpdf', '-r0')

set(cc, 'Units', 'Centimeters');
pos = get(cc, 'Position');
set(cc, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(cc, '10_1_Pitch', '-dpdf', '-r0')

set(dd, 'Units', 'Centimeters');
pos = get(dd, 'Position');
set(dd, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(dd, '10_1_Elevation', '-dpdf', '-r0')