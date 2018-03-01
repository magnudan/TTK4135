load('Data.mat')

Time = data(1,:);
Travel = data(2,:)+pi;
Travel_rate = data(3,:);
Pitch = data(4,:);
Pitch_rate = data(5,:);
Elevation = data(6,:);
Elevation_rate = data(7,:);

h = figure(1);
subplot(511)
stairs(t,u),grid
ylabel('u (V)')
subplot(512)
plot(t, x1, 'b', t, x1, 'bo', Time, Travel', 'r'),grid
ylabel('lambda (rad)')
subplot(513)
plot(t, x2, 'b', t, x2', 'bo', Time, Travel_rate, 'r'),grid
ylabel('r (rad)')
subplot(514)
plot(t, x3, 'b', t, x3, 'bo', Time, Pitch, 'r'),grid
ylabel('p (rad)')
subplot(515)
plot(t, x4, 'b', t, x4', 'bo', Time, Pitch_rate, 'r'),grid
ylabel('p\_dot (rad)')
xlabel('Time (s)')

set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(h, '10_2_Q_is_0,2', '-dpdf', '-r0')