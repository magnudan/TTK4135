load('Data.mat')

Time = data(1,:);
Travel = data(2,:);
Travel_rate = data(3,:);
Pitch = data(4,:);
Pitch_rate = data(5,:);
Elevation = data(6,:);
Elevation_rate = data(7,:);

aa = figure(1);
subplot(211)
stairs(t,u1),grid
ylabel('u_d (V)')
subplot(212)
stairs(t,u2),grid
ylabel('u_s (V)')
xlabel('Time (s)')

bb = figure(2);
subplot(211)
plot(t, x1, 'b', t, x1, 'bo', Time, Travel', 'r'),grid
ylabel('lambda (rad)')
subplot(212)
plot(t, x2, 'b', t, x2', 'bo', Time, Travel_rate, 'r'),grid
ylabel('r (rad)')
xlabel('Time (s)')

cc = figure(3);
subplot(211)
plot(t, x3, 'b', t, x3, 'bo', Time, Pitch, 'r'),grid
ylabel('p (rad)')
subplot(212)
plot(t, x4, 'b', t, x4', 'bo', Time, Pitch_rate, 'r'),grid
ylabel('p\_dot (rad)')
xlabel('Time (s)')

dd = figure(4);
subplot(211)
plot(t, x5, 'b', t, x5', 'bo', Time, Elevation, 'r'),grid
ylabel('e (rad)')
subplot(212)
plot(t, x6, 'b', t, x6' ,'bo', Time, Elevation_rate, 'r'),grid
ylabel('e\_dot (rad)')
xlabel('Time (s)')

set(aa, 'Units', 'Centimeters');
pos = get(aa, 'Position');
set(aa, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(aa, 'Inputs_OL', '-dpdf', '-r0')

set(bb, 'Units', 'Centimeters');
pos = get(bb, 'Position');
set(bb, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(bb, 'Travel_OL', '-dpdf', '-r0')

set(cc, 'Units', 'Centimeters');
pos = get(cc, 'Position');
set(cc, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(cc, 'Pitch_OL', '-dpdf', '-r0')

set(dd, 'Units', 'Centimeters');
pos = get(dd, 'Position');
set(dd, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(dd, 'Elevation_OL', '-dpdf', '-r0')