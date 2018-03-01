load('Data.mat')

Time = data(1,:);
Travel = data(2,:);
Travel_rate = data(3,:);
Pitch = data(4,:);
Pitch_rate = data(5,:);
Elevation = data(6,:);
Elevation_rate = data(7,:);

h = figure(1);
subplot(211)
stairs(t,u1),grid
ylabel('u_d (V)')
subplot(212)
stairs(t,u2),grid
ylabel('u_s (V)')
xlabel('Time (s)')

f = figure(2);
subplot(611)
plot(t, x1, 'b', t, x1, 'bo', Time, Travel', 'r'),grid
ylabel('lambda (rad)')
subplot(612)
plot(t, x2, 'b', t, x2', 'bo', Time, Travel_rate, 'r'),grid
ylabel('r (rad)')
subplot(613)
plot(t, x3, 'b', t, x3, 'bo', Time, Pitch, 'r'),grid
ylabel('p (rad)')
subplot(614)
plot(t, x4, 'b', t, x4', 'bo', Time, Pitch_rate, 'r'),grid
ylabel('p\_dot (rad)')
subplot(615)
plot(t, x5, 'b', t, x5', 'bo', Time, Elevation, 'r'),grid
ylabel('e (rad)')
subplot(616)
plot(t, x6, 'b', t, x6' ,'bo', Time, Elevation_rate, 'r'),grid
ylabel('e\_dot (rad)')
xlabel('Time (s)')

set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(h, 'Inputs_CL', '-dpdf', '-r0')

set(f, 'Units', 'Centimeters');
pos = get(f, 'Position');
set(f, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(f, 'States_CL', '-dpdf', '-r0')