load('Data.mat')

Time = data(1,:);
Travel = data(2,:);
Travel_rate = data(3,:);
Pitch = data(4,:);
Pitch_rate = data(5,:);
%Elevation = data(6,:);
%Elevation_rate = data(7,:);

aa = figure(1);
stairs(t,u),grid
ylabel('u_d (V)')
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

set(aa, 'Units', 'Centimeters');
pos = get(aa, 'Position');
set(aa, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(aa, '10_3_Inputs', '-dpdf', '-r0')

set(bb, 'Units', 'Centimeters');
pos = get(bb, 'Position');
set(bb, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(bb, '10_3_Travel', '-dpdf', '-r0')

set(cc, 'Units', 'Centimeters');
pos = get(cc, 'Position');
set(cc, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', 'PaperSize', [pos(3), pos(4)])
print(cc, '10_3_Pitch', '-dpdf', '-r0')