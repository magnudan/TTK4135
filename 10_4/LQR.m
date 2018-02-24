Q1 = 1; %travel
Q2 = 1; %travel rate
Q3 = 0.1; %pitch
Q4 = 1; %pitch rate
Q5 = 100; %elevation
Q6 = 1; %elevation rate
R1 = 1; %pitch input
R2 = 1; %elevation input
Q_LQR = diag([Q1,Q2,Q3,Q4,Q5,Q6]);  %gain for states
R_LQR = diag([R1,R2]);           %gain for inputs

[K_LQR ,Shit_Mat_LQR ,eig_LQR] = dlqr(A1,B1,Q_LQR,R_LQR);

K_LQR_T = K_LQR.';