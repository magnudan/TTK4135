Q_LQR = diag([10,0.1,10,0.1]);  %gain for states
R_LQR = 1;                     %gain for inputs

[K_LQR ,Shit_Mat_LQR ,eig_LQR] = dlqr(A1,B1,Q_LQR,R_LQR);

K_LQR_T = K_LQR.';