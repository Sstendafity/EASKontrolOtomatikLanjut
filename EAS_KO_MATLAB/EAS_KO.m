clc
close all

%% Inisiasi Variabel

bw = 15020;
bc = 350;
kw = 500000;
kc = 80000;
mc = 2500;
mw = 320;

%% Matriks State Space 
A = [0 1 0 -1; -kc/mc -bc/mc 0 bc/mc; 0 0 0 1; kc/mw bc/mw -kw/mw (-bc-bw)/mw];
B = [0 0; 1/mc 0; 0 -1; -1/mw bw/mw];
C = [1 0 0 0];
D = 0;

%% Respons sistem tanpa LQR

sys = ss(A,B,C,D);
figure(1)
step(0.1*sys);           % Step response
poles = pole(sys);
f = stepinfo(0.1*sys);    % Parameter performansi system

%% LQR dan Grafik Respons

Q = eye(4);             % Matriks identitas 4x4, karena terdapat state variabel 4
R = eye(2);             % Matriks identitas 2x2, karena terdapat 2 input

[K,S,P] = lqr(A,B,Q,R); % LQR [gain, solusi persamaan Riccati S, pole]

sys1 = ss(A-B*K,B,C,D);  % Sistem dengan menggunakan gain (k) dari LQR
figure(2)
step(0.1*sys1)             % Step response sistem dengan LQR
f1 = stepinfo(0.1*sys1);      % Parameter performansi system

%% Analisis sistem dengan menggunakan variasi matriks Q dengan matriks R tetap

figure(3)
step(0.1*sys1)                     %Variasi 1
disp('Variasi Q ke-1')
display(K);
display(S);
display(P);
display(struct2table(f1))
hold on 

for i=2:5
    Q1 = i*Q;                    %Variasi Q
    [K1,S1,P1] = lqr(A,B,Q1,R); 
    disp(['Variasi Q ke-', num2str(i)])
    display(K1);
    display(S1);
    display(P1);
    sys2 = ss(A-B*K1,B,C,D);  
    step(0.1*sys2)             
    f2 = stepinfo(0.1*sys2);
    display(struct2table(f2))
end 

legend('Variasi Q ke-1', 'Variasi Q ke-2', 'Variasi Q ke-3', 'Variasi Q ke-4', 'Variasi Q ke-5')

%% Analisis sistem dengan menggunakan variasi matriks R dengan matriks Q tetap

figure(4)
step(0.1*sys1)                     %Variasi 1
disp('Variasi R ke-1')
display(K);
display(S);
display(P);
display(struct2table(f1))
hold on 

for i=2:5
    R1 = i*R;                    %Variasi Q
    [K1,S1,P1] = lqr(A,B,Q,R1); 
    disp(['Variasi R ke-', num2str(i)])
    display(K1);
    display(S1);
    display(P1);
    sys2 = ss(A-B*K1,B,C,D);  
    step(0.1*sys2)             
    f2 = stepinfo(0.1*sys2);
    display(struct2table(f2))
end 

legend('Variasi R ke-1', 'Variasi R ke-2', 'Variasi R ke-3', 'Variasi R ke-4', 'Variasi R ke-5')

%% Analisis Kestabilan dengan Lyapunov

disp('Uji Kestabilan Sistem Tanpa LQR dengan Lyapunov')
P = lyapunov(A);            % Kestabilan sistem sebelum menggunakan LQR

disp('Uji Kestabilan Sistem dengan LQR dengan Lyapunov')
P_LQR = lyapunov(A-B*K);    % Kestabilan sistem setelah menggunakan LQR

function P = lyapunov(A)
Q_lya = eye(4);
P = lyap(A,Q_lya);

if all(eig(P)>0)
    disp('The system is stable')
    display(eig(P))
else
    disp('The system is unstable')
    display(P)
end
end