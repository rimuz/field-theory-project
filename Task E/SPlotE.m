close all
clc

% Plot of S parameters for dominant modes

c = 299792458;

% Sizes of the apertures
a1 = 15e-2;  % 1.5cm
a2 = pi*10*1e-2; % 3.14 cm
a3 = 7.5e-2;    % 1cm

% Location of apertures
xoff1 = 1e-2;      % 1 cm
xoff2 = 0;         % 0 cm
xoff3 = 5e-3;      % 0.5 cm

% Length of the middle waveguide
delta_z = 3e-2;     % 30 cm

w = (14.99e9:10e6:50e9)*2*pi;
f = w /2/pi;
k = w/c;

n_t1 = pi/a1;
n_t2 = pi/a2;
n_t3 = pi/a3;

fprintf("Resonance frequency a2: %d", (c/2*sqrt((n_t2/pi)^2 + (1/delta_z)^2)));

S11_L = zeros(length(w), 1);
S12_L = zeros(length(w), 1);
S21_L = zeros(length(w), 1);
S22_L = zeros(length(w), 1);
S11_R = zeros(length(w), 1);
S12_R = zeros(length(w), 1);
S21_R = zeros(length(w), 1);
S22_R = zeros(length(w), 1);
S11_LR = zeros(length(w), 1);
S12_LR = zeros(length(w), 1);
S21_LR = zeros(length(w), 1);
S22_LR = zeros(length(w), 1);

for i = 1:length(w)
    n_z1 = ( ( k(i) >= n_t1 ) - 1i*( k(i) < n_t1 ) ) .* sqrt( abs( k(i)^2 - n_t1.^2 ) );
    n_z2 = ( ( k(i) >= n_t2 ) - 1i*( k(i) < n_t2 ) ) .* sqrt( abs( k(i)^2 - n_t2.^2 ) );
    n_z3 = ( ( k(i) >= n_t3 ) - 1i*( k(i) < n_t3 ) ) .* sqrt( abs( k(i)^2 - n_t3.^2 ) );
    
    [S11_L(i), S12_L(i), S21_L(i), S22_L(i)] = aperture_S(n_t1, n_t2, n_z1, n_z2, 1, 1, k(i), k(i), xoff1, a1, a2);
    P = exp(-1i*n_z2*delta_z);
    [S22_R(i), S12_R(i), S21_R(i), S11_R(i)] = aperture_S(n_t3, n_t2, n_z3, n_z2, 1, 1, k(i), k(i), xoff3, a3, a2);
    [S11_PR, S12_PR, S21_PR, S22_PR] = combineLR(0, P, P, 0, S11_R(i), S12_R(i), S21_R(i), S22_R(i));
    [S11_LR(i), S12_LR(i), S21_LR(i), S22_LR(i)] = combineLR(S11_L(i), S12_L(i), S21_L(i), S22_L(i), S11_PR, S12_PR, S21_PR, S22_PR);
end

figure(1)
tiledlayout(3,2);
nexttile;
plot(f, abs(S11_L), "b");
title("S11L");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S12_L), "b");
title("S12L");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S21_L));
title("S21L");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S22_L));
title("S22L");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S11_L).^2 + abs(S12_L).^2);
title("S11^2 + S12^2");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S22_L).^2 + abs(S21_L).^2);
title("S22^2 + S21^2");
xlabel("Frequency [Hz]");


figure(2)
tiledlayout(3,2);
nexttile;
plot(f, real(S11_R), "b");
title("S11R");
xlabel("Frequency [Hz]");

nexttile;
plot(f, real(S12_R), "b");
title("S12R");
xlabel("Frequency [Hz]");

nexttile;
plot(f, real(S21_R));
title("S21R");
xlabel("Frequency [Hz]");

nexttile;
plot(f, real(S22_R));
title("S22R");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S11_R).^2 + abs(S12_R).^2);
title("S11^2 + S21^2");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S22_R).^2 + abs(S21_R).^2);
title("S22^2 + S21^2");
xlabel("Frequency [Hz]");

figure(3);
tiledlayout(3,2);
nexttile;
plot(f, abs(S11_LR), "b");
hold on;
plot(f, imag(S11_LR), "r");
title("S11");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S12_LR), "b");
hold on;
plot(f, imag(S12_LR), "r");
title("S12");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S21_LR), "b");
hold on;
plot(f, imag(S21_LR), "r");
title("S21");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S22_LR), "b");
hold on;
plot(f, imag(S11_LR), "r");
title("S22");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S11_LR).^2 + abs(S12_LR).^2);
title("S11^2 + S21^2");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S22_LR).^2 + abs(S21_LR).^2);
title("S22^2 + S21^2");
xlabel("Frequency [Hz]");