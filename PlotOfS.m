close all
clc

% Plot of S parameters for dominant modes

c = 299792458;

% Sizes of the apertures
a1 = 1.5e-2;  % 1.5cm
a2 = pi*1e-2; % 3.14 cm
a3 = 1e-2;    % 1cm

% Location of apertures
xoff1 = 1e-2;      % 1 cm
xoff2 = 0;         % 0 cm
xoff3 = 5e-3;      % 0.5 cm

% Length of the middle waveguide
delta_z = 3e-2;     % 3 cm

w = 0:10e7:100e9;
k = w/c;

n_t1 = pi/a1;
n_t2 = pi/a2;
n_t3 = pi/a3;

S11_LR = zeros(length(w), 1);
S12_LR = zeros(length(w), 1);
S21_LR = zeros(length(w), 1);
S22_LR = zeros(length(w), 1);

for i = 1:length(w)
    [S11_L, S12_L, S21_L, S22_L] = aperture_S(n_t1, n_t2, 1, 1, k(i), k(i), xoff1, a1, a2);
    P = exp(-1j*n_t2*delta_z);
    [S22_R, S12_R, S21_R, S11_R] = aperture_S(n_t3, n_t2, 1, 1, k(i), k(i), xoff3, a3, a2);
    [S11_PR, S12_PR, S21_PR, S22_PR] = combineLR(0, P, P, 0, S11_R, S12_R, S21_R, S22_R);
    [S11_LR(i), S12_LR(i), S21_LR(i), S22_LR(i)] = combineLR(S11_L, S12_L, S21_L, S22_L, S11_PR, S12_PR, S21_PR, S22_PR);
end

tiledlayout(3,2);
nexttile;
plot(w, real(S11_LR), "b");
title("S11");
xlabel("Frequency [Hz]");

nexttile;
plot(w, real(S12_LR), "b");
title("S12");
xlabel("Frequency [Hz]");

nexttile;
plot(w, real(S21_LR));
title("S21");
xlabel("Frequency [Hz]");

nexttile;
plot(w, real(S22_LR));
title("S22");
xlabel("Frequency [Hz]");

nexttile;
plot(w, abs(S11_LR).^2 + abs(S21_LR).^2);
title("S11^2 + S21^2");
xlabel("Frequency [Hz]");

nexttile;
plot(w, abs(S22_LR).^2 + abs(S21_LR).^2);
title("S22^2 + S21^2");
xlabel("Frequency [Hz]");