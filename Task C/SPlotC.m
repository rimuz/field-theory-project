close all
clc

% Plot of S parameters for dominant modes

c = 299792458;

% Sizes of the apertures
a1 = 15.47e-2;  % 1.5cm
a2 = 75.47*sqrt(1.0001)*1e-2; % 3.14 cm
a3 = 15e-2;    % 1cm
a4 = 19e-2;

% Location of apertures
xoff1 = 30e-2;       % 1 cm
xoff2 = 0;           % 0 cm
xoff3 = 7.5e-2;      % 0.5 cm
xoff4 = 30e-2;

% Length of middle cavity
delta_z = 37.5e-2; % 3cm

w = (14.99e9:10e6:20e9)*2*pi;
f = w /2/pi;
k = w/c;

n_t1 = pi/a1;
n_t2 = pi/a2;
n_t3 = pi/a3;
n_t4 = pi/a4; 

S11_LR = zeros(length(w), 1);
S12_LR = zeros(length(w), 1);
S13_LR = zeros(length(w), 1);
S21_LR = zeros(length(w), 1);
S22_LR = zeros(length(w), 1);
S23_LR = zeros(length(w), 1);
S31_LR = zeros(length(w), 1);
S32_LR = zeros(length(w), 1);
S33_LR = zeros(length(w), 1);

for i = 1:length(w)
    n_z1 = ( ( k(i) >= n_t1 ) - 1i*( k(i) < n_t1 ) ) .* sqrt( abs( k(i)^2 - n_t1.^2 ) );
    n_z2 = ( ( k(i) >= n_t2 ) - 1i*( k(i) < n_t2 ) ) .* sqrt( abs( k(i)^2 - n_t2.^2 ) );
    n_z3 = ( ( k(i) >= n_t3 ) - 1i*( k(i) < n_t3 ) ) .* sqrt( abs( k(i)^2 - n_t3.^2 ) );
    n_z4 = ( ( k(i) >= n_t4 ) - 1i*( k(i) < n_t4 ) ) .* sqrt( abs( k(i)^2 - n_t4.^2 ) );

    [S11_L, S12_L, S21_L, S22_L] = aperture_S(n_t1, n_t2, n_z1, n_z2, 1, 1, k(i), k(i), xoff1, a1, a2);
    P = exp(-1i*n_z2*delta_z);
    [S22_R, S21_R, S12_R, S11_R] = two_apertures_S(n_t3, n_t2, n_t4, n_z3, n_z2, n_z4, 1, 1, 1, k(i), k(i), k(i), xoff3, xoff4, a3, a2, a4);
    [S11_PR, S12_PR, S21_PR, S22_PR] = combineLR(0, P, P, 0, S11_R, S12_R, S21_R, S22_R);
    [S11, S12, S21, S22] = combineLR(S11_L, S12_L, S21_L, S22_L, S11_PR, S12_PR, S21_PR, S22_PR);
    
    S11_LR(i) = S11;
    S12_LR(i) = S12(1, 1);
    S13_LR(i) = S12(1, 2);
    S21_LR(i) = S21(1, 1);
    S22_LR(i) = S22(1, 1);
    S23_LR(i) = S22(1, 2);
    S31_LR(i) = S21(2, 1);
    S32_LR(i) = S22(2, 1);
    S33_LR(i) = S22(2, 2);

end

tiledlayout(4,3);
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
plot(f, abs(S13_LR), "b");
hold on;
plot(f, imag(S13_LR), "r");
title("S13");
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
plot(f, imag(S22_LR), "r");
title("S22");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S23_LR), "b");
hold on;
plot(f, imag(S23_LR), "r");
title("S23");
xlabel("Frequency [Hz]");


nexttile;
plot(f, abs(S31_LR), "b");
hold on;
plot(f, imag(S31_LR), "r");
title("S31");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S32_LR), "b");
hold on;
plot(f, imag(S32_LR), "r");
title("S32");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S33_LR), "b");
hold on;
plot(f, imag(S33_LR), "r");
title("S33");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S11_LR).^2 + abs(S12_LR).^2 + abs(S13_LR).^2);
title("S11^2 + S12^2 + S13^2");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S22_LR).^2 + abs(S21_LR).^2 + abs(S23_LR).^2);
title("S21^2 + S22^2 + S23^2");
xlabel("Frequency [Hz]");

nexttile;
plot(f, abs(S31_LR).^2 + abs(S32_LR).^2 + abs(S33_LR).^2);
title("S31^2 + S32^2 + S33^2");
xlabel("Frequency [Hz]");