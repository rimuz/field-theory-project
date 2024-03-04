close all
clc

% Plot of S parameters for dominant modes

c = 299792458;

% Sizes of the apertures
a1 = pi/2*1e-2;    % 1.571 cm
a2 = 3e-2;         % 3 cm
a3 = pi/3*1e-2;    % 1.047 cm

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
    [S11_L, S12_L, S21_L, S22_L] = aperture_S(n_t1, n_t2, 1, k(i), k(i), xoff1, a1, a2);
    P = exp(-1j*n_t2*delta_z);
    [S22_R, S12_R, S21_R, S11_R] = aperture_S(n_t3, n_t2, 1, k(i), k(i), xoff3, a3, a2);
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


% Calculate S parameters for waveguide with height a1 feeding waveguide
% with height a2, a1 < a2. x1 is the offset of the aperture, n_t1 and n_t2
% the eigenvalues at the two waveguides, M is the number of modes.
function [S11,S12,S21,S22] = aperture_S(n_t1, n_t2, M, k1, k2, x1, a1, a2) 
    n_z1 = sqrt(k1^2 - n_t1.^2);
    n_z2 = sqrt(k2^2 - n_t2.^2);
    x2 = x1 + a1;

    C = zeros(M, M);
    for m = 1:M
        for n = 1:M
            n_p = n_t1(m) + n_t2(n);
            n_m = n_t1(m) - n_t2(n);
            
            C(m,n) = conj(n_z2(n)) / sqrt(abs(n_z1(m))*abs(n_z2(n))*a1*a2) * ...
                (sin(n_m*x2 - n_t1(m)*x1) + sin(n_t2(n)*x1))/n_m - ....
                (sin(n_p*x2 - n_t1(m)*x1) - sin(n_t2(n)*x1))/n_p;
        end
    end
    
    % Compendium page 121
    % Medium is vacuum => Z real => D1 and D2 matrices are identity
    S11 = inv(eye(M) + conj(C)*C.') * (eye(M) - conj(C)*C.');
    S12 = 2 .* inv(eye(M) + conj(C)*C.') * conj(C);
    S21 = 2 .* inv(C.'*conj(C) + eye(M)) * C.';
    S22 = inv(C.'*conj(C) + eye(M)) * (C.'*conj(C) - eye(M));
end

% Compendium page 125
function [S11_LR, S13_LR, S31_LR, S33_LR] = ...
        combineLR(S11_L, S21_L, S12_L, S22_L, S22_R, S23_R, S32_R, S33_R)
    
    N = length(S11_L);
    S11_LR = S11_L + S12_L*inv(eye(N) - S22_R*S22_L)*S22_R*S21_L;
    S13_LR = S12_L*inv(eye(N) - S22_R*S22_L)*S23_R;
    S31_LR = S32_R*inv(eye(N) - S22_L*S22_R)*S21_L;
    S33_LR = S33_R + S32_R*inv(eye(N) - S22_L*S22_R)*S22_L*S23_R;
end