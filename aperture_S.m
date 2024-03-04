
% Calculate S parameters for waveguide with height a1 feeding waveguide
% with height a2, a1 < a2. x1 is the offset of the aperture, n_t1 and n_t2
% the eigenvalues at the two waveguides, M and N are the number of modes.
function [S11,S12,S21,S22] = aperture_S(n_t1, n_t2, M, N, k1, k2, x1, a1, a2) 
    n_z1 = sqrt(k1^2 - n_t1.^2);
    n_z2 = sqrt(k2^2 - n_t2.^2);
    x2 = x1 + a1;

    C = zeros(M, M);
    for m = 1:M
        for n = 1:N
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
    S21 = 2 .* inv(C.'*conj(C) + eye(N)) * C.';
    S22 = inv(C.'*conj(C) + eye(N)) * (C.'*conj(C) - eye(N));
end
