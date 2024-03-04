% Calculate S parameters for waveguide with height a1 feeding waveguide
% with height a2, a1 < a2. x1 is the offset of the aperture, n_t1 and n_t2
% the eigenvalues at the two waveguides, M and N are the number of modes.
function [S11,S12,S21,S22] = two_apertures_S(n_t1, n_t2, n_t3, n_z1, n_z2, n_z3, M, N, O, k1, k2, k3, x1, x1', a1, a2, a3) 
    eta = 376.7287537;
    x2 = x1 + a1;
    
    C = zeros(M, N);
    for m = 1:M
        for n = 1:N
            n_p = n_t1(m) + n_t2(n);
            n_m = n_t1(m) - n_t2(n);
            
            C(m,n) = conj(n_z2(n)) / sqrt(abs(n_z1(m))*abs(n_z2(n))*a1*a2) * ...
                ((sin(n_m*x2 - n_t1(m)*x1) + sin(n_t2(n)*x1))/n_m - ....
                (sin(n_p*x2 - n_t1(m)*x1) - sin(n_t2(n)*x1))/n_p);
        end
    end

    Z1 = (eta*k1)./n_z1;
    Z2 = (eta*k2)./n_z2;
    D1 = diag(Z1./abs(Z1));
    D2 = diag(Z2./abs(Z2));
    
    % Compendium page 121
    S11 = inv(conj(D1) + conj(C)*inv(D2)*C.') * (conj(D1) - conj(C)*inv(D2)*C.');
    S12 = 2 * (inv(conj(D1) + conj(C)*inv(D2)*C.') * conj(C));
    S21 = 2 * (inv(C.'*inv(conj(D1))*conj(C) + D2) * C.');
    S22 = inv(C.'*inv(conj(D1))*conj(C) + D2) * (C.'*inv(conj(D1))*conj(C) - D2);
end
