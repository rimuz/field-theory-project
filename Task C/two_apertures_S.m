% Calculate S parameters for two waveguides with height a1 and a3 feeding 
% waveguide with height a2, a1 < a2 and a3 < a2. x1 and y1 are the offset 
% of the apertures, n_t1, n_t2 and n_t3 the eigenvalues at the three
% waveguides, M, N and P are the number of modes.
function [S11,S12,S21,S22] = two_apertures_S(n_t1, n_t2, n_t3, n_z1, n_z2, n_z3, M, N, P, k1, k2, k3, x1, y1, a1, a2, a3) 
    eta = 376.7287537;
    x2 = x1 + a1;
    y2 = y1 + a2;

    C1 = zeros(M, N);
    for m = 1:M
        for n = 1:N
            n_p = n_t1(m) + n_t2(n);
            n_m = n_t1(m) - n_t2(n);
            
            C1(m,n) = conj(n_z2(n)) / sqrt(abs(n_z1(m))*abs(n_z2(n))*a1*a2) * ...
                ((sin(n_m*x2 - n_t1(m)*x1) + sin(n_t2(n)*x1))/n_m - ...
                (sin(n_p*x2 - n_t1(m)*x1) - sin(n_t2(n)*x1))/n_p);
        end
    end
    
    C2 = zeros(P, N);
    for p = 1:P
        for n = 1:N
            n_p = n_t3(p) + n_t2(n);
            n_m = n_t3(p) - n_t2(n);
            
            C2(p,n) = conj(n_z2(n)) / sqrt(abs(n_z3(p))*abs(n_z2(n))*a3*a2) * ...
                ((sin(n_m*y2 - n_t3(p)*y1) + sin(n_t2(n)*y1))/n_m - ...
                (sin(n_p*y2 - n_t3(p)*y1) - sin(n_t2(n)*y1))/n_p);
        end
    end

    C = [C1; C2];

    Z1 = (eta*k1)./n_z1;
    Z2 = (eta*k2)./n_z2;
    Z3 = (eta*k3)./n_z3;
    D1 = [diag(Z1./abs(Z1)) zeros(M,P); zeros(P,M) diag(Z3./abs(Z3))];
    D2 = diag(Z2./abs(Z2)); 
    

    % Compendium page 121
    S11 = inv(conj(D1) + conj(C)*inv(D2)*C.') * (conj(D1) - conj(C)*inv(D2)*C.');
    S12 = 2 * (inv(conj(D1) + conj(C)*inv(D2)*C.') * conj(C));
    S21 = 2 * (inv(C.'*inv(conj(D1))*conj(C) + D2) * C.');
    S22 = inv(C.'*inv(conj(D1))*conj(C) + D2) * (C.'*inv(conj(D1))*conj(C) - D2);
end
