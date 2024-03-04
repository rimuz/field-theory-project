function n_z = nzfromnt(n_t, k)
    N = length(n_t);
    n_z = zeros(N, 1);

    for n = 1:N
        if k > n_t(n)
            n_z(n) = sqrt(k^2 - n_t(n)^2);
        else
            n_z(n) = -1i*sqrt(n_t(n)^2 - k^2);
        end
    end
end