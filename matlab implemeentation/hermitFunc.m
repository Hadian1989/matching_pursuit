function y = hermitFunc(n, x)
    % Compute the normalization factor
    norm_factor = sqrt(factorial(n) * 2^n * sqrt(pi));
    
    % Define the normalized Hermite function
    y = exp(-x.^2 / 2) .* hermiteH(n, x) / norm_factor;
end