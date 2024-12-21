function Lj = lagrangepoly(x, xk, lambda, j)
% computes Lagrange polynomial Lj(x) at x given j, node distribution xk and
% lambda = L'(xk), where L(x) = prod(x-xk) (Barycentric interpolation)
    if all(abs(xk-x) > 1e-5)
        den = sum(lambda ./ (x - xk), 1);
        num = lambda(j+1) ./ (x - xk(j+1));
        Lj = num ./ den;
    elseif abs(xk(j+1)-x) < 1e-5
        Lj = 1;
    else 
        Lj = 0;
    end
end