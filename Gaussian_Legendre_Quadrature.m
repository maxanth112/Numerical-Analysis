function [weights, points, apprx_ans] = Gaussian_Legendre_Quadrature(func, n, a, b)

    % used matlab code as reference, this substitutes table 4.12
    % simplify n instances for recurrence formula
    nm1 = n - 1;
    n1 = nm1 + 1; 
    n2 = n1 + 1;
    % generate vector n points, spaced with (1 - (-1))/(n - 1)
    xu = linspace(-1, 1, n1)';
    % starting guess
    Pn = cos((2*(0:nm1)'+1)*pi/(2*nm1+2))+(0.27/n1)*sin(pi*xu*nm1/n2);
    % generate a n1xn2 Jacobian matrix of zeros and a cooresponding derivative matrix
    Jf = zeros(n1,n2);
    dJf = zeros(n1,n2);
    tmp = 2;

    while max(abs(Pn - tmp)) > 10^(-5)

        Jf(:, 1) = 1;
        Jf(:, 2) = Pn;
        for i = 2:n1
            Jf(:, i+1) = ( (2*i - 1)*Pn.*Jf(:, i) - (i - 1)*Jf(:, i - 1) )/ i;
        end

        dJf = (n2)*( Jf(:, n1) - Pn.*Jf(:, n2) )./(1 - Pn.^2);      
        tmp = Pn;
        Pn = tmp - Jf(:, n2)./ dJf; 
    end

    % points of the quadrature 
    points = (a*(1 - Pn) + b*(1 + Pn)) / 2;      
    % weights of the gaussian quadrature 
    weights = (b - a)./( (1 - Pn.^2).*dJf.^2)*(n2/n1)^2;

    apprx_ans = 0;
    for i = 1:n
        apprx_ans = apprx_ans + func(points(i))*weights(i);
    end
end



