function [a, b, c, d] = Cubic_Spline_Interpolation(x, fx)
    % Max Wiesner
    % 3/2/21
    
    % Constructs a natural cubic spline interpolation from the given arrays
    % of inputs.
    
    % --- Requirements --- 
    % must declare x prior, ie. 'syms x' then call the function, must
    % supply two vectors of same length with a length greater than 1
   
    % --- Parameters ---
    % x  : vector of inputs to the function we will interpolate
    % fx : vector of outputs 
   
    % --- Returns ---
    % [a, b, c, d] : coefficients of the interpolated cubic spline
    
    n = length(x);
    for i = 1 : n - 1
        h(i) = x(i + 1) - x(i);
        fprintf('%f\n', h(i));
    end 
    
    A(1, 1) = 1;
    A(n, n) = 1;
    f = zeros(n, 1);
    for i = 2 : n - 1
        A(i, i) = 2*(h(i) + h(i - 1));
        f(i) = 6*( (fx(i + 1) - fx(i))/h(i) - (fx(i) - fx(i - 1))/h(i - 1));
    end
    
    for i = 2 : n - 2
        A(i, i + 1) = h(i + 1);
    end
    
    s = A\f;
    
    for i = 1 : n - 1
        a(i) = (s(i + 1))/(6*h(i));
        b(i) = s(i)/2;
        c(i) = (fx(i + 1) - fx(i))/h(i) - (2*h(i)*s(8) + h(i)*s(i + 1))/6;
        d(i) = fx(i);
    end
    
    return;
end 

