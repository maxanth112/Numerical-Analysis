function p = Newton_Raphson_Method(f, p_0, N, e)
    % Max Wiesner
    % 2/15/21
    % 
    % This function iterates using the Newton Raphson method to 
    % to find a root on the given continuous funciton, can be though 
    % of as successive tangent approximation.
    % 
    % --- Requirements --- 
    % must declare x prior, ie. 'syms x' then call the function, p_0
    % cannot be such that f'(p_0) = 0.
    %
    % --- Parameters ---
    % f   : the function to find the root on
    % p_0 : initial starting value
    % N   : maximum amount of iterations 
    % e   : tolerance, maximum amount of error willing to accept 
    %
    % --- Returns ---
    % p   : approximated root of f, to within an error of e
    
    syms x;
    df = diff(f,x);
    df = matlabFunction(df);
    f = matlabFunction(f);        
    if f(p_0) == 0
        fprintf("The supplied starting value is a root. ");
        p = p_0;
        return;
    end
    
    i = 1;
    while i <= N
        p = p_0 - f(p_0)/df(p_0);
        fprintf("\np_%d = %-5.10f", i, p);        
        if abs(p - p_0) < e
            fprintf("\n\np found after %d iterations.", i);
            return;
        end
        
        i = i + 1;
        p_0 = p;
    end
    fprintf("\nThe method failed after %d iterations. ", i);
    return;
end

% example 
% >> sysm x
% >> f = x^3 - 0.165*x + 4*(10^-4)
% >> ans = Newton_Raphson_Method(f,0,10,10^(-5))

% output 
% p_1 = 0.0024242424
% p_2 = 0.0024243288
%
% p found after 2 iterations.
% ans =
%
%   0.0024

