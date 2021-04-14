function [w, t] = Backward_Implicit_Euler_Method(df, a, b, h, alpha)
    
    n = (b - a)/h;
    t = zeros(n, 1);
    w = zeros(n, 1);
    
    t(1) = a;
    w(1) = alpha;

    for i = 1 : n - 1
        t(i + 1) = t(i) + h;
        x = w(i);
        % NR starts 
        for j = 1 : 20
            syms t;
            syms x;         
            dff = matlabFunction(df);
            fi = x - dff(t(i + 1), x)*h - w(i);
            
            dfi = matlabFunction(diff(fi, x));
            
            x = x - fi/dfi(x);           
        end
        w(i + 1) = x;
    end
    
    
    subplot
    fimplicit(@(x, y) 2 + x - log(y/(2*x)) - y)
    hold on 
    plot(t, w, '-r')
    title('Implicit Euler with h=0.005')
    xlabel('x')
    ylabel('y')
    legend('show', 'location', 'best')
end