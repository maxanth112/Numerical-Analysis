function [t, w] = Eulers_Method(df, a, b, h, alpha)
    
    n = (b - a)/h;
    t = zeros(n, 1);
    w = zeros(n, 1);
    
    t(1) = a;
    w(1) = alpha;
    
    for i = 1:n
        t(i + 1) = t(i) + h;
        w(i + 1) = w(i) + h*df(t(i), w(i));
    end

    subplot
    fimplicit(@(x, y) 2 + x - log(y/(2*x)) - y)
    hold on 
    plot(t, w, '--r')
    title('Eulers Method with h=0.01')
    xlabel('x')
    ylabel('y')
    legend('show', 'location', 'best')
end


