function [w, t] = Backward_Implicit_Euler_Method(f, df, a, b, alpha, h)
    
    n = (b - a)/h;
    t = zeros(n, 1);
    w = zeros(n, 1);
    t(1) = a;
    w(1) = alpha;

	wi = alpha;
    x = wi;
    tn = a;
    
    for i = 1 : n
        tn = tn + h;
        for j = 1:25 
            % newtons method 
            x = x - (h * df(tn,x) - 1) \ (h * f(tn,x) - x + wi);
        end
        wi = x;
        w(i + 1) = x;
        t(i+1) = tn;
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