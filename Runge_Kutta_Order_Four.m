function [t, w] = Runge_Kutta_Order_Four(df, a, b, h, alpha)
    
    n = (b - a)/h;
    t = zeros(n, 1);
    w = zeros(n, 1);
    t(1) = a;
    w(1) = alpha;
    K = zeros(4, 1);
    
    for i = 1:n
        K(1) = h*df(t(i), w(i));
        K(2) = h*df(t(i) + (h/2), w(i) + (K(1)/2));
        K(3) = h*df(t(i) + (h/2), w(i) + (K(2)/2));
        K(4) = h*df(t(i) + h, w(i) + K(3));
        
        w(i + 1) = w(i) + (K(1) + 2*K(2) + 2*K(3) + K(4))/6;
        t(i + 1) = t(i) + h;
    end
    
    subplot
    fimplicit(@(x, y) 2 + x - log(y/(2*x)) - y)
    hold on 
    plot(t, w, '--r')
    title('Runge-Kutta Order Four with h=0.0025')
    xlabel('x')
    ylabel('y')
    legend('show', 'location', 'best')
end


