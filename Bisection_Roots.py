import math 

def bisection(f,a,b,N,e):
    ''' Max Wiesner 
        2/15/21
    
    Approximate solution to the root finding problem 
    for a given function f, of the form f(x)=0, on the interval [a,b].
    
    --- Criteria ---
    The function f must be continuous on the interval [a,b], with f(a)
    and f(b) of opposite sign, by the Mean Value Theorem, we thus know 
    there exists at least one root.
       
    --- Parameters ---
    f : function 
    a : lower bound 
    b : upper bound 
    N : maximum number of iterations before force quitting 
    e : maximum error allowed (tolerance)
    
    --- Returns ---
    Approximated root p such that f(p)<e, or real root such that f(p)=0. '''
    
    fa = f(a)
    fb = f(b)
    if fa*fb >= 0:
        print(f'fa = {fa}\nfb = {fb}\n')
        print("Criteria not met for bisection method.")
        return None
    if fa == 0:
        return a
    if fb == 0:
        return b
        
    a_n = a
    b_n = b
    i = 1
    while i < N:
        p_n = (a_n + (b_n - a_n)/2)
        fp = f(p_n)
        print(f'p_{str(i).ljust(2)} = {p_n:8.14f}')
        
        if fp == 0 or (b_n - a_n)/2 < e:
            print(f'\np found after {i} iterations.')
            return p_n
        
        i += 1
        if fa*fp < 0:
            b_n = p_n
        else:
            a_n = p_n
    print(f"\nMethod failed after {N} iterations.")
    return None

# calling with lambda functions 
fx = lambda x : math.tan(math.pi*x) - 6
approx_root = bisection(fx,0,0.5,11,10e-10)