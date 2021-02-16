import math

def fixed_point(g,p_0,N,e):
    ''' Max Wiesner 
        2/15/21
    
    Finds the number p for the function g such that g(p)=p. 
    
    --- Criteria ---
    The function g must be C[a,b] and g(x) in [a,b] for all x in [a,b]. 
    This guarantees at least one fixed point in the interval, if in 
    addition, g'(x) exists on (a,b) and a positive k < 1 exists with 
    |g'(x)|<=k for all x in (a,b), then there exists exactly one fixed 
    point in [a,b].
    
    --- Parameters ---
    g   : function 
    p_0 : starting p value 
    N   : maximum amount of iterations allowed 
    e   : tolerance, maximum error allowed 
    
    --- Returns ---
    fixed point p within the given tolerance. '''
    
    i = 1
    while i <= N:
        p_n = g(p_0)
        
        print(f'p_{str(i).ljust(2)} = {p_n:8.14f}')
        if abs(p_n - p_0) < e:            
            print(f'\np found after {i} iterations.')
            return p_n
        
        i += 1
        p_0 = p_n
    print(f"/nThe method failed after {N} iterations. ")
    return None

# calling with lambda functions 
g = lambda x : 1/math.sqrt(1+x)
fixed_point(g,2,100,10e-5)