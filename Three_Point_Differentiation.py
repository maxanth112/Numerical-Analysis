import sympy as sp
from Lagrange_Interpolating_Polynomials import find_max_save_value

def three_point_differentiation(p : list) -> list:
    ''' Max Wiesner
        3/10/21
    Identifyes the most accurate three-point numerical derivative formula to determine the 
    missing entries of the 2D array of points p that is passed, uses either the endpoint 
    formula or the midpoint formula.
    
    --- Required ---
    the x coordinate in the dictionary of p must be evenly spaced apart, that is for any
    i = 0, ..., n it is that x_(i + 1) - x_i = h for the constant integer h. There must also
    be at least 2 points in the list of p such that we are able to determine h
    
    --- Parameters ---
    p : (points) a list of p [x, fx], in the form of [[x_0, fx_0], ..., [x_n, fx_n]]

    --- Returns ---
    p : (updated) list of list of estimations for the derivative of the given p in the third
        index, calculated by the the optimal method '''


    sorted(p, key = lambda x : x[0])
    n = len(p)
    h = p[1][0] - p[0][0]
    h_dir = -1 if h < 0 else 1   
    p = [i + [0]*5 for i in p]

    print(f'\n[Given]:\nh = {p[1][0]} - {p[0][0]}')
    for i in range(n):
        print(f'f({p[i][0]}) = {p[i][1]}')

    print(f'\n[Three Point Derivative Calcs]:')
    for i in range(n):

        if i == 0 or i == (n - 1): # endpoint
            t_dir = -1 if (h_dir > 0 and i == (n - 1)) else h_dir
            t = -h if t_dir == -1 else h            
            dfx = (1/(2*t))*(-3*p[i][1] + 4*p[i + t_dir][1] - p[i + 2*t_dir][1])
            print(f' f\'({p[i][0]}) : using endpoint => {dfx:1.5f}')
            p[i][3] = [p[i][0], p[i + t_dir][0], p[i + 2*t_dir][0]]

        else: # midpoint 
            dfx = (1/(2*h))*(p[i + h_dir][1] - p[i - h_dir][1])
            print(f' f\'({p[i][0]}) : using midpoint => {dfx:1.5f}')
            p[i][3] = [p[i][0], p[i + h_dir][0], p[i - h_dir][0]]
            print(p[i][3])
        p[i][2] = dfx
       
    return p
 
 

def three_point_error(p : list, func) -> list:
    ''' 
    Calculated the actual errors of the list of points and their approximated derivatives that
    were computed and passed from the above function, also computes the error bound, and prints
    results.
    
    --- Required ---
    list of points must have pre estimated derivatives, and a list of inputs in the current domain
    to try and maximize, or passed through the above function successfully
    
    --- Parameters ---
    p : (points) a list of p [x, fx, [m1, m2, ..], e_dfx], where e_dfx is the estimated derivative, 
        and [m1, ...] is the list of points in the domain to try for the maximizing equation

    --- Returns ---
    p : (updated) list of the original p plus the error bound and actual error for each point '''


    x = sp.Symbol('x')
    n = len(p)
    h = p[1][0] - p[0][0]

    df = sp.diff(func, x)
    df2 = sp.diff(df, x)
    df3 = sp.diff(df2, x)
  
    dfx = sp.lambdify(x, df)
    print(f'\n[Error Calcs]:')
    for i in range(n):

        d = 3 if ( i == 0 or i == (n - 1)) else 6
        print(f' for f\'({p[i][0]}) = {p[i][2]:<1.5f}')     
        [xm, fxm] = find_max_save_value(df3, p[i][3])
        p[i][4] = sp.Abs((h**2)/d)*fxm # store bound
        p[i][5] = sp.Abs(p[i][2] - dfx(p[i][0])) # store actual 

        print(f'  bound  : Pick x = {xm} to maximize |(h^2/3)*{df3}| = {p[i][4]:1.5f}')
        print(f'  actual : |{p[i][1]:<1.5f} - {df}| =  |{p[i][1]:<1.5f} - {dfx(p[i][0]):<1.5f}| = {p[i][5]:<1.5f}')

    # print the final table
    for i in range(n):      
        if i == 0:
            print('\n[Final Table]:')
            print(f' {"-"*52}')
            print(f'|{" x":<8}|{" f(x)":<8}|{" df(x)":<10}| { " e actual":<10}| {" e bound":<10}|')
            print(f' {"-"*52}')
        print(f'|{p[i][0]:<8}|{p[i][1]:<8.5f}| {p[i][2]:<9.5f}|  {p[i][5]:<9.5f}|  {p[i][4]:<9.5f}|')
    print(f' {"-"*52}\n')

    return p




# funciton call examples 
if __name__ == "__main__":
    
    # setting the question 
    x = sp.Symbol('x')
    q = 2

    if q == 1: 
        p = [ [1.1, 9.025013], [1.2, 11.02318], [1.3, 13.46374], [1.4, 16.44465] ]
        fx = sp.exp(2*x)
    if q == 2: 
        p = [ [-3/10, -0.27652], [-2/10, -0.25074], [-1/10, -0.16134], [0, 0] ]
        fx = sp.exp(2*x) - sp.cos(2*x)

    p = three_point_differentiation(p)
    p = three_point_error(p, fx)
    


        

