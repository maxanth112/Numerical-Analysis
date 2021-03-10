import sympy as sp
from Lagrange_Interpolating_Polynomials import find_max_save_value


def numerical_differentiation(p : list, func, m = 'both'):
    ''' Max Wiesner
        3/9/21
    Computes the the forward and backward numerical differentiations of the given points when
    possible, and the error bound, and actual error. Results are formatted and printed in a 
    table.
    
    --- Required ---
    the x coordinate in the dictionary of p must be evenly spaced apart, that is for any
    i = 0, ..., n it is that x_(i + 1) - x_i = h for the constant integer h. There must also
    be at least 2 points in the list of p such that we are able to determine h
    
    --- Parameters ---
    p : (points) a list of p [x, fx], in the form of [[x_0, fx_0], ..., [x_n, fx_n]]
    m : (method) the divided difference m to be performed; 'forward', 'backward', 'both'

    --- Returns ---
    p : (updated points) list of list of estimations for the derivative of the given p, 
        calculated by the choosen m, adding index 2 for forward estimation, and index 3
        for backward estimation, returned as displayed by the functions last print statement '''


    sorted(p, key = lambda x : x[0])
    x = sp.Symbol('x')
    df = sp.diff(func, x)
    dfx = sp.lambdify(x, df)
    df2 = sp.diff(df, x)
    n = len(p)
    h = p[1][0] - p[0][0]
    p = [i + [None]*8 for i in p]

    print(f'\n[Given]:\nh = {p[1][0]} - {p[0][0]}')
    for i in range(n):
        print(f'f({p[i][0]}) = {p[i][1]}')
    h_dir = -1 if h < 0 else 1   
    opt = [{ 'method' : 'Forward', 'index' : 2, 'sign' : 1, 'bound' : n },
           { 'method' : 'Backward', 'index' : 5, 'sign' : -1, 'bound' : -1}]

    print(f'\n[Derivative Calcs]:')
    for i in range(len(opt)):
        m = opt[i]
        print(f'{m["method"]}')
        for j in range(n):
            if (j + m['sign']*h_dir) == m['bound']:
                print(f' f\'({p[j][0]}) = can\'t be calculated since f(x_i {"-" if m["sign"] < 1 else "+"} h) wasn\'t provided.')
                p[j][m['index'] + 1] = None
                p[j][m['index'] + 2] = None
            else: 
                fxh = p[j + h_dir*m['sign']][1]
                dfi = m['sign']*(fxh - p[j][1])/h
                p[j][m['index']] = dfi
                print(f' f\'({p[j][0]}) = {dfi:<1.5f}')


    print(f'\n[Error Calcs]:')
    for i in range(len(opt)):
        m = opt[i]
        print(f'{m["method"]}')

        for j in range(n):
            print(f' for f\'({p[j][0]}) = {p[j][1]:<1.5f}')
            if (j + m['sign']*h_dir) == m['bound']:
                print(f'  error of f\'({p[j][0]}) = can\'t be calculated since f(x_i {"-" if m["sign"] < 1 else "+"} h) wasn\'t provided.')
                p[j][m['index'] + 1] = None
                p[j][m['index'] + 2] = None
            else: 
                [xm, fxm] = find_max_save_value(df2, [p[j][0], p[j + m['sign']][0]])
                p[j][m['index'] + 1] = sp.Abs(h/2)*fxm # store bound
                p[j][m['index'] + 2] = sp.Abs(p[j][m['index']] - dfx(p[j][0])) # store actual 
                print(f'  bound  : Pick x = {xm} to maximize |h/{2}*{df2}| = {p[j][m["index"] + 1]:1.5f}')
                print(f'  actual : |{p[j][1]} - {df}| =  |{p[j][m["index"]]:<1.5f} - {dfx(p[j][0]):<1.5f}| = {p[j][m["index"] + 2]:<1.5f}')

    # print the final table
    for i in range(n):      
        if i == 0:
            print('\n[Final Table]:   [B] = Backward, [F] = Forward')
            print(f' {"-"*107}')
            print(f'|{" x":<8}|{" f(x)":<8}|{" [F] df(x)":<14}| { "[F] e bound":<13}| {"[F] e actual":<13}| {"[B] df(x)":<13}| {"[B] e bound":<13}| {"[B] e actual":<13}|')
            print(f' {"-"*107}')

        na = f'na{" ".ljust(11)}'
        cond_f = p[i][2] != None
        cond_b = p[i][5] != None
        df_f = f'{p[i][2]:<13.5f}' if cond_f else na
        ef_b = f'{p[i][3]:<13.5f}' if cond_f else na
        ef_a = f'{p[i][4]:<13.5f}' if cond_f else na
        df_b = f'{p[i][5]:<13.5f}' if cond_b else na
        eb_b = f'{p[i][6]:<13.5f}' if cond_b else na
        eb_a = f'{p[i][7]:<13.5f}' if cond_b else na
        print(f'|{p[i][0]:<8}|{p[i][1]:<8}| {df_f}| {ef_b}| {ef_a}| {df_b}| {eb_b}| {eb_a}|')
    print(f' {"-"*107}\n')

    return p




# funciton call examples 
if __name__ == "__main__":
    
    # setting the question 
    x = sp.Symbol('x')
    q = 1

    if q == 1: 
        p = [ [1/2, 0.4794], [6/10, 0.5646], [7/10, 0.6442] ]
        fx = sp.sin(x)
    if q == 2: 
        p = [ [-3/10, 1.9507], [-2/10, 2.0421], [-1/10, 2.0601] ]

    numerical_differentiation(p, fx)


        

