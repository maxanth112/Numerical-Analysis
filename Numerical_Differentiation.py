import sympy as sp

def numerical_differentiation(p : list, method):
    ''' Max Wiesner
        3/8/21
    Func Description
    
    --- Required ---
    the x coordinate in the dictionary of p must be evenly spaced apart, that is for any
    i = 0, ..., n it is that x_(i + 1) - x_i = h for the constant integer h. There myst also
    be at least 2 p in the list of p such that we are able to determine h
    
    --- Parameters ---
    p : a list of p [x, fx], in the form of [[x_0, fx_0], ..., [x_n, fx_n]]
    method : the divided difference method to be performed; 'f' => forward or 'b' => backward

    --- Returns ---
    dfx    : list of list of estimations for the derivative of the given p, calculated by the
             choosen methods '''

    sorted(p, key = lambda x : x[0])
    n = len(p)
    h = p[1][0] - p[0][0]
    dfx = [[0, 0] for i in range(n)]

    print(f'[Given]:\nh = {p[1][0]} - {p[0][0]}')
    for i in range(n):
        print(f'f({p[i][0]}) = {p[i][1]}')

    print(f'\n[Calculations]:')
    for i in range(n): 

        fx_i = p[i][1]
        h_dir = -1 if h < 0 else 1         
        if method == 'forward' or method == 'both': 

            if (i + h_dir) == n:
                print(f'[F] f\'({p[i][0]}) = can\'t be calculated using forward method since f(x_i + h) wasn\'t provided.')
            else:
                fx_i_ph = p[i + h_dir][1]
                dfx_i = ( fx_i_ph - fx_i ) / h
                dfx[i][0] = dfx_i            
                print(f'[F] f\'({p[i][0]}) = {dfx_i:<1.5f}')

        if method == 'backward' or method == 'both':

            if (i - h_dir) == -1:
                print(f'[B] f\'({p[i][0]}) = can\'t be calculated using backward method since f(x_i - h) wasn\'t provided.')
            else:
                fx_i_mh = p[i - h_dir][1]
                dfx_i = ( fx_i - fx_i_mh ) / h
                dfx[i][1] = dfx_i
                print(f'[B] f\'({p[i][0]}) = {dfx_i:<1.5f}')

    # print the final table
    for i in range(n):      
        if i == 0:
            print('\n[Final Table]:')
            print(f' {"-"*47}')
            print(f'|{" x":<8}|{" f(x)":<8}|{" [F] df(x)":<14}|{" [B] df(x)":<14}|')
            print(f' {"-"*47}')
        forw = f'{dfx[i][0]:<13.5f}' if dfx[i][0] != 0 else f'na{" ".ljust(11)}'
        back = f'{dfx[i][1]:<13.5f}' if dfx[i][1] != 0 else f'na{" ".ljust(11)}'
        print(f'|{p[i][0]:<8}|{p[i][1]:<8}| {forw}| {back}|')
    print(f' {"-"*47}')

    return dfx


#
def numerical_differentiation_error(p : list, method):
    '''
    Func Description
    
    --- Required ---
    the x coordinate in the dictionary of p must be evenly spaced apart, that is for any
    i = 0, ..., n it is that x_(i + 1) - x_i = h for the constant integer h. There myst also
    be at least 2 p in the list of p such that we are able to determine h
    
    --- Parameters ---
    p : a list of p [x, fx], in the form of [[x_0, fx_0], ..., [x_n, fx_n]]
    method : the divided difference method to be performed; 'f' => forward or 'b' => backward

    --- Returns ---
    dfx    : list of list of estimations for the derivative of the given p, calculated by the
             choosen methods '''




# funciton call examples 
if __name__ == "__main__":
    
    # setting the question 
    x = sp.Symbol('x')
    q = 2

    if q == 1: 
        p = [
            [1/2, 0.4794],
            [6/10, 0.5646],
            [7/10, 0.6442] 
        ]

    if q == 2: 
        p = [
            [-3/10, 1.9507],
            [-2/10, 2.0421],
            [-1/10, 2.0601] 
        ]

    numerical_differentiation(p, 'both')


        

