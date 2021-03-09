import sympy as sp

def divided_difference(points : list, method):
    ''' Max Wiesner
        X/X/21
    Func Description
    
    --- Required ---
    the x coordinate in the dictionary of points must be evenly spaced apart, that is for any
    i = 0, ..., n it is that x_(i + 1) - x_i = h for the constant integer h. There myst also
    be at least 2 points in the list of points such that we are able to determine h
    
    --- Parameters ---
    points : a list of points [x, fx], in the form of [[x_0, fx_0], ..., [x_n, fx_n]]
    method : the divided difference method to be performed; 'f' => forward or 'b' => backward

    --- Returns ---
    dfx    : list of estimations for the derivative of the given points, calculated by the
             choosen method '''

    sorted(points, key = lambda x : x[0])
    n = len(points)
    h = points[1][0] - points[0][0]
    dfx = [0]*n
    msg = "Forward" if method == 'f' else 'Backward'

    print(f'\n[{msg} Difference]:')
    for i in range(n): 

        fx_i = points[i][1]
        h_dir = -1 if h < 0 else 1 
        
        if method == 'f':
            
            if (i + h_dir) == n:
                print(f'f\'({fx_i}) = can\'t be calculated since f(x_i + h) wasn\'t provided.')
            else:
                fx_i_ph = points[i + h_dir][1]
                dfx_i = ( fx_i_ph - fx_i ) / h
                dfx[i] = dfx_i

                print(f'f\'({fx_i}) = {dfx_i}')

        else: 

            
    

   

    return 




# funciton call examples 
if __name__ == "__main__":
    
    # setting the question 
    x = sp.Symbol('x')
    q = 1

    if q == 1: 

        points = [
            [-3/10, 1.9507],
            [-2/10, 2.0421],
            [-1/10, 2.0601] 
        ]

    divided_difference(points, 'f')


        

