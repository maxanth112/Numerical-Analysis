import sympy as sp

class Spline: 
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d


def cubic_spline_interpolation(x : list, fx : list, dfs = 0):
    ''' Max Wiesner
        3/2/21
    Finds the natural/clamped cubic spline interpolation coefficients from the 
    given set of points.
    
    --- Required ---
    The lists to be passed need to be the same length, and they both need to have 
    at least 2 elements.

    --- Parameters ---
    x    : list of x coordinates for the spline interpolation
    fx   : list of the outputs from the previous list of x values
    
    --- Returns ---
    Sx_i : the respective cubic spline approximations '''


    n = len(x)
    if (n < 2 or n != len(fx)):
        raise ValueError('Not enough points were supplied.\n')
    # initialize the working arrays 
    l, u, z = [0]*(n + 1), [0]*(n + 1), [0]*(n + 1)
    b, c, d = [0]*(n + 1), [0]*(n + 1), [0]*(n + 1)
    a = fx
    a += [0] + [0]
    
    # construct the system of equations 
    h = [ (x[i + 1] - x[i]) for i in range(n - 1)] + [0]
    alpha = [0] + [ (3/h[i])*(a[i + 1] - a[i]) - (3/h[i - 1])*(a[i] - a[i - 1]) for i in range(1, n - 1) ] + [0, 0]
  
    if (dfs): # clamped spline ops 1
        alpha[0] = 3*(a[1] - a[0])/h[0] - 3*dfs[0]
        alpha[n - 1] = 3*dfs[1] - 3*(a[n - 1] - a[n - 2])/h[n - 2]
        l[0] = 2*h[0]
        u[0] = 0.5
        z[0] = alpha[0]/l[0]

    # solve for the system of equations 
    for i in range(1, n - 1):
        l[i] = 2*(x[i + 1] - x[i - 1]) - h[i - 1]*u[i - 1]
        u[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i - 1]*z[i - 1])/l[i]
        
    if (dfs): # clamped spline ops 2
        l[n - 1] = h[n - 2]*(2 - u[n - 2])
        z[n - 1] = (alpha[n - 1] - h[n - 2]*z[n - 2])/l[n - 1]
        c[n - 1] = z[n - 1]

    for i in range(n - 2, -1, -1):
        c[i] = z[i] - u[i]*c[i + 1]
        b[i] = (a[i + 1] - a[i])/h[i] - h[i]*(c[i + 1] + 2*c[i])/3
        d[i] = (c[i + 1] - c[i])/(3*h[i])
    
    Sx = [ Spline(a[i], b[i], c[i], d[i]) for i in range(n - 1) ]

    # print the results before returning 
    for j in range(len(Sx)): 
        print('Spline Interpolation:')
        print(f'S(x) =>\n a_{j} = {Sx[j].a}\n b_{j} = {Sx[j].b}')
        print(f' c_{j} = {Sx[j].c}\n d_{j} = {Sx[j].d}\n')

    return Sx


def spline_error_approximation(func, x_e : int, x_in : list, fx : list):
    '''
    Finds the actual error from the above cubic spline interpolation method, and the error
    of its derivative
    
    --- Required ---
    Same requirements as the above function

    --- Parameters ---
    func : the actual function that we estimated 
    x_e  : the value that we will use to estimate the error and compare outputs
    x_in : list of input values for the cubic spline interpolation 
    fx   : list of the outputs from the previous list of x values
    
    --- Returns ---
    void : doesn't return anything, for printing out an analysis of the interpolation '''


    x = sp.Symbol('x')
    Sxi    = cubic_spline_interpolation(x_in, fx)[0] # just fo the first one 
    f_Sxi  = Sxi.a + Sxi.b*(x - x_in[0]) + Sxi.c*(x - x_in[0])**2 + Sxi.d*(x - x_in[0])**3
    
    header = [f'({x_e})', f"'({x_e})"]
    print('Error Evaluation: ')
    for i in range(len(header)):

        ff = sp.lambdify(x, func)
        ss = sp.lambdify(x, f_Sxi)
        if i == 1:
            ff = sp.lambdify(x, sp.diff(f_Sxi, x))
            ss = sp.lambdify(x, sp.diff(func, x))
        print(f'f{header[i]} = {ff(x_e)}')
        print(f'S{header[i]} = {ss(x_e)}')
        print(f'f_err  = |S{header[i]} - f{header[i]}|\n       = {abs(ff(x_e) - ss(x_e))}\n')

    return 



# funciton calls
if __name__ == "__main__":
    
    # setting the question 
    q = 2
    x = sp.Symbol('x')

    if q == 1:    
        func = x*sp.log(x)
        e    = 8.4
        x    = [8.3, 8.6]
        fx   = [17.56492, 18.50515]
        df   = [1.116256, 1.151872]

    if q == 2:
        func = sp.sin(sp.exp(x) - 2)
        e    = 0.9
        x    = [0.8, 1.0]
        fx   = [0.22363362, 0.65809197]
        df   = [2.1691753, 2.0466965]

    # spline = cubic_spline_interpolation(x, fx)
    
    spline = cubic_spline_interpolation(x, fx, df)
    # error = spline_error_approximation(func, e, x, fx)