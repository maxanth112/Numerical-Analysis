import sympy as sp

class Spline: 
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d


def cubic_spline_interpolation(x : list, fx : list):
    ''' Max Wiesner
        3/2/21
    Finds the natural cubic spline interpolation coefficients from the given set of 
    points.
    
    --- Required ---
    The lists to be passed need to be the same length, and they both need to have 
    at least 2 elements.

    --- Parameters ---
    x    : list of x coordinates for the spline interpolation
    fx   : list of the outputs from the previous list of x values
    
    --- Returns ---
    Sx_i : the respective free cubic spline approximations '''


    n = len(x)
    if (n < 2 or n != len(fx)):
        raise ValueError('Not enough points were supplied.\n')

    # construct the system of equations 
    h_i = [ (x[i + 1] - x[i]) for i in range(n - 1)]

    a_i = [ (h_i[i]/(h_i[i] + h_i[i + 1])) for i in range(n - 2)] + [0]
    b_i = [2]*n
    c_i = [0] + [ (h_i[i + 1]/(h_i[i] + h_i[i + 1])) for i in range (n - 2)]
    d_i = [0] + [
            6 * ((fx[i + 1] - fx[i]) / h_i[i] - (fx[i] - fx[i - 1]) / h_i[i - 1]) / (h_i[i] + h_i[i - 1]) 
            for i in range(1, n - 1)
        ] + [0]
    
    # solve the system of equations 
    c_j    = c_i + [0]
    c_j[0] = c_i[0] / b_i[0]
    d_j    = [0]*n
    d_j[0] = d_i[0] / b_i[0]

    for i in range(1, n):
        c_j[i] = c_j[i] / (b_i[i] - c_j[i - 1]*a_i[i - 1])
        d_j[i] = (d_i[i] - d_j[i - 1]*a_i[i - 1]) / (b_i[i] - c_j[i - 1]*a_i[i - 1])
    
    Sx = [0]*n
    Sx[-1] = d_j[-1]
    for i in range(n - 2, -1, -1):
        Sx[i] = d_j[i] - c_j[i]*Sx[i + 1]

    # constructing the final coefficients of S(x_i)
    Sx_i = [
        Spline(fx[i],
                ((fx[i+1] - fx[i])/h_i[i] - (Sx[i + 1] + 2*Sx[i])*h_i[i]**2/6), 
                Sx[i]*h_i[i]*h_i[i]/2, 
                (Sx[i+1] - Sx[i])*h_i[i]**2/6)         
        for i in range(n - 1)]

    # print the results 
    index = 0
    for spline in Sx_i: 
        print('Spline Interpolation:')
        print(f'S(x) =>\n a_{index} = {spline.a}\n b_{index} = {spline.b}')
        print(f' c_{index} = {spline.c}\n d_{index} = {spline.d}\n')
        index += 1

    return Sx_i


def spline_error_approximation(func, x_e : int, x_in : list, fx : list):

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
    q = 1
    x = sp.Symbol('x')

    if q == 1:    
        func = x*sp.log(x)
        e    = 8.4
        x    = [8.3, 8.6]
        fx   = [17.56492, 18.50515]

    if q == 2:
        func = sp.sin(sp.exp(x) - 2)
        e    = 0.9
        x    = [0.8, 1.0]
        fx   = [0.22363362, 0.65809197]

    spline = cubic_spline_interpolation(x, fx)

    error = spline_error_approximation(func, e, x, fx)