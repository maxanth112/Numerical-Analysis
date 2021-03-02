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

    return Sx_i


# funciton calls
if __name__ == "__main__":
    
    # setting the question 
    q = 2
    
    if q == 1:    
        x = [8.3, 8.6]
        y = [17.56492, 18.50515]

    if q == 2:
        x = [0.8, 1.0]
        y = [0.22363362, 0.65809197]

    spline = cubic_spline_interpolation(x, y)
    for i in range(len(spline)): 
        print(f'\nS(x) =>\n a_{i} = {spline[i].a}\n b_{i} = {spline[i].b}')
        print(f' c_{i} = {spline[i].c}\n d_{i} = {spline[i].d}')