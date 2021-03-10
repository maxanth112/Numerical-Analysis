import sympy as sp

class Point: 
    def __init__(self, x, fx):
        self.x = x
        self.fx = fx
        
def lagrange_interpolation(p : list, x_p : float, n : int, p_out : bool = True) -> float:
    ''' Max Wiesner 
        3/1/21
    
    Approximation by Lagrange interpolation by constructing Lagrange Polynomials of degree n,
    given a set of known points.
    
    --- Criteria ---
    The degree n of approximation requires n + 1 unique points to by given. 
       
    --- Parameters ---
    p     : list of known unique points  
    x_p   : the point in question to be approximated 
    n     : the degree of the Lagrange polynomial to be constructed for the approximation
    p_out : disable of enable a printout of the computed polynomial P_x, default is true
    
    --- Returns ---
    x_r : result, the approximated value using Lagrange interpolating polynomial of degree n '''
    
    if (len(p) < n or len(p) < 2): 
        raise ValueError('Not enough known points provided.')
    if (n < 1):
        raise ValueError('Degree has to be at least 1.')

    p = sorted(p, key = lambda point : sp.Abs(x_p - point.x))
    Px = 0.0    
    for i in range(n + 1):        
        L_i = 1
        for j in range(n + 1):
            if i != j:
                L_i *= (x - p[j].x) / (p[i].x - p[j].x) 
        Px += p[i].fx*L_i

    # the final polynomial and answer
    P_x = sp.lambdify(x, Px)
    x_r = P_x(x_p)

    # printing out the constructed Lagrange polynomial and problem specifics
    if p_out:
        print(f'\n[Lagrange Polynomial]:\nP_{n}(x) = {sp.expand(Px)}\nP_{n}({x_p}) = {P_x(x_p)}\n')
        
    return x_r



def lagrange_error_bound(fx, p : list, x_p : float, n : int, p_out : bool = True) -> float: 
    '''  
    Finds a bound for the given function over the interval used in the Lagrange polynomial 
    interpolation, then compares the bound to the actual error using the error formula
    
    --- Parameters ---
    fx    : the function in which we are comparing the Lagrange interpolation with 
    p     : list of known unique points  
    x_p   : the point where the comparason takes place with fx and Px (Lagrange polynomial)
    n     : the degree of the Lagrange polynomial to be constructed for the error comparing 
    p_out : disable of enable a printout of the computed error calculations, default is true
    
    --- Returns ---
    e     : the actual error computed by the Lagrange interpolation polynomial vs. the function fx '''

    if (len(p) < n or len(p) < 2): 
        raise ValueError('Not enough known points provided.')    
    if (n < 1):
        raise ValueError('Degree has to be at least 1.')

    # actual error
    estimated = lagrange_interpolation(p, x_p, n)
    actual = fx(x_p)
    ae = sp.Abs(actual - estimated)

    if (p_out):
        # maximize the left side of the error funciton 
        max_left = fx(x)
        for i in range(n + 1): 
            max_left = sp.diff(max_left, x)

        p = sorted(p, key = lambda point : sp.Abs(x_p - point.x))
        bound = [p[0].x, p[n].x]
        max_left /= sp.factorial(n + 1)
        max_left_o = find_max_save_value(max_left, bound)

        # maximize the right size of the error funciton 
        max_right = 1
        for i in range(n + 1): 
            max_right *= (x - p[i].x) 
        max_right_o = find_max_save_value(max_right, [x_p])
       
        print(f'[Max Error]:\nme = (f^(n + 1))(xi)/(n + 1)!)*|(x - x_j)*(x - x_j+1)...|')
        print(f'   = {sp.Abs(max_left)} * |{max_right}|')   
        print(f'   = {max_left_o[1]*max_right_o[1]:.12}\n')
        print(f'[Actual Error]:\nae = |actual - estimated|\n   = |{actual} - {estimated}|\n   = {ae:.12}\n')

    return ae



def find_max_save_value(fx, values : list) -> list:
    '''
    Returns the point that achieves the maximum value on a given interval and said max value

    --- Parameters ---
    fx     : function that is being maximized
    values : inputs to the function fx

    --- Returns ---
    [max_i, min_i] : maximum point in the first index, maximum value in the next '''

    x = sp.Symbol('x')
    fx = sp.lambdify(x, fx)
    max_v = 0
    max_i = 0
    for i in values:
        curr_v = sp.Abs(fx(i))
        print(f'f({i}) = {curr_v}')
        if (curr_v > max_v):
            max_v = curr_v
            max_i = i

    return [max_i, max_v]



# funciton calls
if __name__ == "__main__":
    
    x = sp.Symbol('x')

    # setting the question 
    q = 2
    degree = 2
    
    if q == 1:
        x_p = 8.4
        known_points = [Point(8.1, 16.9441), Point(8.3, 17.56492), 
                        Point(8.6, 18.50515), Point(8.7, 18.82091)]
        fx = lambda x : sp.log(x)*x

    if q == 2:
        x_p = -1/3
        known_points = [Point(-3/4, -0.07181250), Point(-1/2, -0.02475),
                        Point(-1/4, 0.3349375), Point(0, 1.101)]
        fx = lambda x : x**3 + 4.001*x**2 + 4.002*x + 1.101

    # x_r = lagrange_interpolation(known_points, x_p, degree)
    error = lagrange_error_bound(fx, known_points, x_p, degree)
