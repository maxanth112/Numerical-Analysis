import sympy as sy

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
    
    if (len(p) < (n + 1) or len(p) < 2): 
        print(f'Not enough known points provided for a lagrange interpolation of degree ${n}')
    
    p = close_points_rearrangement(x_p, p)
    Px = 0.0
    polynomial = ['']*(n + 1)
    
    for i in range(n + 1):
        
        L_i = 1
        for j in range(n + 1):
            if i != j:
                L_i *= (x - p[j].x) / (p[i].x - p[j].x)
                polynomial[i] += f'((x - {p[j].x})/({p[i].x} - {p[j].x}))'
                
        Px += p[i].fx*L_i
        polynomial[i] = f'({p[i].fx:.5f})' + polynomial[i]

    # the final polynomial and answer
    P_x = sy.lambdify(x, Px)
    x_r = P_x(x_p)

    # printing out the constructed Lagrange polynomial and problem specifics
    if p_out:
        print('\n[Given points]: ')
        for i in range(len(polynomial)):
            print(f'x_{i} = {p[i].x}   f(x_{i}) = {p[i].fx}')

        print(f'\n[Simplified]:')
        sy.pprint(f'P_{n}(x) = {sy.expand(Px)}')
        
        for i in range(len(polynomial)):
            start = f'\n[Expanded]: \nP_{n}(x) = ' if i == 0 else ' '*len('P_n(x) = ')
            print(f'{start} {polynomial[i]} + ')
        
        print(f'\n[Estimation]:\nP_{n}({x_p}) = {P_x(x_p)}\n')
        
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

    if (len(p) < (n + 1) or len(p) < 2): 
        print(f'Not enough known points provided for a lagrange interpolation of degree ${n}')
    p = close_points_rearrangement(x_p, p)

    # actual error
    estimated = lagrange_interpolation(p, x_p, n)
    actual = fx(x_p)
    ae = sy.Abs(actual - estimated)

    if (p_out):
        # maximize the left side of the error funciton 
        max_left = fx(x)
        for i in range(n + 1): 
            max_left = sy.diff(max_left, x)
        bound = [p[0].x, p[n].x]
        max_left /= sy.factorial(n + 1)
        max_left_o = find_max_save_value(max_left, bound)

        # maximize the right size of the error funciton 
        max_right = 1
        for i in range(n + 1): 
            max_right *= (x - p[i].x) 
        max_right_o = find_max_save_value(max_right, [x_p])
       
        print(f'[Max Error]:\nme = (f^(n + 1))(xi)/(n + 1)!)*|(x - x_j)*(x - x_j+1)...|')
        print(f'   = {sy.Abs(max_left)} * |{max_right}|')   
        print(f'   = {max_left_o[1]*max_right_o[1]}\n')
        print(f'[Actual Error]:\nae = |actual - estimated|\n  = |{actual} - {estimated}|\n  = {ae}\n')

    return ae

def find_max_save_value(fx, values : list) -> list:
    '''
    Returns the point that achieves the maximum value on a given interval and said max value

    --- Parameters ---
    fx     : function that is being maximized
    values : inputs to the function fx

    --- Returns ---
    [max_i, min_i] : maximum point in the first index, maximum value in the next '''
    
    fx = sy.lambdify(x, fx)
    max_v = 0
    max_i = 0
    for i in values:
        curr_v = sy.Abs(fx(i))
        if (curr_v > max_v):
            max_v = curr_v
            max_i = i

    return [max_i, max_v]


def close_points_rearrangement(x_p : float, p : list) -> list:
    '''
    Returns the list of points sorted by closest to the approximating point x_p to furthest, 
    this is used to enhance Lagrange polynomial interpolation when the degree of the polynomial 
    is low, thus we use the optimal points first

    --- Parameters ---
    x_p : approximation point 
    p   : list of points in no specific order

    --- Returns ---
    p_s : the list p sorted by relevancy '''
    
    p_s = sorted(p, key = lambda point : sy.Abs(x_p - point.x))
    return p_s




# funciton calls
if __name__ == "__main__":
    
    # setting the question 
    q = 2
    x = sy.Symbol('x')
    
    if q == 1:
        # finding the interpolating polynomial and the error
        # P(x)
        x_p = 8.4
        degree = 1
        known_points = [Point(8.1, 16.9441), Point(8.3, 17.56492), 
                        Point(8.6, 18.50515), Point(8.7, 18.82091)]
        # error
        fx = lambda x : sy.log(x)*x
        error = lagrange_error_bound(fx, known_points, x_p, degree)

    if q == 2:
        # P(x)
        x_p = -1/3
        degree = 3
        known_points = [Point(-3/4, -0.07181250), Point(-1/2, -0.02475),
                        Point(-1/4, 0.3349375), Point(0, 1.101)]
        # error
        fx = lambda x : x**3 + 4.001*x**2 + 4.002*x + 1.101
        error = lagrange_error_bound(fx, known_points, x_p, degree)
