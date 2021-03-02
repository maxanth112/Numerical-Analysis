import sympy as sy

class Point: 
    def __init__(self, x, fx):
        self.x = x
        self.fx = fx
        
def lagrange_interpolation(p : list, x_p : float, n : int, p_out = True) -> float:
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
    
    if (len(p) < n): 
        print(f'Not enough known points provided for a lagrange interpolation of degree ${degree}')
    
    x = sy.Symbol('x')
    Px = 0.0
    polynomial = ['']*n
    
    for i in range(n):
        
        L_i = 1
        for j in range(n):
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
        print('\nGiven points: ')
        for i in range(len(polynomial)):
            print(f'x_{i} = {p[i].x}   f(x_{i}) = {p[i].fx}')

        print(f'\nWe construct the Lagrange interpolated polynomial of degree {n} as:\n[Simplified]:')
        sy.pprint(f'P(x) = {sy.expand(Px)}')
        
        for i in range(len(polynomial)):
            start = '\n[Expanded]: \nP(x) = ' if i == 0 else ' '*len('P(x) = ')
            print(f'{start} {polynomial[i]} + ')
        
        print(f'\nP({x_p}) = {P_x(x_p)}\n')
        
    return x_r



# funciton calls
if __name__ == "__main__":
    
    x_p = 8.4
    degree = 4
    known_points = [Point(8.1, 16.9441), Point(8.3, 17.56492), 
                    Point(8.6, 18.50515), Point(8.7, 18.82091)]

    fx_p = lagrange_interpolation(known_points, x_p, degree)

