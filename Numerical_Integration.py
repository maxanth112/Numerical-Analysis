import sympy as sp
import numpy as np
from functools import reduce
from Lagrange_Interpolating_Polynomials import find_max_save_value

def trapezoid_approximation(interval : list, func, n = 2) -> float:

    func = sp.lambdify(x, func)
    a, b = interval[0], interval[1]
    delta_x = (b - a)/n
    terms = [func(a)] + [2*func(a + i*delta_x) for i in range(1, n)] + [func(b)]
    
    return sum([delta_x*term for term in terms])/2


def simpsons_approximation(interval : list, func, n = 3) -> float:

    a, b = interval[0], interval[1]
    delta_x = (b - a) / n
    func = sp.lambdify(x, func)
    ans = func(a) + func(b)
    
    for i in range(1, n):
        k = a + i*delta_x        
        if i % 2 == 0:
            ans += 2*func(k)
        else:
            ans += 4*func(k)
    ans *= delta_x/3
    return ans


def integral_error(interval : list,  func, estimated : float, n : int, rule = 'trap'):

    a, b = interval[0], interval[1]
    integral = sp.integrate(func, x)
    fx_integrated = sp.lambdify(x, integral)

    if rule == 'trap':
        df = sp.diff(func, x)
        max_e = (((b - a)**3)/12)*df
        [error_index, error_bound] = find_max_save_value(max_e, [a, b])
    else: 
        df = sp.diff(sp.diff(sp.diff(func, x), x), x)
        max_e = (((b - a)**5)/12)*df
        [error_index, error_bound] = find_max_save_value(max_e, [a] + [a + i*((b - a)/(n - 1)) for i in range(1, n)]+ [b])
 
    actual = fx_integrated(b) - fx_integrated(a)    
    error_actual = actual - estimated
    
    print(f'Given a = {a:>.5f}, b = {b:>.5f}, fx = {func}\n  Compute integral(fx) from a to b => {integral} from a to b => {estimated:>.5f}')
    print(f'Errors:\n  Actual definite integral value: {actual:>.5f}\n  Estimated value: {estimated:>.5f}')
    print(f'  Actual Error = actual integral - estimated = {actual:>.5f} - {estimated:>.5f} = {error_actual:>.5f}  \nError Bound = max[|(((b - a)^3)/12)*f"(xi)|] for xi in [{a:>.5f}, {b:>.5f}]')
    print(f'   choose xi = {error_index:>.5f} => {error_bound:>.5f}')





# funciton call examples 
if __name__ == "__main__":
    
    # setting the question 
    x = sp.Symbol('x')
    q = 1

    if q == 1: 
        interval = [np.exp(1), np.exp(1) + 1]
        fx = 1/(x*sp.log(x))

    # p = trapezoid_approximation(interval, fx, 8)
    p = simpsons_approximation(interval, fx, 3)
    e = integral_error(interval, fx, p, 8, 'simp')


    
