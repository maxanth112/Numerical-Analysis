import math


def false_position(f, p_0, p_1, N, e):
    ''' Max Wiesner 
        2/16/21

    Approximates the root of a given function f using two points, 
    p_0 and p_1, iterating using the the method of false position.

    --- Criteria ---
    The given points p_0 and p_1 must intersect the x axis.

    --- Parameters ---
    f   : the funciton 
    p_0 : initial point (first)
    p_1 : initial point (second)
    N   : maximum amount of iterations allowed 
    e   : tolerance, maximum error allowed 

    --- Returns ---
    approximation of the root p within the given tolerance. '''

    i = 2
    q_0 = f(p_0)
    q_1 = f(p_1)
    while i <= N:
        p = p_1 - ((p_1 - p_0) * q_1)/(q_1 - q_0)
        print(f'p_{str(i).ljust(2)} = {p:8.14f}')

        if abs(p - p_1) < e:
            print(f'\np found after {i} iterations.')
            return p
        i += 1
        q = f(p)

        if q * q_1 < 0:
            p_0 = p_1
            q_0 = q_1
        else:
            p_1 = p
            q_1 = q
    print(f"\nThe method failed after {N} iterations. ")
    return None


# calling with lambda functions
def f(x): return x**7 - 2*x**2 + 3


false_position(f, 1, 2, 20, 10e-5)
