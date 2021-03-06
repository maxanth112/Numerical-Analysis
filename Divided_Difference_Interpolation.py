import sympy as sp

def newton_divided_difference(points : list):
    ''' Max Wiesner 
        3/4/21

    --- Required ---
    Must have at least n + 1 points supplied for a polynomial interpolation of degree n, 
    the points must be unique

    --- Parameters ---
    points : a 2D list of points in the format of [[x1, fx1], [x2, fx2]]

    --- Returns ---
    None   : is used to generate a polynomial of degree n, and print the polynomial to the terminal '''

    x = sp.Symbol('x')
    n = len(points)
    F = [ ( [points[j][1]] + [0]*n ) for j in range(n) ]

    for i in range(1, n):  
        for j in range(n - i):  
            F[j][i] = ((F[j][i - 1] - F[j + 1][i - 1]) / (points[j][0] - points[i + j][0])); 

    Px = F[0][0]
    for i in range(1, n):
        px_i = 1
        for j in range(i):
            px_i *= (x - points[j][0]) 
        Px += F[0][i]*px_i
    print(f'\nGenerated Polynomial:\nP_{n - 1} = {Px}\n')

    return None



# funciton calls
if __name__ == "__main__":

    points = [ [1, 2], [2, 3], [3, 5] ]
    ans = newton_divided_difference(points)
