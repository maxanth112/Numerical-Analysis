
def proterm(i, value, x):  
    pro = 1;  
    for j in range(i):  
        pro = pro * (value - x[j]);  
    return pro;  
  

def dividedDiffTable(x, y): 

    n = len(x)
    F = [ ( [y[j]] + [0]*n ) for j in range(n) ]

    for i in range(1, n):  
        for j in range(n - i):  
            F[j][i] = ((F[j][i - 1] - F[j + 1][i - 1]) / (x[j] - x[i + j])); 
 
    return F; 
  

def applyFormula(value, x, y):  
  
    sum = y[0][0];  
    n = len(x)
  
    for i in range(1, n): 
        sum = sum + (proterm(i, value, x) * y[0][i]);  
      
    return sum;  
  
 
def printDiffTable(y):  
    n = len(y)
    for i in range(n):  
        for j in range(n - i):  
            print(round(y[i][j], 4), "\t",  
                               end = " ");  
  
        print("");  
  
# Driver Code 
y = [0 for j in range(3)];  
x = [0, 0.25, 0.5];  
    
y[0] = 1;  
y[1] = 1.64872;  
y[2] = 2.71828;  

y=dividedDiffTable(x, y);  
printDiffTable(y);  
  
 
value = 7;  

print("\nValue at", value, "is", 
        round