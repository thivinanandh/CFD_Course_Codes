# %% 
import numpy as np
import matplotlib.pyplot as plt
import sys
# %%
## Function for Computing Error Metric 

def ErrorMetric(input_array,N,p):
    """ErrorMetric

    Args:
        input_array (np.array): Input numpy array for which we need to compute Error Metric
        N (int): Dimension of input numpy array
        p (int): value of the pth norm to be computed

    Returns:
        _type_: _description_
    """
    #if its an Matrix, flatten it to vector
    flatenedArray = input_array.reshape(-1,)

    ## loop through all elements in the flattened array to compute the 
    ## p-norm of the array
    sumValue = 0.0
    for value in flatenedArray:
        sumValue += abs(value)**p
    
    return sumValue**(1.0/p)    


## Function to Apply boundary condition to the input Array
def bdyVal(option,w):
    """Boundary Value Function

    Args:
        option (int): Type of boundary function to be returned
        w (float): co-ordinates of the location where the boundary value is needed

    Returns:
        numpyArray  : Numpy array of all the filled boundary values
    """''
    ## setup functions 
    def f1(x):
        return np.exp(np.pi*x)
    def f2(x):
        return -1.0*np.exp(np.pi*x)
    def g1(y):
        return np.cos(np.pi * y)
    def g2(y):
        return np.exp(np.pi)*np.cos(np.pi * y)


    if(option == 1):
        return f1(w)
    elif(option ==2):
        return f2(w)
    elif(option ==3):
        return g1(w)
    elif(option ==4):
        return g2(w)
    else:
        print("Invalid Option")


# %%


def SOR(A,relaxationParameter,maxiter,tolerance,p):
    """SOR Method

    Args:
        A (numpy array): Input numpy array
        relaxationParameter (float): Relax parameter for SOR
        maxiter (int): Maximum iteration for SOR
        tolerance (float): exit tolerance
    """

    residualArr = []
    
    errorOld = 0;
    residual = 100  ## Entry Criteria
    iterations = 0
    size = A.shape[0];
    h = 1/(size-1);
    while(iterations < maxiter and residual >  tolerance):
        for i in range(1,size-1):
            for j in range(1,size-1):
                currUpdate =  0.25 * (A[i,j+1] + A[i,j-1]+ A[i-1,j] + A[i+1,j]);
                A[i,j] = relaxationParameter* currUpdate + (1.0-relaxationParameter)*A[i,j]

        error = ErrorMetric(A,A.shape,2)
        residual = abs(error - errorOld)
        errorOld = error
        if(iterations % 100 ==0):
            print(f"At Iter {iterations}  residual is : {residual}")
        residualArr.append(residual)
        iterations += 1
    return A,residualArr,iterations

# %%
N = 64
relax = 1.8

x = np.linspace(0,1,N)    ## 1d array of x
y = np.linspace(0,1,N)    ## 1D array of y


## Allocate Space for the Matrix , which stores the linear system of Equations 
A = np.zeros((N,N))

## Apply boundary condition

# Boundary Condition - 1 y=0  -- Bottom Row
# Pick the last row using -1 and update every value of x in every column of the row with the function
A[0,:] = bdyVal(1,x)

# Boundary Condition - 2 y=1  -- Top Row
# Pick the first row using 0 and update every value of x in every column of the row with the function
A[-1,:] = bdyVal(2,x)

# Boundary Condition - 3 x=0  -- left Column 
# Pick the first column using 0 and update every value of y in the every row of the column with the function
A[:,0] = bdyVal(3,y)

# Boundary Condition - 4 x=1  -- Right Column
# Pick the last column using -1 and update every value of y in the every row of the column with the function
A[:,-1] = bdyVal(4,y)

##Send the matrix for Solving using SOR 
result, resArray,iter = SOR(A,relax,10000,1e-4,2)

print(f"SOR converged in {iter} iterations")



# %%

h = 0.01
### FOr plotting residual 
relaxataion = [1.0,1.2,1.4,1.6,1.8,1.9]
ITERARRAY = []
RESARRAY = []
for r in relaxataion:
    N = int(1/h) + 1
    A = np.zeros((N,N))
    ## generate the 2D Grid for the points x, y 
    x = np.linspace(0,1,N)    ## 1d array of x
    y = np.linspace(0,1,N)    ## 1D array of y
    A[0,:] = bdyVal(1,x)
    A[-1,:] = bdyVal(2,x)
    A[:,0] = bdyVal(3,y)
    A[:,-1] = bdyVal(4,y)
    print("----------------------------------------")
    print("Running SOR for relaxation value : ", r ," For h : " , h)
    print("----------------------------------------")
    result, resArray,iter = SOR(A,r,10000,1e-4,2)
    ITERARRAY.append(iter)
    RESARRAY.append(resArray)

with open("ITER_RELAX",'wb') as f:
    np.save(f,ITERARRAY)

with open("RESARRAY_RELAX",'wb') as f:
    np.save(f,RESARRAY)

with open("RESULT_RELAX",'wb') as f:
    np.save(f,result)

# %%
markArr = ['v','o',"^",'<','>','s','p','h']
fig,ax = plt.subplots(1,1, figsize=(6,4))
for i,resArr in enumerate(RESARRAY):
    ax.plot(resArr[10:80],'-',label=relaxataion[i],marker=markArr[i], markevery=10)
ax.set_title('Convergence for various Relaxation Values')
ax.set_xlabel("Iterations")
ax.set_ylabel("Residual")
ax.legend()
fig.savefig('convergence.png')

fig,ax = plt.subplots(1,1, figsize=(6,4))

ax.plot(relaxataion,ITERARRAY,'r-o',label="iterations")
ax.set_title('Iteration for various Relaxation Values')
ax.set_xlabel("Relaxation Values")
ax.set_ylabel("Iterations")
ax.legend()
ax.grid()
fig.savefig('iterations_for_relax.png')


# %%

### FOr plotting residual 
hArr = [0.1,0.05,0.01,0.005,0.001]
relax = 1.9
ITERARRAY = []
RESARRAY = []
for h in hArr:
    N = int(1/h) + 1
    x = np.linspace(0,1,N)    ## 1d array of x
    y = np.linspace(0,1,N)    ## 1D array of y
    A = np.zeros((N,N))
    A[0,:] = bdyVal(1,x)
    A[-1,:] = bdyVal(2,x)
    A[:,0] = bdyVal(3,y)
    A[:,-1] = bdyVal(4,y)
    print("----------------------------------------")
    print("Running SOR for H value : ", h, " With Relax : " , relax)
    print("----------------------------------------")
    result, resArray,iter = SOR(A,r,10000,1e-4,2)
    ITERARRAY.append(iter)

x = np.linspace(0,1,N)    ## 1d array of x
y = np.linspace(0,1,N)    ## 1D array of y
X,Y = np.meshgrid(x,y)

with open("ITER_H",'wb') as f:
    np.save(f,ITERARRAY)

with open("RESARRAY_H",'wb') as f:
    np.save(f,RESARRAY)

with open("RESULT_H",'wb') as f:
    np.save(f,result)

exact_A = np.zeros((N,N))
exact_A = np.exp(np.pi * X)*np.cos(np.pi * Y)

fig,ax = plt.subplots(1,2, figsize=(10,5))
ax[0].imshow(exact_A,origin='lower',cmap='jet')
ax[0].axis('off')
ax[0].axis("tight")
ax[0].set_title('Analytical Solution')

ax[1].imshow(result,origin='lower',cmap='jet')
ax[1].axis('off')
ax[1].axis("tight")
ax[1].set_title('Computed Solution')
fig.suptitle(f"Actual vs Computed Solution - SOR Relaxation {relax}")
fig.savefig("Actual vs Computed.png",dpi=300)


fig,ax = plt.subplots(1,2, figsize=(10,5))
ax[0].contourf(X,Y,exact_A)
ax[0].contour(X,Y,exact_A)
ax[0].axis('off')
# ax[0].axis("tight")
ax[0].set_title('Analytical Solution')

ax[1].contourf(X,Y,result)
ax[1].contour(X,Y,result)
ax[1].axis('off')
# ax[1].axis("tight")
ax[1].set_title('Computed Solution')
fig.suptitle(f"Actual vs Computed Solution - SOR Relaxation {relax}")
fig.savefig("Actual vs Computed_contour.png",dpi=300)


fig,ax = plt.subplots( figsize=(6,4))
ax.plot(hArr,ITERARRAY,'r-o',label="iterations")
ax.set_title('Iteration for various h')
ax.set_xlabel("h")
ax.set_ylabel("Iterations")
ax.legend()
ax.grid()
fig.savefig('iterations_for_h.png')