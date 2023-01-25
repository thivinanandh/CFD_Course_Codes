![horizontal line](media/image3.png){width="6.470197944006999in" height="0.11458333333333333in"} 
------------------------------------------------------------------------------------------------

![Placeholder image](media/image19.jpg){width="6.692708880139983in"
height="3.25in"}

Finite Difference Solver for Laplacian equation using Successive Over
Relaxation (SOR) method

09.04.20XX

**─**\
Your Name\
Your Company\
123 Your Street

**[Introduction to PDE and Numerical
Methods](#introduction-to-pde-and-numerical-methods) 2**

> [Laplace Equation](#laplace-equation) 2
>
> [Numerical Solution of PDE using
> FDM](#numerical-solution-of-pde-using-fdm) 3
>
> [Taylor's Series](#taylors-series) 3
>
> [Finite Difference Formulas](#finite-difference-formulas) 3
>
> [First order Derivatives](#first-order-derivatives) 4
>
> [Second Order Derivatives](#second-order-derivatives) 4

**[Problem Solution and Technique](#problem-solution-and-technique) 5**

> **[Problem Equation](#problem-equation) 5**
>
> [Numerical Discretisation of laplace
> equation](#numerical-discretisation-of-laplace-equation) 5
>
> **[Iterative Methods](#iterative-methods) 6**
>
> [Successive Over Relaxation](#successive-over-relaxation) 6

**[Verification of Accuracy](#verification-of-accuracy) 7**

> [Convergence Rate](#convergence-rate) 7
>
> [Validation against Analytical
> Solution](#validation-against-analytical-solution) 8

**[Graphs/ Tables/ Visuals](#graphs-tables-visuals) 8**

> [Effect of Mesh size on
> convergence](#effect-of-mesh-size-on-convergence) 8
>
> [Contour plot](#contour-plot) 9
>
> [Stopping Criteria](#stopping-criteria) 10
>
> [Effect of Relaxation parameters](#effect-of-relaxation-parameters) 10

**[Code](#code) 11**

> [Grid Generation](#grid-generation) 11
>
> [Boundary Conditions](#boundary-conditions) 12
>
> [Error Metric](#error-metric) 12
>
> [SOR](#sor) 12
>
> [Exit criteria](#exit-criteria) 12
>
> [Plotting](#plotting) 13
>
> [Full code](#full-code) 13

**[Result](#result) 17**

Introduction to PDE and Numerical Methods
=========================================

A Equation is called a partial differential equation when they have an
unknown function of two or more variables and its partial derivatives
with respect to these variables The partial differential equations form
the major part of the equations which governs how physical phenomena
behave in the natural world. Some of the examples of these equations are

-   Navier-Stokes equations which governs the fluid flow,

-   Black scholes equation which governs the dynamics of physical
    > objects,

-   Poisson Equation which governs the convection of quantities like
    > heat, concentration etc

-   Maxwell\'s equations - for electromagnetics

-   Schrodinger Equation - for Quantum mechanics ( wave equation )

Unlike the Ordinary differential equations, The analytical solution does
not exist for all kinds of partial differential equations. Due to which
, we need to resort to numerical methods like "Finite Difference
schemes" in order to compute the solution of the problem.

Laplace Equation 
----------------

![](media/image10.png){width="1.3645833333333333in"
height="0.4895833333333333in"}

Laplace equation is the very basic form of the partial differential
equation which is available for analysis. The laplacian operator in
general signifies the diffusion operator.(i,e) how an unstable system
will come to an equilibrium when the external force is removed. With the
external force term ( RHS not zero ), it will be called as the
poission's equation. This can be used to represent various physical
problems like temperature distribution electrostatic potential etc.
These equations upon solved, will give the equilibrium state of the
Ungoya variable (in this case electrostatic potential ) in the entire
domain, given the boundary conditions of the given system.

￼

Numerical Solution of PDE using FDM 
-----------------------------------

A numerical technique for solving differential equations is the finite
difference method. In order to solve it, a continuous function must be
replaced with a collection of discrete values at a specific set of
points, followed by mathematical operations. When solving partial
differential equations, which are equations involving partial
derivatives of a function with respect to two or more variables, this
method is extremely helpful. It is based on the series called as
Taylor's Series

Taylor's Series
---------------

A mathematical series called the Taylor series is used to approximation
functions. It has the name of Brook Taylor, a mathematician who
initially presented the idea in the early 18th century. A polynomial
function of the function being approximated is represented by each term
in the series, which is an infinite sum of terms. The polynomial is set
up so that the first few terms of the series correspond to the
function\'s values at a specific point, while the later terms correct
the approximation. The series will approximate more accurately the more
terms it contains.

![](media/image27.png){width="4.302083333333333in"
height="1.2708333333333333in"}

Finite Difference Formulas
--------------------------

The finite difference approach employs finite difference formulas, which
are mathematical expressions used to approximate function values at
specific points. These formulas can be used to approximate the
derivatives of a function and often involve differences of the function
at different places. The order of the approximating derivative and the
number of points involved will determine the exact form of the
formula.For example, the forward difference formula is used to
approximate the first derivative of a function at a point, and is given
by:

![](media/image20.png){width="2.2604166666666665in" height="0.40625in"}

where f(x) is the function being approximated, \$x\$ is the point at
which the derivative is being approximated, and h is the step size. By
using higher-order finite difference formulas, which contain differences
of the function at more places, might increase the approximation\'s
accuracy.

First order Derivatives
-----------------------

The first order derivatives can be distinguished into

-   Forward Difference

> ![](media/image20.png){width="2.2604166666666665in"
> height="0.40625in"}

-   Backward Difference

> ![](media/image1.png){width="2.2604166666666665in" height="0.40625in"}

-   Central Difference

> ![](media/image9.png){width="2.6041666666666665in" height="0.40625in"}

Second Order Derivatives
------------------------

The second order derivatives can be represented as combination of the
first order formulas. So the most commontly used formula is the

![](media/image18.png){width="3.3958333333333335in"
height="0.4270833333333333in"}

By comparing the first derivatives at two places, one on either side of
the point being approximated, this formula estimates the second
derivative. The approximation\'s accuracy can be increased by reducing
the step size \$h\$.

It is also feasible to approximate the second derivative using
higher-order finite difference formulas. For instance, the third-order
finite difference formula known as the five-point formula can be used to
approximate the second derivative.

![](media/image7.png){width="6.21875in" height="0.4270833333333333in"}

Problem Solution and Technique
==============================

Problem Equation
----------------

We need to Solve the

![](media/image8.png){width="1.3125in" height="0.46875in"}

With the boundary conditions

![](media/image14.png){width="4.286458880139983in"
height="1.4288199912510937in"}

Numerical Discretisation of laplace equation
--------------------------------------------

![](media/image8.png){width="1.3125in" height="0.46875in"}

This equation can be approximated using a collection of discrete values
using the finite difference approach. This is accomplished by
substituting finite difference formulas for the equation\'s derivatives.
For instance, the second derivatives in the equation can be roughly
calculated using the central difference formula, giving the discretized
form shown below:

![](media/image13.png){width="6.4375in" height="0.4270833333333333in"}

This will result in system of equations which can be solved for further
analysis

This equation can be rearranged as

![](media/image16.png){width="5.072916666666667in" height="0.40625in"}

We will solve for this equation in the given code. The C will be the
unknown variable, in our case it should be **"u"**

Iterative Methods
-----------------

Mathematical problems are solved numerically using iterative methods.
They entail continually performing a set of operations on a beginning
value to produce a series of approximations that eventually lead to the
problem\'s solution. These techniques are frequently used to discover
the root of a function or to solve systems of equations.Iterative
techniques come in various varieties, such as fixed point iteration,
Newton\'s method, and Jacobi method. The specific procedures carried out
at each iteration\'s step and the circumstances under which they are
certain to converge to a solution vary between these methods.

Successive Over Relaxation 
--------------------------

It is a variation of the Gauss-Seidel approach, which consists executing
a set of operations to the values of the variables in the system of
equations repeatedly in order to arrive at a series of approximations
that converge to the solution. Omega is a further relaxation parameter
that is included. By modifying the method\'s \"over-relaxation\" to a
desired level, this parameter can increase the convergence rate. For
some kinds of equation systems, the SOR approach can be more effective
and is typically quicker than the Gauss-Seidel method. However, it is
also more sensitive to the relaxation parameter \$omega\$ choice,
therefore it is crucial to pick a suitable value to guarantee
convergence.

The formulation is given by

![](media/image22.png){width="6.5in" height="1.1805555555555556in"}

Where b is the RHS and the A is the matrix coefficients.

Verification of Accuracy 
========================

Convergence Rate
----------------

Rate of convergence is a measure of how fast the difference between the
solution point and its estimates goes to zero. We can use the number of
iterations for different matrix sizes to compute the convergence rate

  [[S.No]{.underline}](http://s.no)   H Value   N      Num Iterations   Log N         Log Iter      Convergence
  ----------------------------------- --------- ------ ---------------- ------------- ------------- -----------------
  1                                   0.1       11     56               3.459431619   5.807354922   
  2                                   0.05      21     89               4.392317423   6.475733431   0.716463372
  3                                   0.01      101    233              6.658211483   7.864186145   0.6127615312
  4                                   0.005     201    436              7.651051691   8.768184325   0.9105172942
  5                                   0.001     1001   10134            9.967226259   13.30691611   **1.959581049**

Order of convergence : 2

Validation against Analytical Solution
--------------------------------------

We have computed the l2 norm of the actual solution vs computed
solution, the error norm is given by **0.07411994497733916 (7e-2)**

![](media/image17.png){width="6.5in" height="1.3472222222222223in"}

Graphs/ Tables/ Visuals
=======================

Effect of Mesh size on convergence
----------------------------------

![](media/image6.png){width="4.054687226596675in"
height="2.7031255468066493in"}

It seems, As the h size decreases the iterations taken increases very
rapidly. The growth is similar to an O(n\^2) or O(n\^3)

  S.No   H Value   Num Iterations
  ------ --------- ----------------
  1      0.10      56
  2      0.05      89
  3      0.01      233
  4      0.005     436
  5      0.001     10134

Contour plot 
------------

Here is the contour plot of actual and the computed Solution

![](media/image24.png){width="5.557292213473316in"
height="2.778646106736658in"}

![](media/image4.png){width="5.567708880139983in"
height="2.783853893263342in"}

Stopping Criteria
-----------------

I have taken the stopping criteria as the tolerance of 1e-5, This is a
very small number for considering that the system has converged , but
not like very small so that the computations will accumulate numerical
errors

**Stopping tolerance : 1e-5**

Effect of Relaxation parameters
-------------------------------

![](media/image28.png){width="6.25in" height="4.166666666666667in"}

![](media/image25.png){width="6.25in" height="4.166666666666667in"}

Code
====

Grid Generation
---------------

![](media/image11.png){width="4.567708880139983in"
height="3.1842202537182853in"}

To generate this grid, we have used the numpy.meshgrid() function to
generate the computational grid.

Boundary Conditions
-------------------

![](media/image2.png){width="4.354166666666667in" height="1.8125in"}

Error Metric
------------

![](media/image5.png){width="6.5in" height="1.4722222222222223in"}

SOR
---

![](media/image15.png){width="6.5in" height="1.3055555555555556in"}

Exit criteria
-------------

I have taken the stopping criteria as the tolerance of 1e-5, and max
iteration as 10000

![](media/image21.png){width="5.135416666666667in"
height="0.20833333333333334in"}

Plotting
--------

![](media/image23.png){width="6.5in" height="3.111111111111111in"}

Full code
---------

**To run just do → python3 \<filename.py\>**

\# %%

import numpy as np

import matplotlib.pyplot as plt

import sys

\# %%

\#\# Function for Computing Error Metric

def ErrorMetric(input\_array,N,p):

\"\"\"ErrorMetric

Args:

input\_array (np.array): Input numpy array for which we need to compute
Error Metric

N (int): Dimension of input numpy array

p (int): value of the pth norm to be computed

Returns:

\_type\_: \_description\_

\"\"\"

\#if its an Matrix, flatten it to vector

flatenedArray = input\_array.reshape(-1,)

\#\# loop through all elements in the flattened array to compute the

\#\# p-norm of the array

sumValue = 0.0

for value in flatenedArray:

sumValue += abs(value)\*\*p

return sumValue\*\*(1.0/p)

\#\# Function to Apply boundary condition to the input Array

def bdyVal(option,w):

\"\"\"Boundary Value Function

Args:

option (int): Type of boundary function to be returned

w (float): co-ordinates of the location where the boundary value is
needed

Returns:

numpyArray : Numpy array of all the filled boundary values

\"\"\"\'\'

\#\# setup functions

def f1(x):

return np.exp(np.pi\*x)

def f2(x):

return -1.0\*np.exp(np.pi\*x)

def g1(y):

return np.cos(np.pi \* y)

def g2(y):

return np.exp(np.pi)\*np.cos(np.pi \* y)

if(option == 1):

return f1(w)

elif(option ==2):

return f2(w)

elif(option ==3):

return g1(w)

elif(option ==4):

return g2(w)

else:

print(\"Invalid Option\")

\# %%

def SOR(A,relaxationParameter,maxiter,tolerance,p):

\"\"\"SOR Method

Args:

A (numpy array): Input numpy array

relaxationParameter (float): Relax parameter for SOR

maxiter (int): Maximum iteration for SOR

tolerance (float): exit tolerance

\"\"\"

residualArr = \[\]

errorOld = 0;

residual = 100 \#\# Entry Criteria

iterations = 0

size = A.shape\[0\];

h = 1/(size-1);

while(iterations \< maxiter and residual \> tolerance):

for i in range(1,size-1):

for j in range(1,size-1):

currUpdate = 0.25 \* (A\[i,j+1\] + A\[i,j-1\]+ A\[i-1,j\] + A\[i+1,j\]);

A\[i,j\] = relaxationParameter\* currUpdate +
(1.0-relaxationParameter)\*A\[i,j\]

error = ErrorMetric(A,A.shape,2)

residual = abs(error - errorOld)

errorOld = error

if(iterations % 100 ==0):

print(f\"At Iter {iterations} residual is : {residual}\")

residualArr.append(residual)

iterations += 1

return A,residualArr,iterations

\# %%

N = 64

relax = 1.8

x = np.linspace(0,1,N) \#\# 1d array of x

y = np.linspace(0,1,N) \#\# 1D array of y

\#\# Allocate Space for the Matrix , which stores the linear system of
Equations

A = np.zeros((N,N))

\#\# Apply boundary condition

\# Boundary Condition - 1 y=0 \-- Bottom Row

\# Pick the last row using -1 and update every value of x in every
column of the row with the function

A\[0,:\] = bdyVal(1,x)

\# Boundary Condition - 2 y=1 \-- Top Row

\# Pick the first row using 0 and update every value of x in every
column of the row with the function

A\[-1,:\] = bdyVal(2,x)

\# Boundary Condition - 3 x=0 \-- left Column

\# Pick the first column using 0 and update every value of y in the
every row of the column with the function

A\[:,0\] = bdyVal(3,y)

\# Boundary Condition - 4 x=1 \-- Right Column

\# Pick the last column using -1 and update every value of y in the
every row of the column with the function

A\[:,-1\] = bdyVal(4,y)

\#\#Send the matrix for Solving using SOR

result, resArray,iter = SOR(A,relax,10000,1e-4,2)

print(f\"SOR converged in {iter} iterations\")

Result
======

We could observe that the

-   relaxation parameter plays a vital role in the convergence of the
    > algorithm

-   Atleast in this case, the higher the relaxation ,the higher the
    > convergence

-   However, if we get very close to 2( 1.9999 ), as per the
    > literatures, the method might diverge

-   The h size has significant impact on the computation time of the
    > code.
