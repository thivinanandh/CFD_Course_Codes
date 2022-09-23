# CFD Code Repo   
---

This repo contains the basic CFD Codes for 1D and 2D Problem (Stationary and time dependent problem)


## 2D Elliptic and Parabolic PDE's
---

$u_t - \epsilon \Delta u = f$

### Heat Transfer Problem 

$u_t - \epsilon \Delta u = f$ <br />

$u = 0$ at $y=1$


<img height="200px" width="300px" src="Images/2d-1.jpeg">


### Periodic Forcing Function

$u_t - \epsilon \Delta u = f$ <br />

$u = 0$ on all Boundaries

$f(x,y,t) = sin( 2\pi x) * sin(2 \pi y) * sin(2 \pi t)$

<img height="200px" width="300px" src="Images/2d_2.jpeg">

<img height="200px" width="300px" src="Images/2d-3.jpeg">


---
---

## 1D hyperbolic PDE
---

$u_t + \alpha u_x  =f$

The stability of the following schemes have beeb studied for various CFL numbers
* Crank Nicholson
* FTBS - Forward in Time , Backward in Space
* BTBS - Backward in Time, Backward in Space (Implicit)
* CTCS - Central in Time Central in Space (Leap Frog)



### **CFL = 0.5**
---


#### Crank Nicholson
<img height="200px" width="300px" src="Images/1d-1.png">

#### FTBS
<img height="200px" width="300px" src="Images/1d-2.png">

#### CTCS
<img height="200px" width="300px" src="Images/1d-3.png">

#### BTBS
<img height="200px" width="300px" src="Images/1d-1.png">

#### Summary :


<img src="Images/1DTable.png" />

### **CFL = 1**
---


#### Crank Nicholson
<img height="200px" width="300px" src="Images/1d-5.png">

#### FTBS
<img height="200px" width="300px" src="Images/1d-6.png">

#### CTCS
<img height="200px" width="300px" src="Images/1d-7.png">

#### BTBS
<img height="200px" width="300px" src="Images/1d-8.png">

#### Summary :


<img src="Images/1d-Table2.png" />

### **CFL = 2**
---


#### Crank Nicholson
<img height="200px" width="300px" src="Images/1d-9.jpeg">

#### FTBS
<img height="200px" width="300px" src="Images/1d-10.png">

#### CTCS
<img height="200px" width="300px" src="Images/1d-11.png">

#### BTBS
<img height="200px" width="300px" src="Images/1d-12.png">

#### Summary :


<img src="Images/1dTable3.png" />


### **Overall Plot**


<img height="200px" width="300px" src="Images/1d-13.png">
<img src="Images/1dTable4.png" />

