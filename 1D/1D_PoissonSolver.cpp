#include <iostream>
#include<cmath>

// POISSON EQUATION 1D Solver with Dirchlet Boundary //

using namespace std;

// Enter your Function Here 
static int beta = 80;
static double pi = 3.14159265358979323;

double function_val(double x)
{
    double u = - sin(pi*beta*x) ;
    return u;
}


double exact_function_val( double x)
{
    double y =  (1/(beta*beta*pi*pi)) * ( sin(beta*pi*x) );
    return y;
}


double* dirichlet_boundary_values ()
{
    static double a[1];
    a[0] = 0.0;
    a[1] = 0.0;
    return a;
}

int main()
{
    int N;
        //*********************************************//
    // TRIDIAGONAL MATRIX   
    // from i =  2 to i = N-2
    //       c*u[i+1] - b*u[i] + a*u[i-1] =  h^2 * d       
    //  [ b1 c1      ]    =  | h*h*d1 |
    //  [ a2 b2 c2   ]    =  | h*h*d2 |
    //  [       b3 c3]    =  | h*h*d3 |
    //*********************************************//
    
    // Declaration of equation parameters
    double* b ;
    b = new double[N]; 
    double* a ;
    a = new double[N];
    double* c ;
    c = new double[N];
    double* d ;
    d = new double[N];
    double* u ;
    u = new double[N];
    
    // Declaration of x and u_exact error
    double* x ;
    x = new double[N];
    double* u_exact ;
    u_exact = new double[N];
    double error;
    
    // Declaration of variables for forward substitution
    double* p ;
    p = new double[N];
    double* q ;
    q = new double[N];
    
    
   cout << " Enter the Value of N : " ;
   cin >> N ;
    
//     cout << "Enter Starting x value : " ;
//     cin >> x[0] ;
//     cout << " Enter the ending  x value : " ;
//     cin >> x[N];
    x[0] = 0;
    x[N] = 1;
    
    
    
    double h =  (x[N] - x[0])/N ; 
    cout<<endl<<"H value is : "<<h<<endl;

    for (int i =1; i<N ; i++)
    {
            x[i] = i*h;
            d[i] = (h*h)*function_val(x[i]);
            b[i] = -2;
            c[i] = 1;
            a[i] = 1;  
            u_exact[i] = exact_function_val(x[i]);
    }
    
    a[1] = 0;
    c[N-1] = 0;

    
    // Due to dirichlet boundry at the ends , we will only have nodes from 1 to N-1
    // b is diagonal element , c - super diagonal ; a - sub diagonal; d = rhs 
    
    // get the Dirichlet boubdary values into Array 
    
    double *bound_val;
    bound_val =  dirichlet_boundary_values();
    u[0] = u_exact[0] = bound_val[0] ;
    u[N] = u_exact[N] = bound_val[1] ;
    
    // **********   Thomas Algorithm *********** //
    
    //   Forward Substitution 
    for (int i = 1 ; i <= N-1 ; i++)
    {
        if(i == 1)
        {
            p[i] = c[i]/b[i];
            q[i] = d[i]/b[i];
        }
        
        double temp = b[i] - (a[i] * p[i-1]);
        p[i] = c[i] / temp;
        q[i] = ( d[i] - (a[i]*q[i-1]) ) / temp;
    }
    // Set value f U[N-1]
    u[N-1] = q[N-1];
    
    //backward Substitution
    for ( int i = N-2 ; i > 0 ; i--)
    {
        u[i] = p[i]*u[i+1] + q[i];
    }
    
//     cout<<endl;
    
    //***** end of Thomas Algorithm ******* //
    
    // Error Calculation //

    cout<<"Beta Value : "<< beta <<endl;
    
    error = 0.0 ;      // initialising error //
     for (int i = 1; i <= N-1 ; i++)
    {
      cout <<u[i] << " \t\t"  << u_exact[i] <<"\t\t" << (u[i] - u_exact[i])<<endl;
       //cout <<u[i] << " , " ; 
       error = error + ((u[i] - u_exact[i])*(u[i] - u_exact[i]));
        
    }
    
    
    double error_norm =  sqrt(error/(N-2));
    //cout<<"Error : " << error<<endl;
//     cout<<"u_exact :  " << u_exact[0]<<endl;
//     cout<< "u : "<<u[0]<<endl;
    cout << " Error Norm : " << error_norm <<endl;
    
//     cout<<"U - Calculated : "<<endl;
//     //print Functions 
//     for (int i = 0; i < N ; i++)
//     {
//        cout <<u[i] << " , ";
//        
//         
//     }
 
   
    
//     cout<<endl;

    
        
//     delete [] b;
//     delete [] a;
//     delete [] c;
//     delete [] d; 
//     delete [] p; 
//     delete [] q;
//     
//   
//     delete u;
//     delete x;
//     delete u_exact;

   
   
    
    
    return 0;
}
