#include<iostream>
#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>


using namespace std;

// Set the value of pi
double const PI = 4.0*atan(1.0);

// Set the value of alpha ; alpha is 0 for poission equation //
double const alpha = 1000.0;

/* *************** parameters to be altered for the equation < Start> ************************ */

// RHS (F) component is hardcoded here 

double force_vector(double x , double y, double hx, double hy)
{
    double a;
    a = -200.0*PI*PI*sin(10.0*PI*x)*cos(10.0*PI*y) ;
	a += 1000.0*( sin(10.0*PI*x)*cos(10.0*PI*y) ) ;
    a *= hx*hx ;
    return a;
            

}

// Dirichlet boundary

long double dirichlet_boundary  (double x , double y , int c )
{
    long double a;
    switch(c)
    {
         // Bottom parallel to X axis  i.e (y=0 ;  u(x,0))
        case 1: a = sin(10*PI*x);
            break;
        //parallel to Y axis - left i.e (x = 1 ; u(1,y))
        case 2: a =  0.0;
            break;
        // top parallel to X axis  i.e (y=1 ;  u(x,1))
        case 3: a = sin(10*PI*x);
            break;
        //parallel to Y axis - left i.e (x = 0 ; u(0,y))
        case 4: a =  0.0;
            break;
        default: a =  0.0;
        break;
    }
    return a;
}

/* ****** parameters to be altered for the Solution < end> **** */

// Function to set the memory for array and initialising the values to Zero

void memset_2d( double**& arr, int row , int col)
{
    arr = new double*[row] ;
    for (int i = 0 ; i < row ; i ++ )
    {
        arr[i] = new double[col];
        for (int j = 0 ; j< col ; j++ )
            arr[i][j] = 0.0;
    }
}

void matrix_multiply( double**& ans_arr , double**& arr_a , double**& arr_b, int row1, int col1,int col2)
{
    for(int i=0;i<row1;i++)
        for (int j=0;j<col2;j++)
            for (int k=0;k<col1;k++)
            ans_arr[i][k]  += (arr_a[i][j] * arr_b[j][k]);
}

void print_matrix(double**& arr1,int row , int col)
{
    for(int i=0;i<row;i++)
    {
        cout<<endl;
        for (int j=0;j<col;j++)
        {
            cout<<arr1[i][j]<<"\t";
        }
    }
    cout<<endl;
}


int main ()
{
// Set the precision flags for cout
    cout.flags( ios::dec | ios::scientific );
    cout.precision(20);
    
    //Declarations
    int Nx,Ny,i,j,k,vec_size;
    double Hx,Hy;
    // Set Initial Values 
    double x0 = 0.0;double xn  = 1.0;double y0 = 0.0;double yn = 1.0;
    
   // Get value from user for the Input 
    cout << "Enter N Value : ";
    cin >> Nx;
    
    Ny = Nx;Hx = Hy;
    Hx = (xn-x0)/double(Nx);
    Hy = (yn-y0)/double(Ny);
    
    // Declarations of variables
    double *x ,*y, *Eig_val, *norm ;
    double **Eig_vec,**uExact,**u,**Eig_vec_trans,**temp1, **force , **F_Tilde , **temp2, **U_Tilde,**temp3;
    double Maximum_error , L2_error;
    x = new double[Nx+1]; 
    y = new double[Ny+1];
    Maximum_error = L2_error = 0.0;
    
    // Set Grid points ( UNiform )
    for ( int i = 0; i <= Nx ; i++)  x[i] = i/double(Nx);   
    for ( int i = 0; i <= Ny ; i++)  y[i] = i/double(Ny);
       
    
    // Memory Allocations for the arrays
    memset_2d(Eig_vec,(Nx-1),(Ny-1) ); memset_2d(Eig_vec_trans,(Nx-1),(Ny-1) ); 
    memset_2d(u,Nx+1,Ny+1);memset_2d(uExact,Nx+1,Ny+1);
    memset_2d(force,Nx-1,Ny-1); memset_2d(F_Tilde,Nx-1,Ny-1);memset_2d(U_Tilde,Nx-1,Ny-1);
    memset_2d(temp1,Nx-1,Ny-1); memset_2d(temp2,Nx-1,Ny-1);memset_2d(temp3,Nx-1,Ny-1);
    
    Eig_val = new double[Nx-1];
    norm = new double[Nx-1];
    
    // Exact Solution for the system
    for(i = 0 ; i <= Nx ; i++)
        for(j = 0 ; j <= Ny ; j++){
            u[i][j] = 0.0 ; 
            uExact[i][j] = sin(10.0*PI*x[i])*cos(10.0*PI*y[j]) ;
	    }
    
    //Compute Eigen Vector and Eigen Value of the Matrix and calculate the force Vectors
    for( i = 1 ; i <Nx ; i++ )
    {
        Eig_val[i-1] = -4.0*sin(PI*0.5*x[i])*sin(PI*0.5*x[i]) ;
        norm[i-1] = 0.0;
        for ( j = 1 ; j < Ny ; j++)
        {
            Eig_vec[i-1][j-1] = sin(i*PI*x[j]);    
            norm[i-1] = norm[i-1]  + (Eig_vec[i-1][j-1]*Eig_vec[i-1][j-1]) ;
            force[i-1][j-1] = -200.0*PI*PI*sin(10.0*PI*x[i])*cos(10.0*PI*y[j]) ;
			force[i-1][j-1] += 1000.0*( sin(10.0*PI*x[i])*cos(10.0*PI*y[j]) ) ;
			force[i-1][j-1] *= Hx*Hx ;
            
        }
    }
//     cout<<"force"<<endl;
//     print_matrix(force,Nx-1,Ny-1);
//    cout<<"Eig Vec - before Normalisation : " <<endl;
//      print_matrix(Eig_vec,Nx-1,Ny-1);
    //Normalising Eigen Vectors 
    for( i = 1 ; i<Nx ; i++)
        for(j = 1;j<Ny;j++)
            Eig_vec[i-1][j-1] = Eig_vec[i-1][j-1]/(sqrt(norm[i-1]));
 
//     cout<<"done"<<endl;
//     
//     cout<<"Eig Vec - after Normalisation : " <<endl;
//     print_matrix(Eig_vec,Nx-1,Ny-1);
    
    //Appending boundary Value to the RHS of Solution - Dirchlet Boundary and finding INverse Eigenvalues
    for( i = 1 ; i<Nx ; i++)
        for(j = 1;j<Ny;j++) 
        {
            Eig_vec_trans[i-1][j-1] = Eig_vec[j-1][i-1];
            if(j==1)   // Bottom parallel to X axis  i.e (y=0 ;  u(x,0))
            {
                force[i-1][j-1] = force[i-1][j-1] - dirichlet_boundary(x[i],y[j-1],1);
                u[i][j-1] = dirichlet_boundary(x[i],y[j-1],1);
            }
            
            if(j==Nx-1)  // top parallel to X axis  i.e (y=1 ;  u(x,1))
            {
                force[i-1][j-1] = force[i-1][j-1] - dirichlet_boundary(x[i],y[j+1],3);
                u[i][j+1] = dirichlet_boundary(x[i],y[j+1],3);
            }
            
            if(i==1) //parallel to Y axis - left i.e (x = 0 ; u(0,y))
            {
                force[i-1][j-1] = force[i-1][j-1] - dirichlet_boundary(x[i-1],y[j],4);
                u[i-1][j] = dirichlet_boundary(x[i-1],y[j],1);
            }
            if(i==Ny-1) //parallel to Y axis - left i.e (x = 1 ; u(1,y))
            {
                force[i-1][j-1] = force[i-1][j-1] - dirichlet_boundary(x[i+1],y[j],2);
                u[i+1][j] = dirichlet_boundary(x[i+1],y[j],1);
            }
        }
    
//      cout<<"done2"<<endl;
//     cout<<" Force : "<<endl;
//     print_matrix(force,Nx-1,Ny-1); 
    
    
    // calculate F-Tilde -->   Eigvec^(-1) * Force * Eig_vec 
    // ( Note : matrix Multiply ( answer, mat1 , mat2 , row(mat1) , col(mat2) , col(mat1)))
    matrix_multiply(temp1,Eig_vec_trans,force,Nx-1,Ny-1,Nx-1);
//    cout<<" temp1  : "<<endl;
//     print_matrix(temp1,Nx-1,Ny-1);
    matrix_multiply(F_Tilde,temp1,Eig_vec,Nx-1, Ny-1,Nx-1);
//     cout<<" F_tilde  : "<<endl;
//     print_matrix(F_Tilde,Nx-1,Ny-1);

    // Calulate U-Tilde --> Force / (sum of eig_values + alpha)
    for(i=1;i<Nx;i++)
        for(j=1;j<Ny;j++)
        {
            U_Tilde[i-1][j-1]  = F_Tilde[i-1][j-1] / (Eig_val[i-1] + Eig_val[j-1] + (alpha*Hx*Hy)) ;
        }
//     cout<<"U Tilde : "<<endl;
//     print_matrix(U_Tilde,Nx-1,Ny-1);
    // Calulate U = eig_vec * u-tilde * eig_vec_trans
   
    matrix_multiply(temp2,Eig_vec,U_Tilde,Nx-1,Ny-1,Nx-1);
//     cout<<"temp2 : "<<endl;
//     print_matrix(temp2,Nx-1,Ny-1);
    
    for(i = 1 ; i < Nx ; i++)
    {
        for (j = 1; j < Ny ; j++)
        {
             u[i][j] = 0.0;
            for ( k = 1; k < Nx ; k++)
            {
                u[i][j] += (temp2[i-1][k-1] * Eig_vec_trans[k-1][j-1]);
            }
        }
    }
    
//     cout<<"u"<<endl;
// 	print_matrix(u,Nx+1,Ny+1);
//     cout<<"uexact"<<endl;
// 	print_matrix(uExact,Nx+1,Ny+1);
     
    ofstream outfile("Output.dat", ios::out) ;
    if(!outfile) {cerr<< "Error: Output file couldnot be opened.\n";}
    
    outfile.flags( ios::dec | ios::scientific);
    outfile.precision(8);
    outfile << "TITLE = Temp_distribution" <<endl;
    outfile<< " VARIABLES = x, y, u_computed, u_exact " <<endl;
    outfile << "Zone T = psi I = " << Ny+1 << " J = " << Nx+1 << endl ;
    for( i = 0 ; i <= Nx; i++)
    {
        for(j = 0 ; j <= Ny ; j++)
        {
            if( fabs(u[i][j] - uExact[i][j])  > Maximum_error)
                Maximum_error = fabs(u[i][j] - uExact[i][j]) ;
            
            L2_error = L2_error + ( ((u[i][j] - uExact[i][j])*(u[i][j] - uExact[i][j])) / ((Nx+1)*(Ny+1)) );
           outfile <<x[i] <<"\t"<<y[j] <<"\t"<< u[i][j] << "\t" << uExact[i][j] <<endl;
        }
    }
     
    
    L2_error = sqrt(L2_error);
    
    cout <<"L2_error : " << L2_error <<"     L-inf error : " <<Maximum_error << endl;
    
    outfile.close() ; 
    // Free up the heap memory 
    delete Eig_val;delete Eig_vec;delete Eig_vec_trans;delete x; delete y;delete temp1;delete temp2;
   delete u;delete uExact;
    //double **Eig_vec,**uExact,**u,**Eig_vec_trans,**temp1, **force , **F_Tilde , **temp2, **U_Tilde;
        
}
