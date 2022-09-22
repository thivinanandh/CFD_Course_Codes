#include<iostream>
#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>


using namespace std;

// Set the value of pi
double const PI = 4.0*atan(1.0);

// Set value of B
double const B1 = 1.0;

// Comment out this line to Eliminate the print statements in the file or un comment it //
//#define print 0

/************************************************************************************************/
/* *************** parameters to be altered for the equation < Start> ************************ */

// RHS (F) component is hardcoded here 

void force_vector(double& a,double x , double y, double hx, double hy, double t)
{
   a = sin(2.0*PI*x)*sin(2.0*PI*y)*sin(2.0*PI*t)  ;
    a *= hx*hy ;
}

// Dirichlet boundary values
 double dirichlet_boundary  (double x , double y , int c )
{
    double a;
    switch(c)
    { 
        case 1: a = 0.0;       // Bottom parallel to X axis  i.e (y=0 ;  u(x,0))
            break;  
        case 2: a =  0.0;               //parallel to Y axis - left i.e (x = 1 ; u(1,y))
            break;
        case 3: a = 0.0;      // top parallel to X axis  i.e (y=1 ;  u(x,1))
            break;
        case 4: a =  0.0;              //parallel to Y axis - left i.e (x = 0 ; u(0,y))
            break;
        default: a =  0.0;
        break;
    }
    return a;
}

/* *********************** parameters to be altered for the problem < end> ************************* */
/************************************************************************************************* */

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

//Matrix Multiplication Function
void matrix_multiply( double**& ans_arr , double** arr_a , double** arr_b, int row1, int col1,int col2)
{
    for(int i=0;i<row1;i++)
        for (int j=0;j<col2;j++){
            ans_arr[i][j] = 0.0;
            for (int k=0;k<col1;k++)
            ans_arr[i][j]  += (arr_a[i][k] * arr_b[k][j]);
        }
}

//Function to print matrix
void print_matrix(double**& arr1,int row , int col)
{
    for(int i=0;i<row;i++)
    {
        cout<<endl;
        for (int j=0;j<col;j++){
            cout<<arr1[i][j]<<"\t";
        }
    }
    cout<<endl;
}

void force_calculator(double**& arr, double*& x, double*& y, int Nx , int Ny , double hx, double hy,double**& u,double t)
{
    double b;
    for(int i = 1 ; i <Nx ; i++ )
        for (int j = 1 ; j < Ny ; j++)
        {
            force_vector(b,x[i],y[j],hx,hy,t);
            arr[i-1][j-1] = b;
            // Boundary Conditions being appended to the force Matrix ( Do not change Values here )
            if(j==1)  {                     // Bottom parallel to X axis  i.e (y=0 ;  u(x,0))
                arr[i-1][j-1] = arr[i-1][j-1] + dirichlet_boundary(x[i],y[j-1],1);
                u[i][j-1] = dirichlet_boundary(x[i],y[j-1],1);
            }
            
            if(j==Nx-1) {                   // top parallel to X axis  i.e (y=1 ;  u(x,1))
                arr[i-1][j-1] = arr[i-1][j-1] + dirichlet_boundary(x[i],y[j+1],3);
                u[i][j+1] = dirichlet_boundary(x[i],y[j+1],3);
            }
            
            if(i==1) {                      //parallel to Y axis - left i.e (x = 0 ; u(0,y))
                arr[i-1][j-1] = arr[i-1][j-1] + dirichlet_boundary(x[i-1],y[j],4);
                u[i-1][j] = dirichlet_boundary(x[i-1],y[j],1);
            }
            if(i==Ny-1) {                   //parallel to Y axis - left i.e (x = 1 ; u(1,y))
                arr[i-1][j-1] = arr[i-1][j-1] + dirichlet_boundary(x[i+1],y[j],2);
                u[i+1][j] = dirichlet_boundary(x[i+1],y[j],1);
            }
        }
}

// function to get stencil of matrix of Dxx and Dyy //
void stencil_2d ( double**& arr, int row , int col)
{
for (int i =0;i<row;i++)
    for(int j=0;j<col;j++){
        if(i==j)
            arr[i][j] = -2.0 ;
        else if ( fabs(i-j) == 1)
            arr[i][j] = 1.0 ;
        else
            arr[i][j] = 0.0;
    }
}

//function to print the outstatement
void print_to_output(ofstream& outfile,double**& u, double**& uExact,double*& x, double*& y, double Nx, double Ny,int time_step_no , double time)
{
    double Maximum_error = 0.0,L2_error = 0.0;
   
    outfile<< " VARIABLES = x, y, u_computed" <<endl;
    outfile << "Zone T = psi I = " << Ny+1 << " J = " << Nx+1 << endl ;
    outfile << " StrandID =1, SolutionTime = "<<time_step_no <<endl;
    for(int i = 0 ; i <= Nx; i++)
        for(int j = 0 ; j <= Ny ; j++)
        {
            if( fabs(u[i][j] - uExact[i][j])  > Maximum_error)
                Maximum_error = fabs(u[i][j] - uExact[i][j]) ;
            
            L2_error = L2_error + ( ((u[i][j] - uExact[i][j])*(u[i][j] - uExact[i][j])) / ((Nx+1)*(Ny+1)) );
            outfile <<x[i] <<"\t"<<y[j] <<"\t"<< u[i][j] <<endl;
        }
    L2_error = sqrt(L2_error);
    outfile<<endl<<endl;
    cout <<"Time : "<<time<<"   L2_error : " << L2_error <<"     L-inf error : " <<Maximum_error << endl;
}


void matrix_multiply_u(double**& ans, double** mat1, double** mat2,int row1, int col2, int col1)
{
     for(int i = 1 ; i < row1 ; i++)
        for (int j = 1; j < col2 ; j++)
        {
            ans[i][j] = 0.0;
            for (int k = 1; k < col1 ; k++)
                ans[i][j] += (mat1[i-1][k-1] * mat2[k-1][j-1]);
        }
}
int main ()
{
// Set the precision flags for cout
    cout.flags( ios::dec | ios::scientific );
    cout.precision(5);
    //Declarations
    int Nx,Ny,i,j,k,vec_size;
    double hx,hy;
    // Set Initial Values 
    double x0 = 0.0;double xn  = 1.0;double y0 = 0.0;double yn = 1.0;
    double time = 0.0;
    double time_step = 0.01;
    int time_step_no = 1;
    
   // Get value from user for the Input 
    cout << "Enter N Value : ";
    cin >> Nx;
//     cout << "Enter time stp Value : ";
//     cin >> time_step;
    
    Ny = Nx;hx = hy;
    hx = (xn-x0)/double(Nx);
    hy = (yn-y0)/double(Ny);
    
    // Set the Alpha value to be 2/(del(t) * b)
    double alpha = (2.0 ) / (time_step * B1);
    
    // Declarations of variables
    double *x ,*y, *Eig_val, *norm,error_norm ;
    double **Eig_vec,**uExact,**u,**Eig_vec_trans,**temp1, **force , **F_Tilde , **temp2, **U_Tilde,**temp3;
    double **force_new , **u_new , **stencil, **RHS,**DxxU, **UDyy, **u_old;
    x = new double[Nx+1]; 
    y = new double[Ny+1];
    
    
    // Set Grid points ( UNiform )
    for ( int i = 0; i <= Nx ; i++)  x[i] = i/double(Nx);   
    for ( int i = 0; i <= Ny ; i++)  y[i] = i/double(Ny);
       
    
    // Memory Allocations for the arrays
    memset_2d(Eig_vec,(Nx-1),(Ny-1) ); memset_2d(Eig_vec_trans,(Nx-1),(Ny-1) ); 
    memset_2d(u,Nx+1,Ny+1);memset_2d(uExact,Nx+1,Ny+1);
    memset_2d(force,Nx-1,Ny-1); memset_2d(F_Tilde,Nx-1,Ny-1);memset_2d(U_Tilde,Nx-1,Ny-1);
    memset_2d(temp1,Nx-1,Ny-1); memset_2d(temp2,Nx-1,Ny-1);memset_2d(temp3,Nx-1,Ny-1);
    
    memset_2d(u_new,Nx+1,Ny+1);
    memset_2d(force_new,Nx-1,Ny-1) ;  memset_2d(stencil,Nx-1,Ny-1);memset_2d(RHS,Nx-1,Ny-1);
    memset_2d(DxxU,Nx-1,Ny-1);memset_2d(UDyy,Nx-1,Ny-1);memset_2d(u_old,Nx+1,Ny+1);
    
    Eig_val = new double[Nx-1];
    norm = new double[Nx-1];
    
    
    // Exact Solution for the system
    for(i = 0 ; i <= Nx ; i++)
        for(j = 0 ; j <= Ny ; j++){
            u[i][j] = 0.0;
            uExact[i][j] = sin(2.0*PI*x[i])*sin(2.0*PI*y[j])*sin(2.0*PI*time) / (4.0 * PI * PI);
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
        }
    }
    
    //Calculate the force vector 
    force_calculator(force,x,y,Nx,Ny,hx,hy,u,time);

    
    //Normalising Eigen Vectors 
    for( i = 1 ; i<Nx ; i++)
        for(j = 1;j<Ny;j++)
            Eig_vec[i-1][j-1] = Eig_vec[i-1][j-1]/(sqrt(norm[i-1]));
    
  
        
    // Eigen Vector Inverse
    for( i = 0 ; i<Nx-1 ; i++)
        for(j = 0;j<Ny-1;j++)
           Eig_vec_trans[i][j] = Eig_vec[j][i]; 
    
    ofstream outfile("Output256.dat", ios::out) ;
    outfile << "TITLE = Temp_distribution" <<endl;
    if(!outfile) {cerr<< "Error: Output file couldnot be opened.\n";}
    outfile.flags( ios::dec | ios::scientific);
    outfile.precision(4);
    
    print_to_output(outfile,u,uExact,x,y,Nx,Ny,time_step_no,time);
    //print_matrix(force,Nx-1,Ny-1);
    //get the stencil of the 2d dirichlet boundary matrix
    stencil_2d(stencil,Nx-1,Ny-1);
    #ifdef print
    cout<<"stencil"<<endl;print_matrix(stencil,Nx-1,Ny-1);
    #endif
    // assiging error to prevent the loop from terminating at the start of the operation
    error_norm = 1000.0;
    ////////////////////////////// Loop for time from 2nd time step to last time step. ////////////////////////////////////////
    while (time_step_no < 500 && fabs(error_norm) > 1.0e-13)
    {

    time = time + time_step; 

    // Calculate force at time step n+1
    force_calculator(force_new,x,y,Nx,Ny,hx,hy,u_new,time);
    
    #ifdef print
    cout<<"Force:"<<endl;print_matrix(force_new,Nx-1,Ny-1);
    #endif
    //Calculate Dxx * u
    for(i = 1 ; i < Nx ; i++)
        for (j = 1; j < Ny ; j++){
            DxxU[i-1][j-1] = 0.0;
            for ( k = 1; k < Nx ; k++)
                DxxU[i-1][j-1] += (stencil[i-1][k-1] * u[k][j]);
        }

    #ifdef print 
    cout<<"DxxU"<<endl;print_matrix(DxxU,Nx-1,Ny-1); 
    #endif
    
        // Calculate u*Dyy
    for(i = 1 ; i < Nx ; i++)
        for (j = 1; j < Ny ; j++){
            UDyy[i-1][j-1] = 0.0;
            for ( k = 1; k < Nx ; k++)
                UDyy[i-1][j-1] += (u[i][k] * stencil[k-1][j-1]);
        }
    #ifdef print
    cout<<"UDyy"<<endl;print_matrix(UDyy,Nx-1,Ny-1);
    #endif

    //Calculate the RHS Value
    for(int i=1;i<Nx;i++)
        for(int j=1;j<Ny;j++)
            RHS[i-1][j-1] = 0.0- (1.0*(force[i-1][j-1] +  force_new[i-1][j-1] + DxxU[i-1][j-1] + UDyy[i-1][j-1] + (alpha * hx*hx*u[i][j])));
    #ifdef print
    cout<<"RHS"<<endl;print_matrix(RHS,Nx-1,Ny-1);
    #endif
    
    // Calculate F_tilde =  Eig_vec' * RHS * Eig_vec
    // ( Note : matrix Multiply ( answer, mat1 , mat2 , row(mat1) , col(mat2) , col(mat1)))
    matrix_multiply(temp1,Eig_vec_trans,RHS,Nx-1,Ny-1,Nx-1);
    matrix_multiply(F_Tilde,temp1,Eig_vec,Nx-1, Ny-1,Nx-1);
    #ifdef print
    cout<<"F_Tilde"<<endl;print_matrix(F_Tilde,Nx-1,Ny-1);
    #endif

    // Calulate U-Tilde --> F_Tilde / (sum of eig_values + alpha)
    for(i=0;i<Nx-1;i++)
        for(j=0;j<Ny-1;j++)
            U_Tilde[i][j]  = F_Tilde[i][j] / (Eig_val[i] + Eig_val[j] - (alpha*hx*hx)) ;
    
    #ifdef print
    cout<<"U_Tilde"<<endl;print_matrix(U_Tilde,Nx-1,Ny-1);
    #endif
    // Calulate U --> eig_vec * u-tilde * eig_vec_trans
    matrix_multiply(temp2,Eig_vec,U_Tilde,Nx-1,Ny-1,Nx-1);

    for(i = 1 ; i < Nx ; i++)
        for (j = 1; j < Ny ; j++){
            u_new[i][j] = 0.0;
            for ( k = 1; k < Nx ; k++)
                u_new[i][j] += (temp2[i-1][k-1] * Eig_vec_trans[k-1][j-1]);   
        }
         
    // Mean Square Error Calculation      
    
    #ifdef print
    cout<<"U_new : "<<endl;
    for( i = 1 ; i < Nx ; i++){
       for( j=1 ; j< Ny ; j++)
        cout<<u_new[i][j]<<"\t";    
         cout<<endl;      
    }
    #endif
        // assign U new to U , and Fnew to F
        // Error Calculation
    error_norm = 0.0;
    for(i = 1 ; i < Nx ; i++)
        for (j = 1; j < Ny ; j++){
         error_norm +=  (( u_new[i][j] - u[i][j] ) * ( u_new[i][j] - u[i][j] ) ) / (pow(Nx-1,2));
          u[i][j] = u_new[i][j];
          force[i-1][j-1] = force_new[i-1][j-1];
          //  cout<<"u_new : "<<u_new[i][j]<<"   u : "<<u[i][j]<<endl;
            
        }
        //error_norm =  error_norm);
        cout<<"Error Norm : "<<error_norm<<"   at iter : "<<time_step_no<<endl;
    
      for(i = 0 ; i <= Nx ; i++)
        for(j = 0 ; j <= Ny ; j++){
            uExact[i][j] = sin(2.0*PI*x[i])*sin(2.0*PI*y[j])*sin(2.0*PI*time) / (16.0 * PI * PI);
	    }
    // Outstream parameters
    time_step_no++;
    print_to_output(outfile,u_new,uExact,x,y,Nx,Ny,time_step_no,time);
    
    // Stopping Criteria - Mean Error 

    }
    ////////////////////////////////////////////// Loop Close ///////////////////////////////////////////////////////////
    outfile.close() ; 
    cout<< " total iterations : "<<time_step_no<<endl;
    cout<< " last error : "<<error_norm<<endl;
    // Free up the heap memory 
   delete Eig_val;delete Eig_vec;delete Eig_vec_trans;delete x; delete y;delete temp1;delete temp2;
   delete u;delete uExact;
}
