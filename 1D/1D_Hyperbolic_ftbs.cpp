#include<iostream>
#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////   Standing Wave Equation solver in 1D    < Periodic Bound Cond >           //////////////////
//////////              u(x,0) =  g(x), u(x+1,t) = u(x,t )                            ///////////////////
//////////              Ut  + C.Ux = 0                                                //////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


// Set the value of pi
// double const PI = 4.0*atan(1.0);
double const PI = 3.141592654;

// Set value of B
double const B1 = 1.0;


// -------------------------------- Start User Input Parameters --------------------------// 
double const x0 = 0.0;
double const xn  = 1.0;
double const time_step = 0.005;
int const max_iteration = 1000;
double const C = 1.0;

// uncomment the below line to Suppress the output parameters
#define print 1;

// Define any one of the Below Parameters 
#define CTCS 1
// #define  FTBS 1


// ------------------------------ End of user input parameters ---------------------------------//

// Function to set the memory for array and initialising the values to Zero

void memset_1d( double*& arr, int row)
{
    arr = new double[row] ;
    for (int i = 0 ; i < row ; i ++ )
        arr[i] = 0.0;
}

void mat_vec_mul(double*& ans, double** mat , double* vec, int Nx )
{
    
}

void validation()
{
    int a  = 0;
    
    #ifdef CTCS
    a += 1; cout <<" CTCS Scheme " <<endl;
    #endif
    
    #ifdef FTBS
    a += 1;cout <<" FTBS Scheme " <<endl;
    #endif
    
    if( a >1){throw std::runtime_error("Error Message : Multiple Methods Selected");}
        
}
//Function to print Array ( 1 d )
void print_array(double* arr1,int row)
{
    for(int i=0;i<row;i++)
            cout<<arr1[i]<<"\t";
    cout<<endl;
}

void print_to_output(ofstream& outfile,double* u, double* uExact,double* x, double Nx,int time_step_no , double time, double cfl)
{
    double Maximum_error = 0.0,L2_error = 0.0;
    outfile<< " Title  =  \" FTBS -  CFL : " <<cfl<<"  Time : "<<time <<"\" "<<endl;
    outfile<< " VARIABLES = x, u_computed" <<endl;
    outfile << "Zone T = uComputed ,  I = " << Nx+1 <<" , C = RED"<< endl ;
    outfile << " StrandID = 1, SolutionTime = "<<time_step_no <<endl;
    //outfile<<"Time : "<<time<<endl;
//  
    for(int i = 0 ; i <= Nx; i++)
    {
        if( fabs(u[i] - uExact[i])  > Maximum_error)
            Maximum_error = fabs(u[i] - uExact[i]) ;
        
        L2_error = L2_error + ( ((u[i] - uExact[i])*(u[i] - uExact[i])) / (Nx+1) );
        outfile <<x[i] <<"\t"<< u[i] <<endl;
    }
     L2_error = sqrt(L2_error);
    outfile<< " Title  =  \" FTBS -  CFL : " <<cfl<<"  Time : "<<time <<"\" "<<endl;
    outfile<< " VARIABLES = x, uExact" <<endl;
    outfile << "Zone T = uExact , I = " << Nx+1 <<" , C = BLUE"<< endl ;
    outfile << " StrandID = 2,  SolutionTime = "<<time_step_no <<endl;
    for(int i = 0 ; i <= Nx; i++)
    {
           outfile <<x[i] <<"\t"<<uExact[i] <<endl;
    }
    //outfile<<" TEXT X = 0.2, Y = 1.1,  T = \"FTBS  -  CFL : " <<cfl<<"  Time : "<<time <<"\" " <<endl;
   
    outfile<<endl<<endl;
    cout <<"Time : "<<time<<"   L2_error : " << L2_error <<"     L-inf error : " <<Maximum_error << endl;
}



// Main Function 

int main ()
{
    // Set the precision flags for cout
    cout.flags( ios::dec | ios::scientific );
    cout.precision(5);
    
    //Declarations
    int Nx,i,j,k,vec_size , count;
    double error,l2_error, iteration_error , h , time;
    int time_step_no = 0;
    time = 0.0;
    cout << "Enter Nx Value : ";
    cin >> Nx;
    
    // Validation for declaration of multiple Variables
    validation();
    
    // Calculate h and Cfl Values
    h = (xn-x0)/double(Nx);
    double cfl = (C*time_step)/(h);
    
    // Declarations of variables 
    double *x, *error_norm, *u , *u_new, *force , *uExact  ;
    
    // Allocation of memory for arrays
    memset_1d(x,Nx+1);memset_1d(error_norm,Nx+1);memset_1d(u,Nx+1);memset_1d(u_new,Nx+1);
    memset_1d(uExact,Nx+1);
    
    #ifdef CTCS
    double *u_old;
    memset_1d(u_old,Nx+1);
    #endif
    
    // Set Grid points ( UNiform ) in X
    for ( int i = 0; i < Nx+1 ; i++)  
        x[i] = i/double(Nx);  
    
    error = 0.0;
    // Set Initial value of U at time step = 0.0 and U_exact 
    for ( i = 0 ; i < Nx+1 ; i++ ){
        u[i] = cos(4*PI*x[i]);
        uExact[i] = cos(4*PI*(x[i] - (C*time)));
        error += pow((u_new[i] - uExact[i]),2) / (Nx+1);
    }
    
    // print variables to output file
    ofstream outfile("Output001.dat", ios::out) ;
    if(!outfile) {cerr<< "Error: Output file couldnot be opened.\n";}
    outfile.flags( ios::dec | ios::scientific);
    outfile.precision(4);
    
    print_to_output(outfile,u,uExact,x,Nx,time_step_no,time,cfl);
    
    count = 0;
    #ifdef CTCS   
    
    // Increase time step by 1
    time =  time + time_step;
    for ( i = 0 ; i < Nx+1 ; i++ ){
        u_old[i] = u[i]; 
        u[i] = 0.0;
        u[i] = cos(4*PI*(x[i] - (C*time)));
        uExact[i] = cos(4*PI*(x[i] - (C*time)));
        error += pow((u_new[i] - uExact[i]),2) / (Nx+1);
    }
    print_to_output(outfile,u,uExact,x,Nx,time_step_no,time,cfl);
    count = count + 1;
    #endif
    
    
    // Initialise Count for Loop Iterations for entry Condition
    
    iteration_error = 10e4;
  
    
    // --------------------- Loop Starts here ------------------- //
    while ( count < max_iteration &&  fabs(iteration_error) > 1.0e-8 )
    {
        time =  time + time_step;

        error = 0.0;
        // Calculate u_New
        
        #ifdef FTBS
        for ( i = 0 ; i< Nx+1 ; i++){  
            uExact[i] = cos(4*PI*(x[i] - (C*time)));
            if(i != 0)
                u_new[i] = u[i] + cfl*(u[i-1] - u[i]);
            else  // Periodic Boundary Condition
               u_new[i] = u[i] + cfl*(u[Nx] - u[i]); 
            iteration_error +=  pow( (u_new[i] - u[i] ),2 ) / (Nx-1);
        }
        #endif
        
        #ifdef CTCS      
        for( i = 0 ; i< Nx+1 ; i++ ){
            uExact[i] = cos( 4 * PI * (x[i] - (C*time)));
            if(i == Nx){
                u_new[i] = u_old[i] + cfl*(u[i-1] - u[0]); }
            else if(i == 0){
               u_new[i] = u_old[i] + cfl*(u[Nx] - u[i+1]); 
            }
                
            else
               u_new[i] = u_old[i] + cfl*(u[i-1] - u[i+1]); 
        } 
        #endif
        
        
        // Print the Values to the Output file
        print_to_output(outfile,u_new,uExact,x,Nx,time_step_no,time,cfl);
        
        // Transfer the Previous Iteration variables to the new ones. 
        for(i=0;i< Nx+1;i++)
        {
        #ifdef CTCS
            u_old[i] = u[i];
        #endif
           u[i] = u_new[i];
        }
       
        
        // Increase the count for Iteration 
        count += 1;
        time_step_no += 1;
        
    }
    // Close the output
    outfile<< " TEXT X=0.3, Y=1.1, T = \" FTBS - CFL : "<< cfl << " at iter :  (&SolutionTime) (Del(t) - "<<time_step<<") \" " <<endl;
    outfile.close() ; 
    cout<< " total iterations : "<<time_step_no  <<endl;
    
    // Delete the used up variables
    //delete u; delete u_new ; delete uExact ; delete x ; delete error_norm;
    
    
    return 0;
    
    
    

}
