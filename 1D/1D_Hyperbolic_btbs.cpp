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


// -------------------------------- Start User Input Parameters --------------------------// 
double const x0 = 0.0;
double const xn  = 1.0;
double const time_step = 0.005 ;   
int const max_iteration = 1000;
double const C = 1.0;
string method;

// uncomment the below line to Suppress the output parameters
#define print 1;

// Define any one of the Below Parameters 
//#define CTCS 1
 #define  FTBS 1
//#define BTBS 1
//#define CRANKNICHOLSON 1


// ------------------------------ End of user input parameters ---------------------------------//

// Function to set the memory for array and initialising the values to Zero

void memset_1d( double*& arr, int row)
{
    arr = new double[row] ;
    for (int i = 0 ; i < row ; i ++ )
        arr[i] = 0.0;
}

void vector_dot(double*& ans, double* vec1 , double* vec2,int Nx)
{
    for(int i = 0 ; i< Nx+1 ; i++)
        ans[i] = vec1[i]*vec2[i];
}

void vector_inner_dot(double& ans, double* vec1 , double* vec2,int Nx)
{
    for(int i = 0 ; i< Nx+1 ; i++)
        ans += vec1[i]*vec2[i];
}
    

void thomasAlgorithm(double*& u , double* a, double* b , double* c , double* d, int Nx)
 {
    double *p,*q ; double temp = 0.0;
    memset_1d(p,Nx+1);memset_1d(q,Nx+1);
    
    p[0] = c[0]/b[0];q[0] = d[0]/b[0];
    
    for (int i = 1 ; i < Nx+1 ; i++){
        temp = b[i] - (a[i]  *p[i-1]);
        p[i] = c[i] / temp;
        q[i] = ( d[i] - (a[i]*q[i-1]) ) / temp;
    }
    // Set value f U[N-1]
    u[Nx] = q[Nx];
    //backward Substitution
    for ( int i = Nx-1 ; i >= 0 ; i--)
        u[i] = -(p[i]*u[i+1]) + q[i]; 
    
    delete p,delete q;
}

void sherman(double*& ans, double* yy, double* zz,double* n,int Nx)
{
    double *rhs;
    double nty = 0.0, ntz = 0.0; 
    memset_1d(rhs,Nx+1); 
    
    vector_inner_dot(nty,n,yy,Nx);
    vector_inner_dot(ntz,n,zz,Nx);
    ntz += 1.0;
    for(int i = 0 ; i<Nx+1 ; i++){
        rhs[i] = 0.0;
        rhs[i] = (nty/ntz) * zz[i];
        ans[i] = yy[i] - rhs[i];
    }
    delete rhs;
}




void validation()
{
    int a  = 0;
    #ifdef CTCS
    a += 1; cout <<" CTCS Scheme " <<endl; method = "CTCS";
    #endif
    
    #ifdef FTBS
    a += 1;cout <<" FTBS Scheme " <<endl; method = "FTBS";
    #endif
    
    #ifdef BTBS
    a += 1;cout <<" BTBS Scheme " <<endl; method = "BTBS";
    #endif
    
    #ifdef CRANKNICHOLSON
    a += 1;cout <<" CRANKNICHOLSON Scheme " <<endl; method = "CRANKNICHOLSON";
    #endif
    
    if( a != 1){throw std::runtime_error("Error Message : Multiple Methods Selected");}
        
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
    outfile<< " Title  =  \" "<< method <<" -  CFL : " <<cfl<<"  Time : "<<time <<"\" "<<endl;
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
    outfile<< " Title  =  \" "<< method <<" -  CFL : " <<cfl<<"  Time : "<<time <<"\" "<<endl;
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
    double error,iteration_error , h , time;
    double *x, *error_norm, *u , *u_new , *uExact  ;
    int time_step_no = 0;
    time = 0.0;
    cout << "Enter Nx Value : ";
    cin >> Nx;
    
    // Validation for declaration of multiple Variables
    validation();
    
    // Calculate h and Cfl Values
    h = (xn-x0)/double(Nx);
    double cfl = (C*time_step)/(h);
    cout<<" cfl : "<<cfl<<endl;
   
    // Allocation of memory for arrays ----- //
    memset_1d(x,Nx+1);memset_1d(error_norm,Nx+1);memset_1d(u,Nx+1);memset_1d(u_new,Nx+1);
    memset_1d(uExact,Nx+1);
    
    #ifdef CTCS
    double *u_old;
    memset_1d(u_old,Nx+1);
    #endif
    
    #if defined BTBS || defined CRANKNICHOLSON
    //for thomas Algorithm
    double *aa, *bb , *cc , *m , *n  , *yy , *zz;
    memset_1d(aa,Nx+1);memset_1d(bb,Nx+1);memset_1d(cc,Nx+1);memset_1d(m,Nx+1);memset_1d(n,Nx+1);
    memset_1d(yy,Nx+1);memset_1d(zz,Nx+1);
    #endif
    
    #ifdef CRANKNICHOLSON
    double *o, *p, *y1_tilde, *y2_tilde , *z1_tilde , *u_tilde_new ;
    memset_1d(o,Nx+1);memset_1d(p,Nx+1);memset_1d(y1_tilde,Nx+1);memset_1d(y2_tilde,Nx+1);
    memset_1d(z1_tilde,Nx+1);memset_1d(u_tilde_new,Nx+1);
    #endif
    
    
    // ---- end of memory Allocation ------//
    
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
    
    string filename = "Output_"  + method + "_" + to_string(cfl) + ".dat";
    // print variables to output file
    ofstream outfile(filename, ios::out) ;
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
    
    #ifdef BTBS
    // Tridiagonal Matrix Creator
    for(i = 0;i<Nx+1;i++){
        aa[i] = -cfl; 
        bb[i] = 1.0+cfl; 
        cc[i] = 0.0;  
    }
    
    aa[0] = 0.0;cc[Nx] = 0.0;
    m[0] = -cfl;n[Nx] = 1.0;
    #endif
    
     #ifdef CRANKNICHOLSON
    // Tridiagonal Matrix Creator
    for(i = 0;i<Nx+1;i++){
        aa[i] = -1.0; 
        bb[i] = 4.0*cfl; 
        cc[i] = 1.0;  
    }
    
    aa[0] = 0.0;cc[Nx] = 0.0;
    m[0] = 1.0; n[Nx] = -1.0;
    o[Nx] = 1.0 ; p[0] = 1.0;
    #endif
    // Initialise Count for Loop Iterations for entry Condition
    iteration_error = 10e4;
    // --------------------- Loop Starts here ------------------- //
    while ( count < max_iteration &&  fabs(iteration_error) > 1.0e-8 )
    {
        time =  time + time_step;
        error = 0.0;
        // Calculate u_New based on any one of the methods

        //------- FTBS --------------//
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
        //------- FTBS(End) --------------//
        
        //------- CTCS --------------//
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
        //------- CTCS(end) --------------//
        
        
        #if defined BTBS || defined CRANKNICHOLSON 
        for( i = 0 ; i< Nx+1 ; i++ )
            uExact[i] = cos( 4 * PI * (x[i] - (C*time)));
        #endif
        
        //------- BTBS --------------//
        #ifdef BTBS
        thomasAlgorithm(yy,aa,bb,cc,u,Nx);
        thomasAlgorithm(zz,aa,bb,cc,m,Nx);
        sherman(u_new,yy,zz,n,Nx);     
        #endif
        //------- BTBS  (end) --------------//
        
        //------- CRANKNICHOLSON --------------//
        #ifdef CRANKNICHOLSON
        // Calculate U_tilde_new
        for ( int i = 0 ; i< Nx+1 ; i++)
        {
            u_tilde_new[i] = 0.0;
            if( i ==0)
                u_tilde_new[i] = (4*cfl*u[i]) - u[i+1] + u[Nx];
            else if ( i == Nx )
                u_tilde_new[i] = (4*cfl*u[i]) - u[0] + u[i-1];
            else
                u_tilde_new[i] = (4*cfl*u[i]) - u[i+1] + u[i-1];
        }

        thomasAlgorithm(y1_tilde,aa,bb,cc,u_tilde_new,Nx);
        //cout<<"y1_tilde : "<<endl;print_array(y1_tilde,Nx+1);
        thomasAlgorithm(y2_tilde,aa,bb,cc,m,Nx);
        //cout<<"y12_tilde : "<<endl;print_array(y2_tilde,Nx+1);
        thomasAlgorithm(z1_tilde, aa,bb,cc, o, Nx);
        //cout<<"z1_tilde : "<<endl;print_array(z1_tilde,Nx+1);

        
        sherman(yy, y1_tilde, z1_tilde , p ,Nx);
        //cout<<"yy : "<<endl;print_array(yy,Nx+1);
        sherman(zz,y2_tilde,z1_tilde,p,Nx);
        //cout<<"zz : "<<endl;print_array(zz,Nx+1);
        sherman(u_new,yy,zz,n,Nx);
        //cout<<"u_new : "<<endl;print_array(u_new,Nx+1);
        #endif
        
        
        
        //------- CRANKNICHOLSON(end) --------------//

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
    outfile<< " TEXT X=0.4, Y=1.0, T = \""<< method <<" -  CFL : "<< cfl << " at iter :  &(SolutionTime) (Del(t) - "<<time_step<<") \" " <<endl;
    outfile.close() ; 
    cout<< " total iterations : "<<time_step_no  <<endl;
    cout<<" CFL : "<<cfl<<endl;
    cout<<"Method : "<<method<<endl;
    
    // Delete the used up variables
    delete u; delete u_new ; delete uExact ; delete x ; delete error_norm;
    
    #ifdef BTBS
    delete aa;delete bb ;delete cc ;delete m ;delete n  ;delete yy ;delete zz;
    #endif
    
    #ifdef CTCS
    delete u_old;
    #endif
    
    
    
    #ifdef CRANKNICHOLSON
    delete o;delete p;delete y1_tilde;delete y2_tilde ;delete z1_tilde;
    #endif
    
    
    return 0;
    
    
    

}
