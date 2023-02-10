// Implementing L2 norm, high res vs low res
// compile with:
// g++ -O3 TEM_2d.cpp -o run_TEM

#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<numeric>
#include<cstdlib>
#include<time.h>
#include<chrono>

typedef std::vector< double > Real1D;
typedef std::vector< Real1D > Real2D;
typedef std::vector< Real2D > Real3D;

Real1D fill_coordinate( double wmin, double wmax, int nw ){
  Real1D w( nw );  
  double dw = ( wmax - wmin ) / ( nw - 1 );
  for ( int i = 0; i < nw; ++i ){
    w[i] = wmin + i * dw;
  }
  return w;
}

double initial_data_function( double x, double y, double z ){
  // set to a simple Gaussian: 
  //double r = sqrt( x*x + y*y + z*z );
  //double sigma = 0.2;
  //return exp( -0.5 * r * r / ( sigma * sigma ) ); 

  int intitial_form_type = 1;

  // Update: set to cos^4 function - has finite support
  double r = sqrt( x*x + y*y + z*z );
  double lambda = 1.2;

  if(intitial_form_type == 0){
    if ( r < ( lambda / 4. ) ){ // nonzero up to pi/2
	    return  cos( 2. * M_PI * r / lambda ); 
    } // else:
  }

  else if (intitial_form_type == 1){
    if ( r < ( lambda / 4. ) ){ // nonzero up to pi/2
      return pow( cos( 2. * M_PI * r / lambda ), 4 ); 
    } // else:
  }
  return 0.; 
}

Real2D create_matrix(double nx, double ny){
	Real1D temp1(  ny , 0.0  );
	Real2D matrix( nx, temp1 );  
	return matrix;
}

Real3D create_3D_matrix(int nx, int ny, int nz){
	Real1D temp1( nz , 0.0 );
	Real2D temp2( ny, temp1);  
	Real3D matrix(nx, temp2);
	return matrix;
}

void fill_initial_data( const Real1D& x, const Real1D& y, const Real1D& z, Real3D& psi ){
  for ( int i = 0; i < x.size(); ++i ){
    for ( int j = 0; j < y.size(); ++j ){
      for ( int k = 0; k < z.size(); ++k ){
        psi[ i ][ j ][ k ] = initial_data_function( x[i], y[j], z[k] );
      }
    }
  }
}

double gaussian_TM( const int nx, const int ny, const double t, const double t_0, const double thalf){
  //double data = exp(- pow(((t - t_0) / thalf),2));
  double data = -1.0;
  return data;
}

Real2D reflecting_matrix(const Real1D& x, const Real1D& y){
  Real2D matrix = create_matrix(x.size(), y.size());
  for ( int i = 1; i < x.size()-1; ++i ){
    for ( int j = 1; j < y.size()-1; ++j ){
      matrix[i][j] = 2.0;
    }
  }
  return matrix; 
}

Real2D plane_wave(const int nx, const int ny, const double t, Real2D matrix) {
  //Real2D matrix = create_matrix(x.size(), y.size());
  
  const double omega = 1e12;
  for ( int j = 1; j < ny-1; ++j ){
    matrix[1][j] = 10000*cos(omega*t);
  }

  return matrix; 
}

void print_Real2D(const Real2D& matrix){
	std::cout<< "----------------> x " << std::endl;
	for(int k = 0; k<matrix.size();k++){
		if(k<matrix.size()-1){std::cout<< "|    ";} else{std::cout << "y    ";}
		for(int i=0; i<matrix[k].size(); i++){
			std::cout << matrix[k][i] << "  ";
			} 
			std::cout<< "" <<std::endl;
	}	
	std::cout<<"y"<< std::endl;
}

Real2D assign_eps_r(const int nx, const int ny ){
  Real2D matrix = create_matrix(nx, ny);
  const int xmin = nx/4; 
  const int xmax = (3*nx)/4;
  const int ymin = ny/4;
  const int ymax = (3*ny)/4;
  for ( int i = 0; i < nx; ++i ){
    for ( int j = 0; j < ny; ++j ){
      if(i > xmin && i < xmax && j > ymin && j < ymax){
        matrix[i][j] = 1.0;
      }
      else{matrix[i][j] = 1.0;}
    }
  }
  return matrix;
}

// Floating dielectic assignment
Real2D assign_eps_r_floating(const Real1D& x, const Real1D& y, const double e_r){
	int nx = x.size();
	int ny = y.size();
  Real2D matrix = create_matrix(nx, ny);
  const double xmax  = x[x.size()-1]/2; 
  const double xmin =  x[0]/2;
  const double ymax  = y[y.size()-1]/2;
  const double ymin =   y[0]/2;
  
  	//std::cout << "xmin = " << xmin << ", xmax = " << xmax << ", ymin = " << ymin << ", ymax = " << ymax <<  std::endl;

  for ( int i = 0; i < x.size(); ++i ){
    for ( int j = 0; j < y.size(); ++j ){
      if(x[i] > xmin && x[i] < xmax && y[j] > ymin && y[j] < ymax){
        matrix[i][j] = e_r;
      }
      else{matrix[i][j] = 1.0;}
    }
  }
  return matrix;
}

double calc_L2_norm(const Real2D& high_res, const Real2D& current, const int high_mult, const int mult, const int nx, const Real1D& x, const Real1D& y){
	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double L2_norm = 0.0;
	int ny = nx;

	for(int i = 0; i < nx; ++i){
		for(int j = 0; j<ny; ++j){
			double high_res_val = high_res[ i*high_mult ][ j*high_mult ];
			double current_val = current[ i*mult ][ j*mult ];
			double diff = high_res_val - current_val;
			L2_norm += diff*diff * dx * dy;
		}
	}
	
	return sqrt(L2_norm);
	
}

void print_symmetry_check(const int nx, const int ny, const Real2D& Ez){
	std::cout << "-------------------------------------------------------" << std::endl;
	std::cout << "Check: Horizontal  Ez[1][ny/2] = Ez[nx-2][ny/2] = Ez["<<1<<"][" << ny/2 <<  "] = Ez[" << nx-2 << "][" << ny/2<< "] = " <<  Ez[1][ny/2]  << " =? " <<  Ez[nx-2][ny/2]  << " ----> Difference = " << Ez[1][ny/2]  - Ez[nx-2][ny/2]  << std::endl;
	std::cout << "Check: Vertical    Ez[nx/2][1] = Ez[nx/2][ny-2] = Ez["<<nx/2<<"][" << 1 <<  "] = Ez[" << nx/2 << "][" << ny-2<< "] = " <<  Ez[nx/2][1] << " =? " <<  Ez[nx/2][ny-2]   << " ----> Difference = " << Ez[nx/2][1] - Ez[nx/2][ny-2]  << std::endl;
	std::cout<< std::endl;
	std::cout << "Check: Diagonal    Ez[1][1]    = Ez[nx-2][ny-2]   = Ez["<<1<<"][" << 1 <<  "] = Ez[" << nx-2 << "][" << ny-2<< "] = " <<  Ez[1][1] << " =? " <<  Ez[nx-2][ny-2]   << " ----> Difference = " << Ez[1][1] - Ez[nx-2][ny-2]  << std::endl;
	std::cout << "Check: Diagonal    Ez[1][ny-2] = Ez[nx-2][1]      = Ez["<<1<<"][" << ny-2 <<  "] = Ez[" << nx-2 << "][" << 1<< "] = " <<  Ez[1][ny-2] << " =? " <<  Ez[nx-2][1]   << " ----> Difference = " << Ez[1][ny-2] - Ez[nx-2][1]  << std::endl;
}

int main(){
  // Constants
  double c   = 299792458;
  double mu  = 4*M_PI*(0.0000001);
  double eps = 8.85418 * pow(10,-12);
  
  // -----------------------------
  // initialize values
  std::vector <int> nx_array{ 3201, 1601, 801, 401, 201, 101, 51 };
  double e_r = 5.0;
  // 1 = Iz, 
  // 2 = Iz*dx
  // 3 = Iz*dx*dx
  int Iz_mode = 3;
  // -----------------------------  
  

  
  
  
  std::vector<Real2D> Ez_final;
  std::vector<Real1D> x_store;
  std::vector<Real1D> y_store;
  std::vector<int> current_multiple;
  std::vector<double> dt_vals;
  std::fstream L2_log("L2_log.dat", std::ios::out);
  
	
	for(int n = 0 ; n < nx_array.size(); ++n ){
		
	// Grid 
    // Set the domain and resolution:
	int nx = nx_array[n];
	std::cout << "nx = " << nx << std::endl;
    double xmin =  -0.1;
    double xmax =  0.1;
    double dx = ( xmax - xmin ) / ( nx - 1 );
    Real1D x = fill_coordinate( xmin, xmax, nx ); 
	x_store.push_back(x);	

    // go ahead and set y same as x:
    double ymin = xmin;
    double ymax = xmax;
    int ny = nx;
    double dy = dx;
    Real1D y = fill_coordinate( ymin, ymax, ny );
	y_store.push_back(y);	

	int t_mult = (nx_array[0]-1)/(nx_array[n]-1);
	int mult    = (nx_array[n]-1)/(nx_array[nx_array.size()-1]-1);
	current_multiple.push_back(mult);
	std::cout << "current multiple =  " << mult << std::endl;
	
    // Courant condition:
    double dt = 1.0 /(c *  sqrt( 1./(dx*dx) + 1./(dy*dy)));
    double t = -dt;
	//int nx_mult = (nx_array[0]-1)/(nx_array[nx_array.size()-1]-1);
    double t_final = dt*100*current_multiple[0]/t_mult;
    int num_steps = int( t_final / dt );
	std::cout << "dt = " << dt << ", t_final = " << t_final << ", num_steps = " << num_steps << std::endl;  
    int report_total = 100;
    int report_period = num_steps / report_total;
    int step = -1;
    
    // Variables for periodic saving of field data
    //int slice_count = 20;
    int num_slices = 30;
    int slice_period = int ((num_steps+1)/num_slices);
    int report_count = 1;
    std::fstream log2D_function( "log_2d.dat", std::ios::out );
    std::fstream log1D_function( "log_1d.dat", std::ios::out );
    
    // Ez field data
    double thalf = 20;
    double t_0 = thalf * 3;
     
    // initialize fields
    Real2D Hy = create_matrix(nx, ny);
    Real2D Hx = create_matrix(nx, ny);
    Real2D Ez = create_matrix(nx, ny);
    Real2D reflecting = reflecting_matrix(x, y);
    
    // Initalize dielectics
    Real2D eps_r = assign_eps_r_floating(x, y, e_r);
    
    
	  
    while ( t < t_final ){
    	t += dt;
    	step += 1;
	    //std::cout << "Percent done = " << (t/t_final)*100 << " %" <<std::endl;
	    
	    // Set E field, t = t
	    int simulation_type = 3;
		if(Iz_mode>0){
			if(Iz_mode == 1){simulation_type = 3;}
			if(Iz_mode == 2){simulation_type = 4;}
			if(Iz_mode == 3){simulation_type = 5;}
		} 
		
	    // E field excitation
	    // Type 0 - Gaussian Pulse source at center
	    // Type 1 - Plane wave source excitation from the side
	    // Type 2 - Single pulse plane wave source
	    // Type 3 - Single Pulse at center
	    if(simulation_type == 0){Ez[nx/2][ny/2] = gaussian_TM(nx, ny, t, t_0, thalf);}
	    else if (simulation_type == 1){Ez = plane_wave(nx, ny, t, Ez);}
	    else if (simulation_type == 2){if(step == 0){Ez = plane_wave(nx, ny, t, Ez);}}
		else if (simulation_type == 3){if(step == 0){Ez[nx/2][ny/2] = gaussian_TM(nx, ny, t, t_0, thalf);}}
		else if (simulation_type == 4){if(step == 0){Ez[nx/2][ny/2] = gaussian_TM(nx, ny, t, t_0, thalf)*dx;}}
	    else if (simulation_type == 5){if(step == 0){Ez[nx/2][ny/2] = gaussian_TM(nx, ny, t, t_0, thalf)*dx*dx;}}
	    
	    
	    // Solve for H fields, t = t + dt/2
	    for ( int i = 0; i < nx; ++i ){
	    	for ( int j = 0; j < ny; ++j ){

	    	  // Account for edge cases
	      	if(i == nx-1 || j == ny-1){
	      	  // Do not update at all for corner
	      	  if(i == nx-1 && j == ny-1){
	      	    Hx[i][j] = Hx[i][j];
	      	    Hy[i][j] = Hy[i][j];
	      	  }
	      	  // Final column
	      	  else if(j == ny-1){
	        	  // Do not update Hx
	        	  Hx[i][j] = Hx[i][j];
	        	  // Update Hy 
	        	  Hy[i][j] = Hy[i][j] + (dt/(mu*dx))      * (Ez[i+1][j] - Ez[i][j]);
	        	}
	        	//Final row
	        	else if(i == nx-1){
	        	  // Update Hx
	        	  Hx[i][j] = Hx[i][j] + (dt/(mu*dy))      * (Ez[i][j] - Ez[i][j+1]); 
	        	  // Do not update Hy
	        	  Hy[i][j] = Hy[i][j];
	        	}
	        }
	        // For all other grid points, update normally
	      	else{	          
	        	Hx[i][j] = Hx[i][j] + (dt/(mu*dy))      * (Ez[i][j] - Ez[i][j+1]); 
			      Hy[i][j] = Hy[i][j] + (dt/(mu*dx))      * (Ez[i+1][j] - Ez[i][j]);
			    }
		    }
	    }
	    	    
	    // Solve for E fields, t = t + dt
	    for(int i = 1; i < nx-1; ++i){
		    for(int j = 1; j < ny-1; ++j){
			  Ez[i][j] = Ez[i][j] + (dt/(eps*eps_r[i][j])) * ((Hy[i][j] - Hy[i-1][j])/dx + ( Hx[i][j-1] - Hx[i][j] )/dy );
		    }
	    }
	    
	    	    
	    // Log data as specified
	    if (step % slice_period == 0 ){
			//print_symmetry_check(nx, ny, Ez);
           for ( int i = 0; i < nx; ++i ){
              for ( int j = 0; j < ny; ++j ){
          	     log2D_function << x[i] << " " << y[j] << " " << Ez[i][j] << std::endl;
          	     if( j == int (ny/2) ){
          	        log1D_function << x[i] << " " << " " << Ez[i][j] << std::endl;
          	     } 	     
          	  }
          	  log2D_function << std::endl;
           } 
           log2D_function << std::endl;
           log2D_function << std::endl;   
           log1D_function << std::endl;
           log1D_function << std::endl;    
           //slice_count += 1;
       }
    }
	
	Ez_final.push_back(Ez);
	
	}
	
	for(int i =0; i < nx_array.size(); ++i){
		//print_Real2D(matrix);
		Real2D Ez = Ez_final[i];
		int nx = nx_array[i];
		int ny = nx;
		int low_res_nx = nx_array[nx_array.size()-1];
		int mult = current_multiple[i];
		int high_res_mult = current_multiple[0];
		Real1D x = x_store[i];
		Real1D y = y_store[i];
		
		//std::cout << "nx = " << nx << ", current multiple = " << mult << std::endl;
		if(i>0){
			double L2_norm = calc_L2_norm(Ez_final[0], Ez, high_res_mult, mult, low_res_nx, x, y);
			L2_log << nx_array[i] << " " << L2_norm << std::endl;
			
		}

	}
	
  
return 0;

}
