// Kernel density estimation by Tim Nugent (c) 2014
// Based on Philipp K. Janert's Perl module:
// http://search.cpan.org/~janert/Statistics-KernelEstimation-0.05 
// Multivariate stuff from here:
// http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/spm/spmhtmlnode18.html

#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <iostream>
#include <vector> 
#include <map>
#include "kde.h"
#include <QDebug>

void KDE::calc_bandwidth(){

	for(curr_var = 0; curr_var < data_matrix.size(); curr_var++){
		if(bandwidth_map[curr_var] == -1.0){
			switch(bandwidth_opt_type){
				case 2: 
					optimal_bandwidth();
					break;
				case 3: 
					optimal_bandwidth_safe();
					break;	
				default:
					default_bandwidth();						
					break;
			}
		}
	}
}

void KDE::set_bandwidth_opt_type(int x){

	if(kernel == 1){
    	bandwidth_opt_type = x;
    }	
}

void KDE::set_kernel_type(int x){

    kernel = x;
    if(kernel != 1){
    	extension = 0;
    	bandwidth_opt_type = 1;
    }
}

double KDE::pdf(double x){
	vector<double> tmp;
	tmp.push_back(x);
	return(pdf(tmp));
}	

double KDE::pdf(double x, double y){
	vector<double> tmp;
	tmp.push_back(x);
	tmp.push_back(y);
	return(pdf(tmp));
}	

double KDE::pdf(vector<double>& data){

	calc_bandwidth();
  	double d = 0.0;
  	for(unsigned int i = 0; i < data_matrix[0].size(); i++){
  		double a = 1.0;
  		for(curr_var = 0; curr_var < data_matrix.size(); curr_var++){
  			switch(kernel){
  				case 2: a *= box_pdf(data[curr_var],data_matrix[curr_var][i],bandwidth_map[curr_var]);break;
  				case 3: a *= epanechnikov_pdf(data[curr_var],data_matrix[curr_var][i],bandwidth_map[curr_var]);break;
  				default: a *= gauss_pdf(data[curr_var],data_matrix[curr_var][i],bandwidth_map[curr_var]);break;
  			}
      }
  		d += a;
  	}  	
  	return(d/count_map[0]);
}	

double KDE::cdf(double x){
	vector<double> tmp;
	tmp.push_back(x);
	return(cdf(tmp));
}	

double KDE::cdf(double x, double y){
	vector<double> tmp;
	tmp.push_back(x);
	tmp.push_back(y);
	return(cdf(tmp));
}	

double KDE::cdf(vector<double>& data){

	calc_bandwidth();
  	double d = 0.0;
  	for(unsigned int i = 0; i < data_matrix[0].size(); i++){
  		double a = 1.0;
  		for(curr_var = 0; curr_var < data_matrix.size(); curr_var++){
  			switch(kernel){
  				case 2: a *= box_cdf(data[curr_var],data_matrix[curr_var][i],bandwidth_map[curr_var]);break;
  				case 3: a *= epanechnikov_cdf(data[curr_var],data_matrix[curr_var][i],bandwidth_map[curr_var]);break;
  				default: a *= gauss_cdf(data[curr_var],data_matrix[curr_var][i],bandwidth_map[curr_var]);break;
  			}
  		}
  		d += a;
  	}  	
  	return(d/count_map[0]);

}	

void KDE::default_bandwidth(){

	if(!count_map[curr_var]){
		std::cout << "No data!" << std::endl;
		exit(1);	
	}
  	double x  = sum_x_map[curr_var]/count_map[curr_var];
  	double x2 = sum_x2_map[curr_var]/count_map[curr_var];
  	double sigma = sqrt(x2 - (x*x));
  	double b = sigma*(pow((3.0*count_map[curr_var]/4.0),(-1.0/5.0)));

  	if(bandwidth_opt_type == 1) bandwidth_map[curr_var] = b;
  	default_bandwidth_map[curr_var] = b;
}

// Secant method
void KDE::optimal_bandwidth(int maxiters, double eps){

	if(!count_map[curr_var]){
		std::cout << "No data!" << std::endl;
		exit(1);	
	}
	default_bandwidth();
  	double x0 = default_bandwidth_map[curr_var];
  	double y0 = optimal_bandwidth_equation(x0,get_min(curr_var),get_max(curr_var),data_matrix[curr_var]);
  	double x = 0.8*x0;
	double y = optimal_bandwidth_equation(x,get_min(curr_var),get_max(curr_var),data_matrix[curr_var]);
  	int i = 0;

  	while(i < maxiters){
	    x -= y*(x0-x)/(y0-y);
	    y = optimal_bandwidth_equation(x,get_min(curr_var),get_max(curr_var),data_matrix[curr_var]);
	    if(abs(y) < eps*y0){
	    	break;
	    }
	    i++;
  	}
  	bandwidth_map[curr_var] = x;
}

double KDE::optimal_bandwidth_equation(double w, double min, double max, vector<double>& data){

	double alpha = 1.0/(2.0*sqrt(M_PI));
  	double sigma = 1.0;
  	double n = count_map[curr_var];
  	double q = stiffness_integral(w,min,max,data);
	return w - pow(((n*q*pow(sigma,4))/alpha ),(-1.0/5.0));
}

// Calculates the integral over the square of the curvature using 
// the trapezoidal rule. decreases step size by half until the relative
// error in the last step is less eps.

double KDE::stiffness_integral(double w, double mn, double mx, vector<double>& data){

	double eps = 1e-4;
	double n = 1;
	double dx = (mx-mn)/n;
	double curv_mx = curvature(mx,w,data);
	double curv_mn = curvature(mn,w,data);
 	double yy = 0.5*((curv_mx*curv_mx)+(curv_mn*curv_mn))*dx;
  	double maxn = (mx-mn)/sqrt(eps);

  	maxn = maxn > 2048 ? 2048 : maxn;

	for(int n = 2; n <= maxn; n *= 2){
		dx /= 2.0;
		double y = 0.0;
		for(int i = 1; i <= n-1; i +=2){
			curv_mn = pow(curvature(mn + i*dx, w, data),2);
			y += curv_mn;
		}	
		yy = 0.5*yy + y*dx;
		if(n > 8 && abs(y*dx-0.5*yy) < eps*yy){
			break;
		}
	}
  	return(yy);
}

// Bisection method
void KDE::optimal_bandwidth_safe(double eps){

	if(!count_map[curr_var]){
		std::cout << "No data!" << std::endl;
		exit(1);	
	}
	default_bandwidth();
	double x0 = default_bandwidth_map[curr_var]/count_map[curr_var];
	double x1 = 2*default_bandwidth_map[curr_var];	
	double y0 = optimal_bandwidth_equation(x0,min_map[curr_var],max_map[curr_var],data_matrix[curr_var]);
	double y1 = optimal_bandwidth_equation(x1,min_map[curr_var],max_map[curr_var],data_matrix[curr_var]);
	
	if(y0 * y1 >= 0){
	//	std::cerr << "Interval [ f(x0=$x0)=$y0 : f(x1=$x1)=$y1 ] does not bracket root." << std::endl;
	}

	double x = 0.0, y = 0.0;
	int i = 0;

	while(abs(x0 - x1) > eps*x1){
		i += 1;
		x = (x0 + x1 )/2;
		y = optimal_bandwidth_equation(x,min_map[curr_var],max_map[curr_var],data_matrix[curr_var]);

		if(abs(y) < eps*y0){
			break;
		}
		if(y * y0 < 0){
			x1 = x;
			y1 = y;
		}else{
			x0 = x;
			y0 = y;
		}
	}
	bandwidth_map[curr_var] = x;
}

double KDE::curvature(double x, double w, vector<double>& data){

  	double y = 0.0;
  	for(auto it = data.begin(); it != data.end(); it++){
  		y += gauss_curvature(x,*it,w);
  	}  	
  	return(y/count_map[curr_var]);

}

void KDE::add_data(double x){
	vector<double> tmp;
	tmp.push_back(x);
	add_data(tmp);
}

void KDE::add_data(double x, double y){
	//std::cout << x << "," << y << std::endl;
	vector<double> tmp;
	tmp.push_back(x);
	tmp.push_back(y);
	add_data(tmp);
}

void KDE::add_data(vector<double>& x){

	if(!data_matrix.size()){
		for(size_t i = 0; i < x.size(); i++){
			vector<double> tmp;
			tmp.push_back(x[i]);
			data_matrix.push_back(tmp);
			sum_x_map[i] = x[i];
			sum_x2_map[i] = x[i]*x[i];
			count_map[i] = 1;
			min_map[i] = x[i];
			max_map[i] = x[i];
			bandwidth_map[i] = -1.0;
		}
	}else{
		if(x.size() != data_matrix.size()){
			std::cout << "Number of variables doesn't match!" << std::endl;
		}else{
			for(size_t i = 0; i < x.size(); i++){
				data_matrix[i].push_back(x[i]);
				sum_x_map[i] += x[i];
				sum_x2_map[i] += x[i]*x[i];
				count_map[i]++;			
				min_map[i] = x[i] < min_map[i] ? x[i] : min_map[i];
				max_map[i] = x[i] > max_map[i] ? x[i] : max_map[i];
				bandwidth_map[i] = -1.0;
			}				
		}
	}	
}

double KDE::gauss_cdf(double x, double m, double s){

	// Abramowitz Stegun Normal Cumulative Distribution Function	
	double z = abs(x - m)/s;
  	double t = 1.0/(1.0 + 0.2316419*z);
  	double y = t*( 0.319381530 + t*(-0.356563782 + t*(1.781477937 + t*(-1.821255978 + t*1.330274429 ))));
	if(x >= m){
		return 1.0 - gauss_pdf(x,m,s)*y*s;
	}else{
		return gauss_pdf(x,m,s)*y*s;
	}

}

double KDE::gauss_pdf(double x, double m, double s){
	
	double z = (x - m)/s;
	return exp(-0.5*z*z)/(s*sqrt( 2.0*M_PI));

}

double KDE::gauss_curvature(double x, double m, double s){
	
	double z = (x - m)/s;
  	return ((z*z) - 1.0)*gauss_pdf(x,m,s)/(s*s);

}

double KDE::box_pdf(double x, double m, double s){

	if(x < m-0.5*s || x > m+0.5*s){
		return (0.0);
	}else{
		return (1.0/s);
	}
}	
double KDE::box_cdf(double x, double m, double s){

	if(x < m-0.5*s){
		return (0.0);
	}
	if( x > m+0.5*s){
		return 1.0;
	}
	return (x-m)/s + 0.5;

}	

double KDE::epanechnikov_pdf(double x, double m, double s){

	double z = (x-m)/s;
	if(fabs(z) > 1.0){
		return (0.0);
	}
	return 0.75*(1.0-(z*z))/s;

}	

double KDE::epanechnikov_cdf(double x, double m, double s){

	double z = (x-m)/s;
	if(z < -1.0){
		return (0.0);
	}
	if(z >  1.0){
		return (1.0);
	}
	return 0.25*(2.0 + 3.0*z - (z*z*z));
}	
