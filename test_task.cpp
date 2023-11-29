#include <iostream>
#include <vector>
#include <complex>

using namespace std;

const int M = 16; // Number of points in the constellation 
const double a = 1.; 
const double pi = 3.14159265358979323846;
const int N = 300; // Number of grid points used for computation of the entropy H(y). Default N = 300 is enough to reach required numerical accuracy

template <typename T> T sqr(T a) {
    return a*a;
}

struct Point { // Structure class of the constellation points
    complex<double> c;
    
    double prob(complex<double> y, double sigma);
};

double Point::prob(complex<double> y, double sigma) {
	/* 2D Gaussian probability density about the constellation point,
	i.e. the product of 1D Gaussian probability distributions */	
    return exp(-.5*sqr(abs(y-c)/sigma))/(2.*pi*sqr(sigma));
};


double abs(const Point &p) {
    return abs(p.c);
} 

double pr(complex<double> z, vector<Point> points, double sigma) {
	/* 2D Gaussian probability distribution about the constellation point,
	i.e. the product of 1D Gaussian probability distributions */	
    double Pr = 0.;
    for (auto &c : points) {
        Pr += c.prob(z, sigma);		
    };	
	return Pr;
};

double entropy(vector<Point> points, double sigma) {
	// Integral for entropy H(Y) on a 2D Cartesian grid
	
    double H = 0.; // Entropy to be computed
	double step = (4.*a + 7.*sigma)/N; // Distance between two next grip points
		for (int ix = -N; ix <= N; ix++) {
			double x = step*ix;
			for (int iy = -N; iy <= N; iy++) {
				double y = step*iy;
				complex<double> z(x, y);
				H += pr(z, points, sigma)*log2(pr(z, points, sigma));
		    }
	    }
	
	return -H*sqr(step);
}
	
int main() {
       
    vector<Point> points;
	points.resize(M);
	 	
	int key {1}; // Choose constellation: 1 or 2
	switch(key)
	{
        case 1: // QAM constellation
            for (int i = 0; i < 4; i++) {
		        for (int j = 0; j < 4; j++) {
			    complex<double> p(a*(double)(2*i-3),a*(double)(2*j-3));
		        points[4*i + j] = {p};
		        }
	        }
			cout << "QAM-16 constellation selected." << endl;
            break;
			
        case 2: // AQAM constellation
	
            for (int i = 0; i < 4; i++) {
		        for (int j = 0; j < 4; j++) {
			    double rad = 1. + (pow(2.,.5)-1.)*i;
			    complex<double> p(cos(.25*pi*(2.*j+(i%2))),sin(.25*pi*(2.*j+(i%2))));
		        points[4*i + j] = {2.*a*p*rad};
		        }
	        }
			cout << "AQAM-16 constellation selected." << endl;
            break;
    }
	
    double E = 0.; // Mean energy
	double sigma = 0.;
	double max = 0; // max |c_k|
	for (int i = 0; i < points.size(); i++) {
        const auto &c = points[i];
		if (abs(c) >= max) {
			max = abs(c);
		};
		E += sqr(abs(c));
	}
    E = E/M;
	
	auto PARP = 10.*log10(sqr(max)/E);
	cout << "E =" << "\t" << E << endl;
	cout << "PARP =" << "\t" << PARP << " dB" << endl << endl;
	
	for (int i = -2; i <= 40; i++) { // Loop over SNR
		sigma = pow(.5*E, .5)*pow(10., -(double) i/40.); // Evaluation of sigma via SNR
		double Inf = log2((double)M) - log2(2.*pi*exp(1.)*sqr(sigma)) + entropy(points, sigma)/(double)M; // Finally, mutual information 
	    cout << .5*i << "\t" << Inf << "\t" << log2(1.+ pow(10., (double) i/20.)) << endl;
	}
	
    return 0;
}
