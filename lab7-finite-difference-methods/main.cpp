
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>


using Real = float;

Real fabs(Real v) {
	if(v < 0) return -v; else return v;
}


float SolveOnce(std::vector<Real>& v, const float h, const bool gaussseidel) {
	Real px = v[0];
	Real maxErr = 0, x = h;
	for(int i=1; i<v.size()-1; ++i, x+=h) {
		Real py = v[i];
		const Real y = v[i] = (
					(1.0-5.0*h/2.0)*v[i+1] +
					(1.0+5.0*h/2.0)*px -
					10.0*h*h*x
				) / (2.0 - 10.0*h*h);
		maxErr = std::max(maxErr, fabs((y - py)/y));
		px = gaussseidel ? y : py;
	}
	return maxErr;
}

void Init(std::vector<Real>& values) {
	values.resize(11);
	for(int i=0; i<values.size(); ++i) {
		values[i] = i*10;
	}
}

int Solve(std::vector<Real>& v, const float h, const bool gaussseidel) {
	Init(v);
	for(int i=0; i<10000; ++i) {
		Real err = SolveOnce(v, h, gaussseidel);
		if(err <= 0.0005)
			return i+1;
	}
	return -1;
}

void Print(std::vector<Real>& v, Real h) {
	printf("\n");
	for(int i=0; i<v.size(); ++i) {
		printf(" y(%.1f) = %f\n", i*h, v[i]);
	}
}



int main() {
	Real h = 0.1;
	std::vector<Real> v;
	
	int steps = Solve(v, h, false);
	printf("Jacobi [%i]:", steps);
	Print(v, h);
	
	steps = Solve(v, h, true);
	printf("Gauss-Seidel [%i]:", steps);
	Print(v, h);
	
	return 0;
}

