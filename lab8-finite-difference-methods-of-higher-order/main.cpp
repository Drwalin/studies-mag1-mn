
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using Real = double;

Real fabs(Real v) {
	if(v < 0) return -v; else return v;
}

Real SolveOnce(std::vector<Real>& v, const Real h) {
	Real px = v[0];
	Real maxErr = 0, x = h;
	int n = v.size()-1;
	for(int i=1; i<n; ++i, x+=h) {
		Real py = v[i];
		if(i==1) {
			v[i] = (
					-v[i+4] +
					(6+5*h) * v[i+3] +
					(-14-30*h) * v[i+2] +
					(4+90*h) * v[i+1] +
					(-10-15*h) * v[i-1] +
					120*h*h*x
				) / (-15 + 50*h + 120*h*h);
		} else if(i==n-1) {
			v[i] = (
					-v[i-4] +
					(6-5*h) * v[i-3] +
					(-14+30*h) * v[i-2] +
					(4-90*h) * v[i-1] +
					(-10+15*h) * v[i+1] +
					120*h*h*x
				) / (-15 - 50*h + 120*h*h);
		} else {
			v[i] = (
					(1-5*h) * v[i+2] +
					(-16+40*h) * v[i+1] +
					(-16-40*h) * v[i-1] +
					(1+5*h) * v[i-2] +
					120*h*h*x
				) / (-30 + 120*h*h);
		}
		maxErr = std::max(maxErr, fabs((v[i] - py)/v[i]));
		px = v[i];
	}
	return maxErr;
}

void Init(std::vector<Real>& values, Real h) {
	values.resize(21);
	for(int i=0; i<values.size(); ++i) {
		values[i] = 100*i/(values.size()-1);
	}
}

int Solve(std::vector<Real>& v, const Real h) {
	Init(v, h);
	for(int i=0; i<10000; ++i) {
		Real err = SolveOnce(v, h);
		if(err <= 0.0000005)
			return i+1;
	}
	return -1;
}

void Print(std::vector<Real>& v, Real h) {
	printf("\n");
	for(int i=0; i<v.size(); ++i) {
		printf(" y(%.2f) = %f\n", i*h, v[i]);
	}
}

int main() {
	Real h = 0.05;
	std::vector<Real> v;
	
	int steps = Solve(v, h);
	printf("Gauss-Seidel [%i]:", steps);
	Print(v, h);
	
	return 0;
}

