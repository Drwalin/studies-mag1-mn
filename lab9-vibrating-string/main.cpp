
#include <cmath>
#include <map>
#include <cstdio>

using Real = float;
using P = std::pair<int, int>;

class Matrix : public std::map<P, Real> {
public:
	using parent = std::map<P, Real>;
	
	Real& operator[](P p) {
		_width = std::max(std::max(p.first, p.second), _width);
		return parent::operator[](p);
	}
	Real operator[](P p) const {
		return std::map<P, Real>::at(p);
	}
	Real operator()(int x, int y) const {
		return this->operator[]({x, y});
	}
	Real& operator()(int x, int y) {
		return this->operator[]({x, y});
	}
	
	void clear() {
		parent::clear();
		_width = 0;
	}
	
	int width() const { return _width; }
	
	int _width;
};


Real Det(const Matrix& A, Real lambda) {
	int width = A.width();
	Real D1 = 0, D2 = 0, D = 0;
	D1 = A(0,0) - lambda;
	D2 = (A(1,1) - lambda) * D1 - A(0,1)*A(1,0);
	for(int i=2; i<=width; ++i) {
		D = (A(i,i)-lambda)*D2 - A(i-1,i)*A(i,i-1)*D1;
		D1 = D2;
		D2 = D;
	}
	return D;
}

void ConstructMatrix(Matrix& A, Real T, Real L, Real h, int N, Real(*u)(Real x)) {
	A.clear();
	auto ui = [=](int i){
		return u(L * (i+1) / N);
	};
	for(int i=0; i<N-1; ++i) {
		A(i,i) = 2.0f*T/ui(i);
		if(i+1<N-1) {
			A(i+1, i) = -T/ui(i);
			A(i, i+1) = -T/ui(i+1);
		}
	}
}

Real Bisection(int steps, const Matrix& A, Real minl, Real maxl) {
	Real solmin, solmax;
	solmin = Det(A, minl);
	solmax = Det(A, maxl);
// 	printf(" solmin = %2.2f, solmax = %2.2f\n", (float)solmin, (float)solmax);
	for(int i=0; i<steps; ++i) {
		const Real mid = (minl+maxl) * 0.5;
		const Real sol = Det(A, mid);
		if(solmin * sol > 0) {
			solmin = sol;
			minl = mid;
		} else {
			solmax = sol;
			maxl = mid;
		}
// 		printf("Det(lambda=%4.4f) = %4.4f\n", (float)mid, (float)sol);
	}
	Real mid = (minl+maxl) * 0.5;
// 	printf("Det(lambda=%2.2f) = %2.2f\n", (float)mid, (float)Det(A, mid));
	return mid;
}

Real AnalyticalSolutionBase(Real L, Real T, Real(*u)(Real x)) {
	return M_PI * sqrt(T/u(0)) / L;
}

Real FindSolution(Matrix& A, Real h, Real L, Real T, Real N, int n,
		Real(*u)(Real x), Real(*u0)(Real x)) {
	Real analyticalSol = AnalyticalSolutionBase(L, T, u0);
	Real ln = pow(n*analyticalSol*h, 2);
	Real lnnext = (pow((n+1)*analyticalSol*h, 2)*0.2 + 0.8*ln);
	Real lnprev = (pow((n-1)*analyticalSol*h, 2)*0.2 + 0.8*ln);
	
	
	float margin = 0.0001;
	
	Real lambda = Bisection(100, A, std::min(lnprev, lnnext)+margin, std::max(lnprev, lnnext)-margin);
// 	printf(" ln(i) = %f  <=>  %f (found)\n", n, (float)ln, (float)lambda);
	Real w = sqrt(lambda) / h;
	return w;
}

int main() {
	auto u1 = [](Real x)->Real{ return 1.0f; };
	auto u2 = [](Real x)->Real{ return 1.0f + x - 0.5f; };
	
	Matrix A;
	Real h = 0.01;
	Real L = 1;
	Real T = 1;
	Real N = 100;
	ConstructMatrix(A, T, L, h, N, u1);
	Real(*u)(Real x) = u1;
	
	printf(" solutions for constant u:\n");
	for(int i=1; i<=5; ++i) {
		float soli = FindSolution(A, h, L, T, N, i, u, u1);
		printf("    sol %i = %2.2f\n            %2.2f\n", i, soli, i*(float)AnalyticalSolutionBase(L, T, u));
	}
	
	u = u2;
	ConstructMatrix(A, T, L, h, N, u2);
	printf(" solutions for non-const u:\n");
	for(int i=1; i<=5; ++i) {
		float soli = FindSolution(A, h, L, T, N, i, u2, u1);
		printf("    sol %i = %2.2f\n", i, soli);
	}
	
	
	
	
	
	
// 	Real lambda = Bisection(100, A, -5, 5);
// 	Real w = sqrt(lambda) / h;
// 	
// 	printf(" solution omega = %2.2f\n", (float)w);

	
	
	
	
	
	return 0;
}

