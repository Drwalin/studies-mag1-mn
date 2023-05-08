
#include <cstdio>
#include <cmath>
#include <cstdlib>

using Real = double;
struct Pair {
	Real x, y;
};
Real fabs(Real v) {
	if(v < 0) return -v; else return v;
}



Real FValue(Pair p) {
	return p.y*p.y + 1;
}
Real FValue(Real x, Real y) {
	return FValue({x, y});
}



Pair SolverRKF(Pair old, Real& dx) {
	Real h = dx;
	const Real x0 = old.x;
	const Real y0 = old.y;
	constexpr Real epsilon = 0.00001;
	Real f0 = FValue(old);
	Real f1 = FValue(
			x0 + h/4,
			y0 + h*f0/4);
	Real f2 = FValue(
			x0 + 3*h/8,
			y0 + 3*h*f0/32 + 9*h*f1/32);
	Real f3 = FValue(
			x0 + 12*h/13,
			y0 + 1932*h*f0/2197 - 7200*h*f1/2197 + 7296*h*f2/2197);
	Real f4 = FValue(
			x0 + h,
			y0 + 439*h*f0/216 - 8*h*f1 + 3680*h*f2/513 - 845*h*f3/4104);
	Real f5 = FValue(
			x0 + h/2,
			y0 - 8*h*f0/27 + 2*h*f1 - 3544*h*f2/2565 + 1859*h*f3/4104 - 11*h*f4/40);
	Real error = h*(f0/360 - 128*f2/4275 - 2197*f3/75240 + f4/50 + 2*f5/55);
	dx = 0.9*h*pow(fabs(h)*epsilon/fabs(error), 0.25);
	if(fabs(error) > h*epsilon) {
		return SolverRKF(old, dx);
	}
	Pair next = {
		x0+h,
		y0 + h*(16*f0/135 + 6656*f2/12825 + 28561*f3/56430 - 9*f4/50 + 2*f5/55)
	};
	return next;
}


Pair SolveEuler(Pair p, Real& dx) {
	Pair r;
	r.x = p.x + dx;
	Real dfdx = p.y*p.y + 1;
	r.y = p.y + dfdx*dx;
	return r;
}

Pair SolveRungeKnutta(Pair p, Real& dx) {
	Real _dx = dx*0.5f;
	Pair m = SolveEuler(p, _dx);
	Pair r;
	r.x = p.x + dx;
	Real dfdx = m.y*m.y + 1;
	r.y = p.y + dfdx*dx;
	return r;
}




Pair* SolveAdamsBashforth2(Pair p0, Real dx, Real maxx) {
	size_t num = (maxx-p0.x)/dx+1;
	Pair *p = new Pair[num+1];
	p[0] = p0;
	for(int i=1; i<=num; ++i) {
		p[i].x = p[i-1].x + dx;
	}
	p[1].y = p[0].y + dx * FValue(p[0].x + dx/2, p[0].y + dx/2*FValue(p[0]));
	for(int i=2; i<=num; ++i) {
		Pair p1=p[i-2], p2=p[i-1];
		Real y = p2.y + dx/2*(
				3*FValue(p2) - FValue(p1)
			);
		p[i].y = y;
	}
	return p;
}

Pair* SolveAdamsBashforth4(Pair p0, Real dx, Real maxx) {
	size_t num = (maxx-p0.x)/dx+1;
	Pair *p = new Pair[num+1];
	p[0] = p0;
	for(int i=1; i<=num; ++i) {
		p[i].x = p[i-1].x + dx;
	}
	for(int i=1; i<4; ++i) {
		Pair o = p[i-1];
		Real x0 = o.x;
		Real y0 = o.y;
		Real f0, f1, f2, f3;
		f0 = FValue(o);
		f1 = FValue(x0 + dx/2, y0 + dx*f0/2);
		f2 = FValue(x0 + dx/2, y0 + dx*f1/2);
		f3 = FValue(x0 + dx, y0 + dx*f2);
		p[i].y = 
			y0 + dx/6*(f0 + 2*f1 + 2*f2 + f3);
	}
	for(int i=4; i<=num; ++i) {
		Pair pn=p[i-1], pn1=p[i-2], pn2=p[i-3], pn3=p[i-4];
		Real y = pn.y + dx/24*(
				55*FValue(pn) - 59*FValue(pn1) + 37*FValue(pn2) - 9*FValue(pn3)
			);
		p[i].y = y;
	}
	return p;
}




using SolutionType = Pair (*)(Pair, Real&);

void PrintSolutions(Real x0, Real y0, Real _dx, Real maxx, int solutionsCount, SolutionType solutions[]) {
	printf("y(x_0 = %f) = %f\n", (float)x0, (float)y0);
	printf(" dx = %f\n", (float)_dx);
	
	Pair* p = new Pair[solutionsCount];
	Real* dx = new Real[solutionsCount];
	for(int i=0; i<solutionsCount; ++i) {
		dx[i] = _dx;
	}
	for(;;) {
		printf(" y( %f ) = ", (float)p[0].x);
		for(int i=0; i<solutionsCount; ++i) {
			printf("\t%f", (float)p[i].y);
			p[i] = solutions[i](p[i], dx[i]);
		}
		printf("\n");
		if(p[0].x > maxx)
			break;
	}
	delete[] p;
	delete[] dx;
	printf("\n\n");
}

void Solve(Real x0, Real y0, Real maxx, int solutionsCount, SolutionType solutions[]) {
	PrintSolutions(x0, y0, 0.2, maxx, solutionsCount, solutions);
// 	PrintSolutions(x0, y0, 0.1, maxx, solutionsCount, solutions);
// 	PrintSolutions(x0, y0, 0.05, maxx, solutionsCount, solutions);
// 	PrintSolutions(x0, y0, 0.025, maxx, solutionsCount, solutions);
}

int main() {
	Pair p0{0, 0};
	Pair *a, *b;
	
	a = SolveAdamsBashforth2(p0, 0.05, 1);
	b = SolveAdamsBashforth4(p0, 0.05, 1);
	
	int count = 1/0.05+1;
	
	for(int i=0; i<count; ++i) {
		printf(" y(%f) = \t%f    |    \t%f\n", a[i].x, a[i].y, b[i].y);
	}
	
	delete[] a;
	delete[] b;
}

