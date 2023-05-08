
#include <cstdio>
#include <cmath>

using Real = long double;

struct Pair {
	Real x, y;
};

Pair SolveEuler(Pair p, Real dx) {
	Pair r;
	r.x = p.x + dx;
	Real dfdx = p.y*p.y + 1;
	r.y = p.y + dfdx*dx;
	return r;
}

Pair SolveRungeKnutta(Pair p, Real dx) {
	Pair m = SolveEuler(p, dx*0.5f);
	Pair r;
	r.x = p.x + dx;
	Real dfdx = m.y*m.y + 1;
	r.y = p.y + dfdx*dx;
	return r;
}

using SolutionType = Pair (*)(Pair, Real);

void PrintSolutions(Real x0, Real y0, Real dx, Real maxx, int solutionsCount, SolutionType solutions[]) {
	printf("y(x_0 = %f) = %f\n", (float)x0, (float)y0);
	printf(" dx = %f\n", (float)dx);
	
	Pair* p = new Pair[solutionsCount];
	for(; p[0].x < maxx+0.00001;) {
		printf(" y( %f ) = ", (float)p[0].x);
		for(int i=0; i<solutionsCount; ++i) {
			printf("\t%f", (float)p[i].y);
			p[i] = solutions[i](p[i], dx);
		}
		printf("\n");
	}
	printf("\n\n");
}

void Solve(Real x0, Real y0, int solutionsCount, SolutionType solutions[]) {
	PrintSolutions(x0, y0, 0.2, 1, solutionsCount, solutions);
	PrintSolutions(x0, y0, 0.1, 1, solutionsCount, solutions);
	PrintSolutions(x0, y0, 0.05, 1, solutionsCount, solutions);
	PrintSolutions(x0, y0, 0.025, 1, solutionsCount, solutions);
}

int main() {
	Real x0 = 0, y0 = 0;
	SolutionType arr[] = {SolveEuler, SolveRungeKnutta};
	Solve(x0, y0, 2, arr);
	return 0;
}

