#include "cfd_solver.hpp"

int main() {
	using cfdsolv = CFD::CFDSolver<double>;
	using bedge = cfdsolv::BoundaryEdge_t;
	using bcond = cfdsolv::BoundaryCondition_t;
	cfdsolv::numeric_parameters np{ 0.9, 1.7, 0.001, 500 };
	cfdsolv::problem_data pd{ 100, 0, 0 };
	cfdsolv::time_parameters tp{ 0.5, 0.0, 0.05, 100, 10 };
	const std::size_t IMAX = 100;
	const std::size_t JMAX = 20;	
	cfdsolv::problem_geometry pg { 10, 2, IMAX, JMAX };
	using vec_t = cfdsolv::VectorType;
	vec_t u(IMAX + 2, JMAX + 2, 1);
	for(std::size_t j = 0u; j < JMAX/2u + 1u; ++j)
		for(std::size_t i = 0u; i < IMAX + 1u; ++i)
			u(i, j) = 0;

	vec_t v(IMAX + 2, JMAX + 2, 0);
	vec_t p(IMAX + 2, JMAX + 2, 0);
	cfdsolv cfdsv(pd, pg, tp, np, u, v, p);
	cfdsv.SetBoundaryCondition(bedge::left, bcond::outflow);
	cfdsv.SetBoundaryCondition(bedge::right, bcond::outflow);
	cfdsv.SetBoundaryCondition(bedge::bottom, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::top, bcond::noslip);

	using boolvec_t = LinearAlgebra::Vector2D<bool>;
	boolvec_t myobs(IMAX + 1, JMAX + 1, true);
	for(std::size_t j = 1u; j < JMAX/2u + 1u; ++j)
		for(std::size_t i = 1u; i < JMAX/2u + 1u; ++i)
			myobs(i, j) = false;

	cfdsv.SetObstacles(myobs);
	cfdsv.Solve();
}
