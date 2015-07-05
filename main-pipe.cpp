#include "cfd_solver.hpp"

int main() {
	using cfdsolv = CFD::CFDSolver<double>;
	using bedge = cfdsolv::BoundaryEdge_t;
	using bcond = cfdsolv::BoundaryCondition_t;
	cfdsolv::numeric_parameters np{ 0.7, 1.7, 0.001, 100 };
	cfdsolv::problem_data pd{ 1000, 0, 0 };
	cfdsolv::time_parameters tp{ 0.5, 0.0, 0.05, 100, 10 };
	const std::size_t IMAX = 60;
	const std::size_t JMAX = 60;	
	cfdsolv::problem_geometry pg { 12, 12, IMAX, JMAX };
	using vec_t = cfdsolv::VectorType;
	vec_t u(IMAX + 2, JMAX + 2, 1);
	for(std::size_t j = 0u; j < JMAX/2u + 1u; ++j)
		for(std::size_t i = 0u; i < IMAX + 1u; ++i)
			u(i, j) = 0;

	vec_t v(IMAX + 2, JMAX + 2, 0);
	vec_t p(IMAX + 2, JMAX + 2, 0);
	cfdsolv cfdsv(pd, pg, tp, np, u, v, p);
	cfdsv.SetBoundaryCondition(bedge::left, bcond::outflow);
	cfdsv.SetBoundaryCondition(bedge::right, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::bottom, bcond::outflow);
	cfdsv.SetBoundaryCondition(bedge::top, bcond::noslip);

	using boolvec_t = LinearAlgebra::Vector2D<bool>;
	boolvec_t myobs(IMAX + 1, JMAX + 1, true);
	for(std::size_t i = 1u; i < IMAX - JMAX/3; ++i) {
		for(std::size_t j = 1u; j < JMAX/3u + 1; ++j)
			myobs(i, j) = false;
		for(std::size_t j = (2u * JMAX)/3u + 1; j < JMAX + 1; ++j)
			myobs(i, j) = false;
	}

	cfdsv.SetObstacles(myobs);
	cfdsv.Solve();
}
