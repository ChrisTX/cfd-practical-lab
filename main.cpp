#include "cfd_solver.hpp"

int main() {
	using cfdsolv = CFD::CFDSolver<double>;
	using bedge = cfdsolv::BoundaryEdge_t;
	using bcond = cfdsolv::BoundaryCondition_t;
	cfdsolv::numeric_parameters np{ 0.5, 1.7, 0.001, 100 };
	cfdsolv::problem_data pd{ 10, 0, 0 };
	cfdsolv::time_parameters tp{ 0.5, 0.0, 0.02, 2.0, 2.0 };
	cfdsolv::problem_geometry pg { 10, 10, 50, 50 };
	using vec_t = cfdsolv::VectorType;
	vec_t u(52, 52, 0.);
	vec_t v(52, 52, 0.);
	vec_t p(52, 52, 0.);
	cfdsolv cfdsv(pd, pg, tp, np, u, v, p);
	cfdsv.SetBoundaryCondition(bedge::left, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::right, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::bottom, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::top, bcond::inflow, [](double x, double y) -> auto {
		return std::make_pair(1., 0.);
	});
	cfdsv.Solve();
}
