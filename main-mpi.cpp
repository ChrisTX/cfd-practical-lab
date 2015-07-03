#include "cfd_solver_MPI.hpp"
#include <mpi.h>
#include <stdexcept>

int main(int argc, char** argv) {
	if( MPI_Init(&argc, &argv) != MPI_SUCCESS )
		throw std::runtime_error{"MPI_Init failed!"};

	using cfdsolv = CFD::CFDSolverMPI<double>;
	using bedge = cfdsolv::BoundaryEdge_t;
	using bcond = cfdsolv::BoundaryCondition_t;

	const auto sizeI = 100;
	const auto sizeJ = 100;

	cfdsolv::numeric_parameters np{ 0.5, 1.7, 0.0001, 1000 };
	cfdsolv::problem_data pd{ 10, 0, 0 };
	cfdsolv::time_parameters tp{ 0.5, 0.0, 0.02, 2.0, 2.0 };
	cfdsolv::problem_geometry pg { 10, 10, sizeI, sizeJ };

	int numnodes;
	if( MPI_Comm_size(MPI_COMM_WORLD, &numnodes) != MPI_SUCCESS)
		throw std::runtime_error{"MPI_Comm_size failed!"};

	int dims[2] = { 0, 0 };
	if( MPI_Dims_create(numnodes, 2, dims) != MPI_SUCCESS )
		throw std::runtime_error{"MPI_Dims_create failed!"};

	int isperiodic[2] = { 0, 0 };
	MPI_Comm commcart;
	if( MPI_Cart_create(MPI_COMM_WORLD, 2, dims, isperiodic, 1, &commcart) != MPI_SUCCESS )
		throw std::runtime_error{"MPI_Cart_create failed!"};

	int mpirank;
	if( MPI_Comm_rank(commcart, &mpirank) != MPI_SUCCESS )
		throw std::runtime_error{"MPI_Comm_rank failed!"};

	cfdsolv::MPI_information mpii{ mpirank, commcart, static_cast<cfdsolv::size_type>( dims[0] ), static_cast<cfdsolv::size_type>( dims[1] ) };
	
	cfdsolv cfdsv(mpii, pd, pg, tp, np, [](double x, double y) -> double { return 0.; }, [](double x, double y) -> double { return 0.; }, [](double x, double y) -> double { return 0.; });
	cfdsv.SetBoundaryCondition(bedge::left, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::right, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::bottom, bcond::noslip);
	cfdsv.SetBoundaryCondition(bedge::top, bcond::inflow, [](double x, double y) -> auto {
		return std::make_pair(1., 0.);
	});
	cfdsv.Solve();

	if( MPI_Finalize() != MPI_SUCCESS )
		throw std::runtime_error{"MPI_Finalize failed!"};
}
