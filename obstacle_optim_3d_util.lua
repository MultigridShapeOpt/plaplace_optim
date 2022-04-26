
-- Create util namespace
util = util or {}


util.oo = util.oo or {}

--change to baseSolver=LU(),
function util.oo.linear_solver(domainDisc, approxSpace, vrb)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			="gs",
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = SuperLU(),--"lu",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap                 = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 3000,		-- number of iterations
			absolute	= 1e-10,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 0.0,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose	= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	LinSolver=util.solver.CreateSolver(LinSolverDesc)
	return LinSolver
end

function util.oo.ns_solver(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			=ComponentGaussSeidel(0.01, {"p"}, {0}, {0.9}),
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = SuperLU(),--"superlu",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 20000,		-- number of iterations
			absolute	= 1e-14,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-12,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	NonLinSolverDesc = 
	{
	  type    = "newton",
	  convCheck = {
			type		= "standard",
			iterations	= 50,		-- number of iterations
			absolute	= 1e-12,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 0.0,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   },
	  lineSearch  = 
	  {
	  type      = "standard",
		maxSteps      = 50,
		lambdaStart   = 1,
		lambdaReduce  = 0.9,
		acceptBest    = true,
		checkAll      = false,
		verbose		  = false
	  },
	  linSolver	= LinSolverDesc
	}
	non_LinSolver=util.solver.CreateSolver(NonLinSolverDesc)
	return non_LinSolver
end

function util.oo.adjoint_ns_solver(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			=ComponentGaussSeidel(0.01, {"h"}, {0}, {0.9}),
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = SuperLU(),   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 20000,		-- number of iterations
			absolute	= 1e-12,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 0.0,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	LinSolver=util.solver.CreateSolver(LinSolverDesc)
	return LinSolver
end

function util.oo.gmg(approxSpace, smoother, numPreSmooth, numPostSmooth,
						 cycle_str, baseSolver, baseLev, bRAP)

	local gmg = GeometricMultiGrid(approxSpace)
	
	gmg:set_base_level(baseLev)
	gmg:set_base_solver(baseSolver)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	gmg:set_smoother(smoother)
	gmg:set_cycle_type(cycle)
	gmg:set_num_presmooth(numPreSmooth)
	gmg:set_num_postsmooth(numPostSmooth)
	gmg:set_rap(bRAP)
	
	return gmg

end

function util.oo.linear_solver_damping(domainDisc, approxSpace, vrb)

	gmg = GeometricMultiGrid(approxSpace)
	smthr=GaussSeidel();smthr:set_damp(1.83)
	gmg:set_base_level(0)
	gmg:set_base_solver(SuperLU())
	gmg:set_gathered_base_solver_if_ambiguous(false)
	gmg:set_smoother(smthr)
	gmg:set_cycle_type("V")
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_rap(true)
	gmg:set_discretization(domainDisc)

	LinSolver=BiCGStab();LinSolver:set_preconditioner(gmg)
	LinSolver:set_convergence_check(ConvCheck(2000,1.0e-12,0.0,true))

	return LinSolver;
end
