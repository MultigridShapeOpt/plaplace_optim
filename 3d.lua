PluginRequired("PLaplacian")
PluginRequired("FluidOptim")
-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")
ug_load_script("util/conv_rates_static.lua")
ug_load_script("util/solver_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("obstacle_optim_3d_util.lua")

ProfileLUA(true)
PrintBuildConfiguration()

package.path = package.path .. ";../../apps/plaplace_optim/lua-matrix/?.lua"
local matrix = require "matrix"

--Constants
dim = 3
numPreRefs=0
vorder=2
porder=1
diameter=6
m=dim+1
flowNames="v1,v2,v3"
pressureName="p"
adjointFlowNames="q1,q2,q3"
adjointPressure="h"
flowOutputFile="vtk_flows"
adjointFlowOutputFile="vtk_adjointFlows"

nodalFlows="flows"
adjointNodalFlows="adjointFlows"


--Simulation Parameters
numRefs = util.GetParamNumber("-numRefs", 3)--refinements after grid partitioning/distribution
numSteps = util.GetParamNumber("-numSteps", 400)
visc = util.GetParamNumber("-visc", 0.02)--medium viscosity
stab = util.GetParamNumber("-stab", 0.0)--stabilization for navier stokes
stabType = util.GetParamNumber("-stabType",0.0)
epsilon_nnl = util.GetParamNumber("-nnl_eps",1.0e-8)--small number for the plaplacian singularity
epsilon_elliptic = util.GetParamNumber("-lin_eps",1.0e-8)--small number for the plaplacian singularity
step_length=util.GetParamNumber("-step_length",1.0)--small number for Lu sensitivity
p_init = util.GetParamNumber("-p_init",2.0)--initial value for plaplacian p
p_max = util.GetParamNumber("-p_max",4.8)--initial value for plaplacian p
p_increase = util.GetParamNumber("-p_inc", 0.19)--increase interval for p
step_control = util.GetParamNumber("-control", 0.05)--control for p term
lineSearchParam=util.GetParamNumber("-line_search", 0.0001)
gridName = util.GetParam("-grid", "./grids/box_3D_elongated.ugx")
bReadInitialGuess = util.GetParamNumber("-restart", -1)
--Newton solver settings
nsMaxIts = util.GetParamNumber("-nsMaxIts", 30)
nsTol = util.GetParamNumber("-nsTol", 1e-6)

lambda_vol = util.GetParamNumber("-lambda_vol", 0.0)
lambda_x   = util.GetParamNumber("-lambda_x", 0.0)
lambda_y   = util.GetParamNumber("-lambda_y", 0.0)
lambda_z   = util.GetParamNumber("-lambda_z", 0.0)

--Boolean parameters
bNewtonOutput = util.GetParamBool("-bNewtonOutput", false)--output newton method per step convergence information
bOutputPProblem  = util.GetParamBool("-bOutputPProblem",false)--output VTK visualization files p-laplacian
bOutputFlows  = util.GetParamBool("-bOutputFlows",true)--output VTK visualization files flow and pressure
bOutputPressure  = util.GetParamBool("-bOutputPressure",false)--output VTK visualization files flow and pressure
bOutputAdjoints  = util.GetParamBool("-bOutputAdjoints",false)--output VTK visualization files of flow and pressureadjoints
bDebugOutput = util.GetParamBool("-bDebugOutput",false)--output VTK visualization files of Lu and RHS of big problem
bDebugNodalPositions = util.GetParamBool("-bDebugNodalPositions",false)
bDebugSensitivity = util.GetParamBool("-bDebugSensitivity",false)--output VTK files for J'
bDoNothing = util.GetParamBool("-bDoNothing",true)--do nothing on flow outlet
bOutputIntermediateUp=util.GetParamBool("-bOutputIntermediateUp",false)--output VTK files for intermediate u_p files
bActivateProfiler = util.GetParamBool("-bActivateProfiler",true)

print("THE PARAMETERS USED FOR EXECUTION ARE: ")
print("grid: "..gridName)
print("numPreRefs:   ".. numPreRefs)
print("numRefs:      ".. numRefs)
print("numSteps:     ".. numSteps)
print("press.grad. stabilization:".. stab)
print("average.based stab:".. stabType)
print("viscosity:".. visc)
print("velocity order "..vorder)
print("pressure order "..porder)
print("epsilon nonlinear:     ".. epsilon_nnl)
print("epsilon linear:     ".. epsilon_elliptic)
print("p_init:     ".. p_init)
print("p_max:     ".. p_max)
print("p_inc:     ".. p_increase)
print("step_control:     ".. step_control)
print("step_length:     ".. step_length)
print("line search parameter: "..lineSearchParam)
print("Newton solver tolerance set to"..nsTol)
print("restarted at step"..bReadInitialGuess)



p_current = p_init
step_control_init = step_control
step_length_init = step_length
function GenerateString(n)
	local s = tostring(n)
	local is =nil
	local ds =nil
	is,ds = string.match(s, "([^.]*)%.([^.]*)")

	if n % 2 ==0 or n % 3 == 0 or n % 4 == 0 then
		is = s
		ds = 0
	end	
	return is.."_"..ds
end



-- initialize ug with the world dimension and the algebra type
InitUG(dim, AlgebraType("CPU", 1));

-- load grid into domain
dom = Domain()
LoadDomain(dom, gridName)
--dom = util.CreateDomain(gridName, 0, {})
print("Loaded domain from " .. gridName)
number_elements=dom:domain_info():num_surface_elements()

-- This balancing setup makes sense for structured grids with uniform refinement
balancerDesc = {
	partitioner = {
		name = "dynamicBisection",
		verbose = false,
		enableXCuts = true,
		enableYCuts = true,
		enableZCuts = true,
		longestSplitAxis = false,
		clusteredSiblings = true,
		balanceThreshold = 0.9,
		numSplitImprovements = 10
	},


	hierarchy = {
		name                            = "standard",   --- ["standard", "lessRedists", "noRedists"]
		minElemsPerProcPerLevel         = 4,
		maxRedistProcs                  = 24576, 
		qualityRedistLevelOffset        = 1,
		intermediateRedistributions     = true,
		allowForRedistOnEachLevel       = true, -- DANGEROUS: 'true' only useful for non-adaptive runs with static partitioning. Leads to errors else.
		{-- level 0
		upperLvl = 0,
		maxProcs = 1
	},
	{-- levels 1
		upperLvl = 1,
		maxProcs = 192
	},
	{-- levels 2
		upperLvl = 2,
		maxProcs = 768
	},
	{-- levels 3
		upperLvl = 3,
		maxProcs = 6144
	},
	{-- levels 4
		upperLvl = 4,
		maxProcs = 24576
	},
	{-- levels 5
		upperLvl = 5,
		maxProcs = 32768
	}	
	},
}

util.refinement.CreateRegularHierarchy(dom, numRefs, true, balancerDesc)
--]]

print(dom:domain_info():to_string())

--FUNCTION TO CHECK PARALLEL STORAGE TYPE
function is_consistent(gf)
	if gf:has_storage_type_consistent() then
		print("GRID FUNCTION IS CONSISTENT")
		return true;
	end
	return false
end
function is_additive(gf)
	if gf:has_storage_type_additive() then
		print("GRID FUNCTION IS ADDITIVE")
		return true;
	end
	return false
end
function parallel_type_is(gf) 
	if is_consistent(gf) then
		return true;
	elseif is_additive(gf) then
		return true;
	else
		print("GRID FUNCTION IS UNDEFINED OR UNIQUE")
	end
	return false;
end
--NAVIER STOKES 
--DIRICHLET BOUNDARY VALUES AT INLET
function InletVelocities(x,y,z,t)
	local s=math.sqrt(y*y + z*z)*math.pi/diameter
	return math.max(0.0, math.cos(s))
	--return 0
end
function SquareInletVelocities(x,y,z,t)
	local s1=math.abs(y)*math.pi/diameter
	local s2=math.abs(z)*math.pi/diameter
	return math.max(0.0,math.cos(s1)*math.cos(s2))
	--return 0
end

print("NAVIER STOKES: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
NavierStokes_ApproxSpace = ApproximationSpace(dom)
NavierStokes_ApproxSpace:add_fct(flowNames,"Lagrange", vorder)
NavierStokes_ApproxSpace:add_fct(pressureName,"Lagrange",porder) 
NavierStokes_ApproxSpace:init_levels()
NavierStokes_ApproxSpace:init_top_surface()

print("NAVIER STOKES: Approx. Space:")
NavierStokes_ApproxSpace:print_statistic()

NavierStokes_ElemDisc = IncompressibleNavierStokes("v1,v2,v3,p ", "outer")
--useful settings
NavierStokes_ElemDisc:set_nonlinear(true)
NavierStokes_ElemDisc:set_picard(false)
NavierStokes_ElemDisc:set_kinematic_viscosity(visc)
NavierStokes_ElemDisc:set_stabilization(stab)
NavierStokes_ElemDisc:set_stabilization_type(stabType)
--************BOUNDARY CONDITIONS*********--
NavierStokes_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
NavierStokes_Dirich:add("InletVelocities","v1","inlet")
NavierStokes_Dirich:add(0,"v2","inlet")
NavierStokes_Dirich:add(0,"v3","inlet")
--************WALL BOUNDARY**************--
NavierStokes_Dirich:add(0,"v1","wall")
NavierStokes_Dirich:add(0,"v2","wall")
NavierStokes_Dirich:add(0,"v3","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
NavierStokes_Dirich:add(0,"v1","obstacle_surface")
NavierStokes_Dirich:add(0,"v2","obstacle_surface")
NavierStokes_Dirich:add(0,"v3","obstacle_surface")
if not bDoNothing then
	print("Output flows set")
	NavierStokes_Dirich:add("InletVelocities","v1","outlet")
	NavierStokes_Dirich:add(0,"v2","outlet")
	NavierStokes_Dirich:add(0,"v3","outlet")
end
-- Domain Discretization 
NavierStokes_DomainDisc = DomainDiscretization(NavierStokes_ApproxSpace)
NavierStokes_DomainDisc:add(NavierStokes_ElemDisc)
NavierStokes_DomainDisc:add(NavierStokes_Dirich)
print("NAVIER STOKES: create GridFunctions, Matrix Operator, and GlobalGridFunctionGradientDatas")
v = GridFunction(NavierStokes_ApproxSpace);v:set(0.0)
v_temp = GridFunction(NavierStokes_ApproxSpace);v_temp:set(0.0)
NavierStokes_DomainDisc:adjust_solution(v)
NavierStokes_DomainDisc:adjust_solution(v_temp)
v1_gradient_global=GlobalGridFunctionGradientData(v,"v1")
v2_gradient_global=GlobalGridFunctionGradientData(v,"v2")
v3_gradient_global=GlobalGridFunctionGradientData(v,"v3")
v1_value_global=GlobalGridFunctionNumberData(v,"v1")
v2_value_global=GlobalGridFunctionNumberData(v,"v2")
v3_value_global=GlobalGridFunctionNumberData(v,"v3")
---Pressure
p_value_global=GlobalGridFunctionNumberData(v,"p")
--Nonlinear Matrix Operator
navier_Op = AssembledOperator()
navier_Op:set_discretization(NavierStokes_DomainDisc)
print("NAVIER STOKES: all set")


print("ADJOINT FLOW: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
AdjointFlow_ApproxSpace = ApproximationSpace(dom)
AdjointFlow_ApproxSpace:add_fct(adjointFlowNames,"Lagrange",vorder)
AdjointFlow_ApproxSpace:add_fct(adjointPressure,"Lagrange",porder) 
AdjointFlow_ApproxSpace:init_levels()
AdjointFlow_ApproxSpace:init_top_surface()
print("ADJOINT FLOW : Approx. Space:")
AdjointFlow_ApproxSpace:print_statistic()

AdjointFlow_ElemDisc = NavierStokesAdjoint("q1,q2,q3, h ", "outer")
AdjointFlow_ElemDisc:set_kinematic_viscosity(visc)
AdjointFlow_ElemDisc:set_stabilization(stab)
AdjointFlow_ElemDisc:set_stabilization_type(stabType)
--Set Imports
--VELOCITY GRADIENT
AdjointFlow_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
AdjointFlow_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
AdjointFlow_ElemDisc:set_velocity_vector_d3(v3_gradient_global);

--VELOCITY VECTOR
AdjointFlow_ElemDisc:set_velocity_d1(v1_value_global)
AdjointFlow_ElemDisc:set_velocity_d2(v2_value_global)
AdjointFlow_ElemDisc:set_velocity_d3(v3_value_global)

--************BOUNDARY CONDITIONS*********--
AdjointFlow_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","inlet")
AdjointFlow_Dirich:add(0,"q2","inlet")
AdjointFlow_Dirich:add(0,"q3","inlet")
--************WALL BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","wall")
AdjointFlow_Dirich:add(0,"q2","wall")
AdjointFlow_Dirich:add(0,"q3","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","obstacle_surface")
AdjointFlow_Dirich:add(0,"q2","obstacle_surface")
AdjointFlow_Dirich:add(0,"q3","obstacle_surface")
--TODO:see if do nothing is ok, might be adjustable think about it...
--************OUTLET BOUNDARY**************--
--AdjointFlow_Dirich:add(0,"q1","outlet")
--AdjointFlow_Dirich:add(0,"q2","outlet")

-- Domain Discretization 
AdjointFlow_DomainDisc = DomainDiscretization(AdjointFlow_ApproxSpace)
AdjointFlow_DomainDisc:add(AdjointFlow_ElemDisc)
AdjointFlow_DomainDisc:add(AdjointFlow_Dirich)

print("ADJOINT FLOW: create GridFunctions and GlobalGridFunctionGradientDatas")
q = GridFunction(AdjointFlow_ApproxSpace);q:set(0.0)
r_q = GridFunction(AdjointFlow_ApproxSpace);r_q:set(0.0)
AdjointFlow_DomainDisc:adjust_solution(q)
q1_gradient_global=GlobalGridFunctionGradientData(q,"q1")
q2_gradient_global=GlobalGridFunctionGradientData(q,"q2")
q3_gradient_global=GlobalGridFunctionGradientData(q,"q3")
q1_value_global=GlobalGridFunctionNumberData(q,"q1")
q2_value_global=GlobalGridFunctionNumberData(q,"q2")
q3_value_global=GlobalGridFunctionNumberData(q,"q3")
---Pressure
h_value_global=GlobalGridFunctionNumberData(q,"h")
---Pressure
A_adjFlow = AssembledLinearOperator(AdjointFlow_DomainDisc)
print("ADJOINT FLOW SYSTEM: all set:")

print("PLAPLACE DERIVATIVE: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
PLaplaceDerivative_ApproxSpace = ApproximationSpace(dom)
PLaplaceDerivative_ApproxSpace:add_fct("u1,u2,u3","Lagrange",0)
PLaplaceDerivative_ApproxSpace:init_levels()
PLaplaceDerivative_ApproxSpace:init_top_surface()
print("PLAPLACE DERIVATIVE: Approx. Space:")
PLaplaceDerivative_ApproxSpace:print_statistic()
--Saving of grid functions
grid_positions = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace)
--Solutions
delta_u_p = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);delta_u_p:set(0.0)
u_p = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);u_p:set(0.0);
u_p_negative = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);u_p_negative:set(0.0)--to revert transformation
sigma = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);sigma:set(0.0);--solution of rhs of schur complement

--RHS
Lu = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);Lu:set(0.0);
MinusLu_BdeltaLambda = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);MinusLu_BdeltaLambda:set(0.0);
u_zeros = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);u_zeros:set(0.0)
u1_gradient_global=GlobalGridFunctionGradientData(u_p,"u1")
u2_gradient_global=GlobalGridFunctionGradientData(u_p,"u2")
u3_gradient_global=GlobalGridFunctionGradientData(u_p,"u3")
u1_value_global=GlobalGridFunctionNumberData(u_p,"u1")
u2_value_global=GlobalGridFunctionNumberData(u_p,"u2")
u3_value_global=GlobalGridFunctionNumberData(u_p,"u3")
--THIS IS THE A MATRIX OF THE BIG PROBLEM ON THE LHS

PLaplaceDerivative_ElemDisc = PLaplaceDerivative("u1,u2,u3", "outer");PLaplaceDerivative_ElemDisc:set_quad_order(1)
PLaplaceDerivative_ElemDisc:set_lambda_vol(lambda_vol)
PLaplaceDerivative_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
PLaplaceDerivative_ElemDisc:set_p(p_current)
PLaplaceDerivative_ElemDisc:set_epsilon(epsilon_nnl)
PLaplaceDerivative_ElemDisc:set_elliptic_epsilon(epsilon_elliptic)
PLaplaceDerivative_ElemDisc:set_step_length(step_length)
--Set Imports
--DEFORMATION VECTOR
PLaplaceDerivative_ElemDisc:set_deformation_d1(u1_value_global)
PLaplaceDerivative_ElemDisc:set_deformation_d2(u2_value_global)
PLaplaceDerivative_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
PLaplaceDerivative_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
PLaplaceDerivative_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
PLaplaceDerivative_ElemDisc:set_deformation_vector_d3(u3_gradient_global)

--RHS definitions
PLaplaceDerivativeRHS_ElemDisc = PLaplaceDerivativeRHS("u1,u2,u3", "outer");PLaplaceDerivativeRHS_ElemDisc:set_quad_order(1)
PLaplaceDerivativeRHS_ElemDisc:set_kinematic_viscosity(visc)

PLaplaceDerivativeRHS_ElemDisc:set_lambda_vol(lambda_vol)
PLaplaceDerivativeRHS_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
PLaplaceDerivativeRHS_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
PLaplaceDerivativeRHS_ElemDisc:set_p(p_current)
PLaplaceDerivativeRHS_ElemDisc:set_epsilon(epsilon_nnl)
PLaplaceDerivativeRHS_ElemDisc:set_step_length(step_length)
--Set Imports
--DEFORMATION VECTOR
PLaplaceDerivativeRHS_ElemDisc:set_deformation_d1(u1_value_global)
PLaplaceDerivativeRHS_ElemDisc:set_deformation_d2(u2_value_global)
PLaplaceDerivativeRHS_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
--VELOCITY GRADIENT
PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d1(v1_gradient_global)
PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d2(v2_gradient_global)
PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
--VELOCITY VECTOR
PLaplaceDerivativeRHS_ElemDisc:set_velocity_d1(v1_value_global)
PLaplaceDerivativeRHS_ElemDisc:set_velocity_d2(v2_value_global)
PLaplaceDerivativeRHS_ElemDisc:set_velocity_d3(v3_value_global)
--PRESSURE
PLaplaceDerivativeRHS_ElemDisc:set_pressure(p_value_global)

--ADJOINT VELOCITY VECTOR
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
--ADJOINT VELOCITY GRADIENT
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
--ADJOINT PRESSURE
PLaplaceDerivativeRHS_ElemDisc:set_adjoint_pressure(h_value_global)

--************BOUNDARY CONDITIONS*********--
PLaplaceDerivative_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
PLaplaceDerivative_Dirich:add(0,"u1","inlet")
PLaplaceDerivative_Dirich:add(0,"u2","inlet")
PLaplaceDerivative_Dirich:add(0,"u3","inlet")
--************WALL BOUNDARY**************--
PLaplaceDerivative_Dirich:add(0,"u1","wall")
PLaplaceDerivative_Dirich:add(0,"u2","wall")
PLaplaceDerivative_Dirich:add(0,"u3","wall")
--************OUTLET BOUNDARY**************--
PLaplaceDerivative_Dirich:add(0,"u1","outlet")
PLaplaceDerivative_Dirich:add(0,"u2","outlet")
PLaplaceDerivative_Dirich:add(0,"u3","outlet")
-- Domain Discretization 
PLaplaceDerivative_DomainDisc = DomainDiscretization(PLaplaceDerivative_ApproxSpace)
PLaplaceDerivative_DomainDisc:add(PLaplaceDerivative_ElemDisc)
PLaplaceDerivative_DomainDisc:add(PLaplaceDerivative_Dirich)
PLaplaceDerivative_DomainDisc:add(PLaplaceDerivativeRHS_ElemDisc)

print("PLAPLACE DERIVATIVE: create GridFunctions and GlobalGridFunctionGradientDatas")
PLaplaceDerivative_DomainDisc:adjust_solution(sigma)
PLaplaceDerivative_DomainDisc:adjust_solution(u_p)
A_u_Hessian = AssembledLinearOperator(PLaplaceDerivative_DomainDisc)

print("PLAPLACE DERIVATIVE: all set:")

print("WE MUST CREATE A DISCRETIZATION FOR EACH CONSTRAINT FOR ITS M-SOLVE")
--B-MATRIX VECTORS
Bvol=AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);Bvol:set(0.0);
Bx = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);Bx:set(0.0);
By = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);By:set(0.0);
Bz = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);Bz:set(0.0);

print("FOR M=1, Bvol")
--Volume derivative
BVolume_ElemDisc = VolumeConstraintSecondDerivative("u1,u2,u3", "outer");BVolume_ElemDisc:set_quad_order(1)
BVolume_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
BVolume_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
BVolume_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
BVolume_ElemDisc:set_deformation_d1(u1_value_global)
BVolume_ElemDisc:set_deformation_d2(u2_value_global)
BVolume_ElemDisc:set_deformation_d3(u3_value_global)
--Bvol domain discretization
Bvol_DomainDisc = DomainDiscretization(PLaplaceDerivative_ApproxSpace)
Bvol_DomainDisc:add(BVolume_ElemDisc)
Bvol_DomainDisc:add(PLaplaceDerivative_ElemDisc)
Bvol_DomainDisc:add(PLaplaceDerivative_Dirich)
Avol = AssembledLinearOperator(Bvol_DomainDisc)
t_vol = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);t_vol:set(0.0);
Bvol_DomainDisc:adjust_solution(t_vol)
print("FOR M=2, Bx")
--Barycenter on X direction
XBarycenter_ElemDisc = XBarycenterConstraintSecondDerivative("u1,u2,u3", "outer");XBarycenter_ElemDisc:set_quad_order(1)
XBarycenter_ElemDisc:set_index(1)
XBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
XBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
XBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
XBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
XBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
XBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
--Bx domain discretization
Bx_DomainDisc = DomainDiscretization(PLaplaceDerivative_ApproxSpace)
Bx_DomainDisc:add(XBarycenter_ElemDisc)
Bx_DomainDisc:add(PLaplaceDerivative_ElemDisc)
Bx_DomainDisc:add(PLaplaceDerivative_Dirich)
Ax = AssembledLinearOperator(Bx_DomainDisc)
t_x = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);t_x:set(0.0);
Bx_DomainDisc:adjust_solution(t_x)
print("FOR M=3, By")
--Barycenter on Y direction
YBarycenter_ElemDisc = XBarycenterConstraintSecondDerivative("u1,u2,u3", "outer");YBarycenter_ElemDisc:set_quad_order(1) 
YBarycenter_ElemDisc:set_index(2)
YBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
YBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
YBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
YBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
YBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
YBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
--By domain discretization
By_DomainDisc = DomainDiscretization(PLaplaceDerivative_ApproxSpace)
By_DomainDisc:add(YBarycenter_ElemDisc)
By_DomainDisc:add(PLaplaceDerivative_ElemDisc)
By_DomainDisc:add(PLaplaceDerivative_Dirich)
Ay = AssembledLinearOperator(By_DomainDisc)
t_y = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);t_y:set(0.0);
By_DomainDisc:adjust_solution(t_y)
print("FOR M=4, Bz")
--Barycenter on Z direction
ZBarycenter_ElemDisc = XBarycenterConstraintSecondDerivative("u1,u2,u3", "outer");ZBarycenter_ElemDisc:set_quad_order(1) 
ZBarycenter_ElemDisc:set_index(3)
ZBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
ZBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
ZBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
ZBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
ZBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
ZBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
--By domain discretization
Bz_DomainDisc = DomainDiscretization(PLaplaceDerivative_ApproxSpace)
Bz_DomainDisc:add(ZBarycenter_ElemDisc)
Bz_DomainDisc:add(PLaplaceDerivative_ElemDisc)
Bz_DomainDisc:add(PLaplaceDerivative_Dirich)
Az = AssembledLinearOperator(Bz_DomainDisc)
t_z = AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);t_z:set(0.0);
Bz_DomainDisc:adjust_solution(t_z)

print("LARGE PROBLEM DOMAIN DISCRETIZATION")
--Create a RHS element discretization
LargeProblemRHS_ElemDisc = LargeProblemRHS("u1,u2,u3", "outer");LargeProblemRHS_ElemDisc:set_quad_order(1)
LargeProblemRHS_ElemDisc:set_kinematic_viscosity(visc)
LargeProblemRHS_ElemDisc:set_lambda_vol(lambda_vol)
LargeProblemRHS_ElemDisc:set_lambda_barycenter(lambda_x,lambda_y,lambda_z)
LargeProblemRHS_ElemDisc:set_p(p_current)
LargeProblemRHS_ElemDisc:set_epsilon(epsilon_nnl)
LargeProblemRHS_ElemDisc:set_step_length(step_length)
--Set Imports
--DEFORMATION VECTOR
LargeProblemRHS_ElemDisc:set_deformation_d1(u1_value_global)
LargeProblemRHS_ElemDisc:set_deformation_d2(u2_value_global)
LargeProblemRHS_ElemDisc:set_deformation_d3(u3_value_global)
--DEFORMATION GRADIENT
LargeProblemRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
LargeProblemRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
LargeProblemRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
--VELOCITY GRADIENT
LargeProblemRHS_ElemDisc:set_velocity_vector_d1(v1_gradient_global)
LargeProblemRHS_ElemDisc:set_velocity_vector_d2(v2_gradient_global)
LargeProblemRHS_ElemDisc:set_velocity_vector_d3(v3_gradient_global)
--VELOCITY VECTOR
LargeProblemRHS_ElemDisc:set_velocity_d1(v1_value_global)
LargeProblemRHS_ElemDisc:set_velocity_d2(v2_value_global)
LargeProblemRHS_ElemDisc:set_velocity_d3(v3_value_global)
--PRESSURE
LargeProblemRHS_ElemDisc:set_pressure(p_value_global)
--ADJOINT VELOCITY VECTOR
LargeProblemRHS_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
LargeProblemRHS_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
LargeProblemRHS_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
--ADJOINT VELOCITY GRADIENT
LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
--ADJOINT PRESSURE
LargeProblemRHS_ElemDisc:set_adjoint_pressure(h_value_global)

--Large Problem Domain Discretization
LargeProblem_DomainDisc = DomainDiscretization(PLaplaceDerivative_ApproxSpace)
LargeProblem_DomainDisc:add(PLaplaceDerivative_ElemDisc)
LargeProblem_DomainDisc:add(LargeProblemRHS_ElemDisc)
LargeProblem_DomainDisc:add(PLaplaceDerivative_Dirich)
A_Large = AssembledLinearOperator(LargeProblem_DomainDisc)

--TODO:JPRIME DISCRETIZATION
	Jprime_ElemDisc=JPrime("u1,u2,u3","outer");Jprime_ElemDisc:set_quad_order(1)
	Jprime_ElemDisc:set_kinematic_viscosity(visc)
	--VELOCITY GRADIENT
	Jprime_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	Jprime_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
	Jprime_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	--VELOCITY VECTOR
	Jprime_ElemDisc:set_velocity_d1(v1_value_global)
	Jprime_ElemDisc:set_velocity_d2(v2_value_global)
	Jprime_ElemDisc:set_velocity_d3(v3_value_global)
	--PRESSURE
	Jprime_ElemDisc:set_pressure(p_value_global)
	--ADJOINT VELOCITY VECTOR
	Jprime_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
	--ADJOINT VELOCITY GRADIENT
	Jprime_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
	--ADJOINT PRESSURE
	Jprime_ElemDisc:set_adjoint_pressure(h_value_global)

	Jprime_DomainDisc=DomainDiscretization(PLaplaceDerivative_ApproxSpace)
	Jprime_DomainDisc:add(Jprime_ElemDisc)
	Sensitivity=AdvancedGridFunction(PLaplaceDerivative_ApproxSpace)

print("CREATION OF ELEMENT AND DOMAINDISCS FOR THE MATRIX B")
BTranspose_sigma = matrix(m,1)
L_lambda = matrix(m,1)--format from package lua-matrix
Lambda = matrix(m,1)--setting of initial values, default is 0
Lambda[1][1]=lambda_vol
Lambda[2][1]=lambda_x
Lambda[3][1]=lambda_y
Lambda[4][1]=lambda_z
S = matrix(m,m)
inverse_S = matrix(m,m)
rhs = matrix(m,1)
DeltaLambda=matrix(m,1)

ReferenceVolume = VolumeDefect(u_p,0, "outer", "u1,u2,u3",6,false,1,false)--used as a constant
vBarycenterDefect = BarycenterDefect(u_p,"u1,u2,u3","outer",6)
print("The reference volume is: "..ReferenceVolume)

print("CREATION OF SOLVERS")
NavierStokes_Solver=util.oo.ns_solver(NavierStokes_DomainDisc,NavierStokes_ApproxSpace)
--NavierStokes_Solver = CreateNSSolver(NavierStokes_ApproxSpace,"p","bicgstab",NavierStokes_DomainDisc)
AdjointFlow_Solver = util.oo.adjoint_ns_solver(AdjointFlow_DomainDisc,AdjointFlow_ApproxSpace)
--AdjointFlow_Solver = CreateAdjointFlowSolver(AdjointFlow_ApproxSpace, "h", "bicgstab", AdjointFlow_DomainDisc)
PLaplaceDerivative_Solver = util.oo.linear_solver(PLaplaceDerivative_DomainDisc,PLaplaceDerivative_ApproxSpace,false)
Bvol_Solver = util.oo.linear_solver(Bvol_DomainDisc,PLaplaceDerivative_ApproxSpace,false)
Bx_Solver = util.oo.linear_solver(Bx_DomainDisc,PLaplaceDerivative_ApproxSpace,false)
By_Solver = util.oo.linear_solver(By_DomainDisc,PLaplaceDerivative_ApproxSpace,false)
Bz_Solver = util.oo.linear_solver(Bz_DomainDisc,PLaplaceDerivative_ApproxSpace,false)
LargeProblem_Solver = util.oo.linear_solver(LargeProblem_DomainDisc,PLaplaceDerivative_ApproxSpace,true)

print("NAVIER STOKES SOLVER IS:")
print(NavierStokes_Solver:config_string())
print("ADJOINT FLOWS SOLVER IS:")
print(AdjointFlow_Solver:config_string())

vtkWriter = VTKOutput()
--Storage of data
vStep = {}
vDrag = {}
vShapeDerivative = {}
vSupNorm = {}
vStepSize = {}
vStepLength = {}
--Solver data storage per step
vTotalLinearIterations = {}
vLargeSolverIterations = {}
vBvolSolverIterations = {}
vBxSolverIterations = {}
vBySolverIterations = {}
vBzSolverIterations = {}
vRHSSolver = {}
vNonLinearIterationsPerP = {}
sum_total_iterations = 0
sum_largesolver_iterations = 0
sum_bvolsolver_iterations = 0
sum_bxsolver_iterations = 0
sum_bysolver_iterations = 0
sum_bzsolver_iterations = 0
sum_rhssolver_iterations = 0
sum_newtonsteps = 0



step = 0

if bReadInitialGuess > 0 then 
	print("RESTART ACTIVATED")
	print("STARTING FROM STEP:"..bReadInitialGuess)
	initial_nodal_positions = "nodal_positions_"..bReadInitialGuess
	ReadFromFile(grid_positions, initial_nodal_positions)
	SetDomainCoordinatesFromGF(grid_positions,"u1,u2,u3")
	print("LOADED NODAL COORDINATES FOR STEP"..bReadInitialGuess)
	--we load coordinates for the last calculated step, therefore we must start at restart + 1
	step = bReadInitialGuess+1; 
	
	ReferenceVolume = VolumeDefect(u_p,0, "outer", "u1,u2,u3",6,false,1,false)--used as a constant
	vBarycenterDefect = {0,0,0}
	print("The reference volume is: "..ReferenceVolume)

	u1_gradient_global=GlobalGridFunctionGradientData(u_p,"u1")
	u2_gradient_global=GlobalGridFunctionGradientData(u_p,"u2")
	u3_gradient_global=GlobalGridFunctionGradientData(u_p,"u3")
	u1_value_global=GlobalGridFunctionNumberData(u_p,"u1")
	u2_value_global=GlobalGridFunctionNumberData(u_p,"u2")
	u3_value_global=GlobalGridFunctionNumberData(u_p,"u3")

	--Reset the global data from Navier Stokes
	v1_gradient_global=GlobalGridFunctionGradientData(v,"v1")
	v2_gradient_global=GlobalGridFunctionGradientData(v,"v2")
	v3_gradient_global=GlobalGridFunctionGradientData(v,"v3")
	v1_value_global=GlobalGridFunctionNumberData(v,"v1")
	v2_value_global=GlobalGridFunctionNumberData(v,"v2")
	v3_value_global=GlobalGridFunctionNumberData(v,"v3")
	---Pressure
	p_value_global=GlobalGridFunctionNumberData(v,"p")

	--Reset the global data from Adjoint Navier Stokes
	q1_gradient_global=GlobalGridFunctionGradientData(q,"q1")
	q2_gradient_global=GlobalGridFunctionGradientData(q,"q2")
	q3_gradient_global=GlobalGridFunctionGradientData(q,"q3")
	q1_value_global=GlobalGridFunctionNumberData(q,"q1")
	q2_value_global=GlobalGridFunctionNumberData(q,"q2")
	q3_value_global=GlobalGridFunctionNumberData(q,"q3")
	---Pressure
	h_value_global=GlobalGridFunctionNumberData(q,"h")
	--Set anew in AdjointNavierStokes
	AdjointFlow_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	AdjointFlow_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
	AdjointFlow_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	AdjointFlow_ElemDisc:set_velocity_d1(v1_value_global)
	AdjointFlow_ElemDisc:set_velocity_d2(v2_value_global)
	AdjointFlow_ElemDisc:set_velocity_d3(v3_value_global)
	--Set anew in PLaplace Matrix for Large Problem
	PLaplaceDerivative_ElemDisc:set_deformation_d1(u1_value_global)
	PLaplaceDerivative_ElemDisc:set_deformation_d2(u2_value_global)
	PLaplaceDerivative_ElemDisc:set_deformation_d3(u3_value_global)
	PLaplaceDerivative_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	PLaplaceDerivative_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	PLaplaceDerivative_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	--Set anew in PLaplace RHS for Large Problem
	LargeProblemRHS_ElemDisc:set_deformation_d1(u1_value_global)
	LargeProblemRHS_ElemDisc:set_deformation_d2(u2_value_global)
	LargeProblemRHS_ElemDisc:set_deformation_d3(u3_value_global)
	LargeProblemRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	LargeProblemRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	LargeProblemRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	LargeProblemRHS_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	LargeProblemRHS_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
	LargeProblemRHS_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	LargeProblemRHS_ElemDisc:set_velocity_d1(v1_value_global)
	LargeProblemRHS_ElemDisc:set_velocity_d2(v2_value_global)
	LargeProblemRHS_ElemDisc:set_velocity_d3(v3_value_global)
	LargeProblemRHS_ElemDisc:set_pressure(p_value_global)
	LargeProblemRHS_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
	LargeProblemRHS_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
	LargeProblemRHS_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
	LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
	LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
	LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
	LargeProblemRHS_ElemDisc:set_adjoint_pressure(h_value_global)
	--Set anew in PLaplace RHS
	PLaplaceDerivativeRHS_ElemDisc:set_deformation_d1(u1_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_deformation_d2(u2_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_deformation_d3(u3_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
	PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	PLaplaceDerivativeRHS_ElemDisc:set_velocity_d1(v1_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_velocity_d2(v2_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_velocity_d3(v3_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_pressure(p_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
	PLaplaceDerivativeRHS_ElemDisc:set_adjoint_pressure(h_value_global)
	--Volume linear system
	BVolume_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	BVolume_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	BVolume_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	BVolume_ElemDisc:set_deformation_d1(u1_value_global)
	BVolume_ElemDisc:set_deformation_d2(u2_value_global)
	BVolume_ElemDisc:set_deformation_d3(u3_value_global)
	--Bx linear system
	XBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	XBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	XBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	XBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
	XBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
	XBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
	--By linear system
	YBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	YBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	YBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	YBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
	YBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
	YBarycenter_ElemDisc:set_deformation_d3(u3_value_global)	
	--Bz linear system
	ZBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
	ZBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
	ZBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
	ZBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
	ZBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
	ZBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
	--JPrime vector
    Jprime_ElemDisc:set_kinematic_viscosity(visc)
	Jprime_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
	Jprime_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
	Jprime_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
	Jprime_ElemDisc:set_velocity_d1(v1_value_global)
	Jprime_ElemDisc:set_velocity_d2(v2_value_global)
	Jprime_ElemDisc:set_velocity_d3(v3_value_global)
	Jprime_ElemDisc:set_pressure(p_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
	Jprime_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
	Jprime_ElemDisc:set_adjoint_pressure(h_value_global)
end
if bDebugNodalPositions then
	print("Nodal positions saved to vtk")
	vtkWriter:clear_selection()
	vtkWriter:select_nodal("u1,u2,u3","up")
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print("grid_positions", grid_positions,step,step, false)
end
p_solver_failure = false--make sure iteration starts
if bOutputPProblem then
	print("u_p final is saved to ..."..GenerateString(p_max))
	vtkWriter:clear_selection()
	vtkWriter:select_nodal("u1,u2,u3","up")
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print("up_final_"..GenerateString(p_max), u_p,step,step, false)
end
print("SOLVE PHASE: NON-LINEAR SOLUTION OF THE NAVIER STOKES PROBLEM")
NavierStokes_Solver:init(navier_Op)----TODO:navierLine
if NavierStokes_Solver:prepare(v) == false then print ("NavierStokes: Newton solver prepare failed at step "..step.."."); exit(); end 			
if NavierStokes_Solver:apply(v) == false then print ("NavierStokes: Newton solver failed at step "..step.."."); exit(); end 

VecScaleAssign(v_temp,1.0,v)--store as initial guess

drag_old = 0.5 * visc * Drag(u_zeros, v, "v1,v2,v3","u1,u2,u3","outer",3)
vDrag[step] = 0.5 * visc * Drag(u_zeros, v, "v1,v2,v3","u1,u2,u3","outer",3)

if bOutputFlows and (bReadInitialGuess == -1) then
	print("Flows are saved to '" .. flowOutputFile .. "'...")
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(flowNames, nodalFlows)
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print(flowOutputFile, v, step,step, false)
end
if bOutputPressure and bReadInitialGuess == -1 then
	--PRESSURE
	print("Pressure field is saved")
	vtkWriter:clear_selection()
	vtkWriter:select_nodal("p", "pressure")
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print("pressureOutputFile", v,step,step, false)
end

steps_without_line_search = 0

--Begin of optimization loop
print("+++++++++++++++++++ BEGIN OPTIMIZATION LOOP  +++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
while (step <= numSteps) do
	vStep[step] = step
	steps_without_line_search = steps_without_line_search + 1
	print("+++++++++++++++++++ NEW OPTIMIZATION STEP "..step.." +++++++++++++++++++")
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

	print("SOLVE PHASE: LINEAR SOLUTION OF THE ADJOINT FLOWS PROBLEM")
	AdjointFlow_DomainDisc:assemble_linear(A_adjFlow, r_q)--TODO:adjointLine
	AdjointFlow_Solver:init(A_adjFlow, q)
	if AdjointFlow_Solver:apply(q, r_q) == false then print("Adjoint flows solver failed, at step "..ns_i.." p "..p_current);exit(); end
	if bOutputAdjoints then
		print("Adjoints is saved to '" .. adjointFlowOutputFile .. "'...")
		vtkWriter:clear_selection()
		vtkWriter:select_nodal(adjointFlowNames, adjointNodalFlows)
		--vtkWriter:select_all(false)	-- write all variables
		vtkWriter:print(adjointFlowOutputFile, q,step,step, false)
		--PRESSURE
		vtkWriter:clear_selection()
		vtkWriter:select_nodal("h", "adjoint_pressure")
		--vtkWriter:select_all(false)	-- write all variables
		vtkWriter:print("adjointPressureOutputFile", q,step,step, false)
	end
	--TODO:J' gets assembled here 
    Jprime_DomainDisc:assemble_defect(Sensitivity,u_p)
    SetZeroAwayFromSubset(Sensitivity,"u1,u2,u3","obstacle_surface")
    --Sensitivity:change_storage_type_to_consistent()
    if bDebugSensitivity then
        print("Sensitivity is printed")
        vtkWriter:clear_selection()
        vtkWriter:select_nodal("u1,u2,u3", "J'")
        --vtkWriter:select_all(false)	-- write all variables
        vtkWriter:print("senstivity", Sensitivity, step,step, false)
    end	--variables store totals for each step, set to zero upon starting a new step
	sum_total_iterations = 0
	sum_largesolver_iterations = 0
	sum_bvolsolver_iterations = 0
	sum_bxsolver_iterations = 0
	sum_bysolver_iterations = 0
	sum_bzsolver_iterations = 0
	sum_rhssolver_iterations = 0
	sum_newtonsteps = 0
	--here we set stuff to zero for a new step

	while(true) do--step size control loop
		p_current = p_init
		--L_lambda is zero at start of new step
		L_lambda[1][1] = 0.0
		L_lambda[2][1] = 0.0
		L_lambda[3][1] = 0.0
		L_lambda[4][1] = 0.0
		 u_p:set(0.0);
		Lambda[1][1]=0.0;Lambda[2][1]=0.0;Lambda[3][1]=0.0;Lambda[4][1]=0.0;
		PLaplaceDerivative_ElemDisc:set_lambda_vol(0.0)
		PLaplaceDerivative_ElemDisc:set_lambda_barycenter(0.0,0.0,0.0)				
		PLaplaceDerivativeRHS_ElemDisc:set_lambda_vol(0.0)
		PLaplaceDerivativeRHS_ElemDisc:set_lambda_barycenter(0.0,0.0,0.0)		
		LargeProblemRHS_ElemDisc:set_lambda_vol(0.0)
		LargeProblemRHS_ElemDisc:set_lambda_barycenter(0.0,0.0,0.0)
		--Update p-value on the rhs vectors and matrices
		PLaplaceDerivative_ElemDisc:set_p(p_current)
		PLaplaceDerivativeRHS_ElemDisc:set_p(p_current);LargeProblemRHS_ElemDisc:set_p(p_current);
		
		if p_solver_failure == true then--step had failed, set bool back to to false
			p_solver_failure=false--restart boolean
			print("+++++++++++++++++++ AT OPTIMIZATION STEP "..step.." +++++++++++++++++")
			print("++++++++++++++ THE STEP IS REPEATED BECAUSE OF SOLVER FAILURE+++++++")
			print("+++++++++++++++THE ADJOINTS ARE NOT CALCULATED +++++++++++++++++++++")
		end
		while (true) do--start p loop
			PLaplaceDerivative_ElemDisc:set_control(math.pow(step_control,p_current-1.0));--This is A
			PLaplaceDerivativeRHS_ElemDisc:set_control(1.0/math.pow(step_control,p_current-1.0));--This is Lu (For small problem RHS)
			LargeProblemRHS_ElemDisc:set_control(1.0/math.pow(step_control,p_current-1.0));--This is Lu - B'delta_lambda (Big Problem RHS)
			PLaplaceDerivative_ElemDisc:set_step_length(math.pow(step_length, p_current-1.0))
			PLaplaceDerivativeRHS_ElemDisc:set_step_length(math.pow(step_length, p_current-1.0));LargeProblemRHS_ElemDisc:set_step_length(math.pow(step_length, p_current-1.0));
			--START NEWTON SOLVER METHOD
			ns_i=1
			--Store convergence of Newton's method and residual for the current step
			vNS_Step = {}
			vNS_NormDeltaUp = {}
			vNS_NormDeltaLambda= {}
			vNS_NormSum = {}
			vNS_Lu2Norm = {}
			vNS_LuEuclideanNorm = {}
			vNS_EucNormDeltaUp = {}
			vNS_EucNormLargeRHS= {}
			vNS_L_lambdaNorm = {}
			vNS_RHS_Solver_Iterations = {}
			vNS_LargeSolver_Iterations = {}
			vNS_BxSolver_Iterations = {}
			vNS_BySolver_Iterations = {}
			vNS_BzSolver_Iterations = {}
			vNS_BvolSolver_Iterations = {}
			--String for file generation
			p_as_string = GenerateString(p_current)
			
			print("******************START: NEWTON'S METHOD******************")
			print("+++++++++++++++OPTIMIZATION STEP "..step.."+++++++++++++++")
			ns_i = 1--set newton step counter

			while(ns_i <= nsMaxIts) do--nonlinear solver of relaxed p-problem
				--for convergence files
				vNS_Step[ns_i]=ns_i
				--restart certain data structures
				MinusLu_BdeltaLambda:set(0.0);Lu:set(0.0)
                delta_u_p:set(0.0);sigma:set(0.0);
				rhs[1][1]=0.0;rhs[2][1]=0.0;rhs[3][1]=0.0;rhs[4][1]=0.0;
				DeltaLambda[1][1]=0.0;DeltaLambda[2][1]=0.0;DeltaLambda[3][1]=0.0;DeltaLambda[4][1]=0.0;
				BTranspose_sigma[1][1]=0.0;BTranspose_sigma[2][1]=0.0;BTranspose_sigma[3][1]=0.0;BTranspose_sigma[4][1]=0.0;

				--at set(0.0) the gfs are all storage types, inclduing PST_CONSISTENT
				Bvol:set(0.0);Bx:set(0.0);By:set(0.0);Bz:set(0.0);
				--after assemble_defect the gfs are PST_ADDITIVE
				Bvol_DomainDisc:assemble_defect(Bvol,u_p);Bx_DomainDisc:assemble_defect(Bx,u_p);By_DomainDisc:assemble_defect(By,u_p);Bz_DomainDisc:assemble_defect(Bz,u_p)
				--these correct the minus sign introduced by the assemble_defect routine
				VecScaleAssign(Bvol,-1.0,Bvol)
				VecScaleAssign(Bx,-1.0,Bx)
				VecScaleAssign(By,-1.0,By)
				VecScaleAssign(Bz,-1.0,Bz)

				B_vector = {Bvol, Bx, By, Bz}
				print("****************** NEWTON'S METHOD ITERATION 	#:"..ns_i.."*****************************")
				print("****************** OPTIMIZATION STEP:"..step.." FOR P = "..p_as_string.."*******************")
				print("****************** CONTROL VALUE IS: "..step_control.." *******************")
				--We want to solve S=-L_lambda+B'inv(A)Lu
				--Calculate rhs = -L_lambda+BTranspose_sigma, solve A.sigma=Lu for sigma
				print("1.- SOLVE LINEAR PROBLEM OF RHS: A.sigma=(-Lu)")	
				PLaplaceDerivative_DomainDisc:adjust_solution(sigma)--PST_CONSISTENT
                PLaplaceDerivative_DomainDisc:assemble_jacobian(A_u_Hessian, u_p)
                PLaplaceDerivative_DomainDisc:assemble_defect(Lu,u_p)--Lu is additive
				VecScaleAdd2(Lu,1.0,Lu,1.0,Sensitivity)
                if Lu:has_storage_type_additive() == false then print("CATASTROPHIC FAILURE::RHS NOT ADDITIVE");exit();end
				PLaplaceDerivative_Solver:init(A_u_Hessian, sigma)
				if PLaplaceDerivative_Solver:apply(sigma, Lu) == false then print("A.sigma=Lu, solver failed, at step "..ns_i.." p "..p_as_string);p_solver_failure = true; break; end
                VecScaleAssign(sigma,-1.0,sigma)
				Lu:change_storage_type_to_consistent()--this ruins Lu
				if bDebugOutput then
					print("Lu is saved to ConsistentLu_step_"..step.."_p_"..p_as_string)
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2,u3","up")
					vtkWriter:print("ConsistentLu_step_"..step.."_p_"..p_as_string, Lu,ns_i,ns_i, false)
				end				
				--Perform B.delta_lambda
				--sigma:consistent while B_vector[i] is additive
				--in vector prod <left, right> left changes of storage
				for i = 1, m do
					BTranspose_sigma[i][1]=VecProd(B_vector[i], sigma);
				end
				rhs[1][1]= -(1.0/math.pow(step_control, p_current-1.0))*L_lambda[1][1] - BTranspose_sigma[1][1]
				rhs[2][1]= -(1.0/math.pow(step_control, p_current-1.0))*L_lambda[2][1] - BTranspose_sigma[2][1]
				rhs[3][1]= -(1.0/math.pow(step_control, p_current-1.0))*L_lambda[3][1] - BTranspose_sigma[3][1]
				rhs[4][1]= -(1.0/math.pow(step_control, p_current-1.0))*L_lambda[4][1] - BTranspose_sigma[4][1]
				print("THE RHS OF REDUCED PROBLEM IS:")
                rhs:print()

				Bvol:set(0.0)
                Bvol_DomainDisc:assemble_defect(Bvol,u_p)
				Bvol_DomainDisc:assemble_jacobian(Avol,u_p)
				Bvol_Solver:init(Avol, t_vol)
				if Bvol_Solver:apply(t_vol, Bvol) == false then print("Solver for B[1] failed");p_solver_failure = true;break; end
                VecScaleAssign(Bvol,-1.0,Bvol)	
				S[1][1]= VecProd(Bvol,t_vol)
				S[2][1]= VecProd(Bx,t_vol)
				S[3][1]= VecProd(By,t_vol)
				S[4][1]= VecProd(Bz,t_vol)
				
				--X-Barycenter Constraint
				Bx:set(0.0)
                Bx_DomainDisc:assemble_defect(Bx,u_p)
				Bx_DomainDisc:assemble_jacobian(Ax,u_p)			
				Bx_Solver:init(Ax, t_x)
				if Bx_Solver:apply(t_x, Bx) == false then print("Solver for B[2] failed");p_solver_failure = true; break; end	
                VecScaleAssign(Bx,-1.0,Bx)

				S[1][2]= VecProd(Bvol,t_x)
				S[2][2]= VecProd(Bx,t_x)
				S[3][2]= VecProd(By,t_x)
				S[4][2]= VecProd(Bz,t_x)

				--Y-Barycenter Constraint
				By:set(0.0)
                By_DomainDisc:assemble_defect(By,u_p)
				By_DomainDisc:assemble_jacobian(Ay,u_p)
				By_Solver:init(Ay, t_y)
				if By_Solver:apply(t_y, By) == false then print("Solver for B[3] failed");p_solver_failure = true; break; end	
                VecScaleAssign(By,-1.0,By)
				S[1][3]= VecProd(Bvol,t_y)
				S[2][3]= VecProd(Bx,t_y)
				S[3][3]= VecProd(By,t_y)
				S[4][3]= VecProd(Bz,t_y)

				--Z-Barycenter Constraint
				Bz:set(0.0)
                Bz_DomainDisc:assemble_defect(Bz,u_p)
				Bz_DomainDisc:assemble_jacobian(Az,u_p)
				Bz_Solver:init(Az, t_z)
				if Bz_Solver:apply(t_z, Bz) == false then print("Solver for B[3] failed");p_solver_failure = true; break; end	
                VecScaleAssign(Bz,-1.0,Bz)
				S[1][4]= VecProd(Bvol,t_z)
				S[2][4]= VecProd(Bx,t_z)
				S[3][4]= VecProd(By,t_z)
				S[4][4]= VecProd(Bz,t_z)

				print("THE SCHUR COMPLEMENT MATRIX S IS:")	
				S:print()
				inverse_S = S:invert()
				print("THE INVERSE OF SCHUR COMPLEMENT MATRIX S IS:")	
				inverse_S:print()
				print("<S,inv(S)> RETURNS IDENTITY MATRIX:")	
				Iden = matrix(m,m);Iden = S:mul(inverse_S);Iden:print();
				--solve reduced problem by direct inverse
				DeltaLambda=matrix(m,1)
				DeltaLambda=inverse_S:mul(rhs)
				print("DeltaLambda IS:")
				DeltaLambda:print()
				lhs = matrix(m,1)
				lhs= S:mul(DeltaLambda)
				diff = matrix(m,1)
				diff=lhs:sub(rhs)				
				print("REDUCED PROBLEM RESIDUAL AND EUC.NORM: ")
				diff:print()
				print("NORM OF REDUCED PROBLEM: "..math.sqrt(diff[1][1]*diff[1][1]+diff[2][1]*diff[2][1]+diff[3][1]*diff[3][1]+diff[4][1]*diff[4][1]))
				if (math.sqrt(diff[1][1]*diff[1][1]+diff[2][1]*diff[2][1]+diff[3][1]*diff[3][1]+diff[4][1]*diff[4][1]) > 1.0e-10) then print("Reduced problem failed with norm > 1.0e-10"); end; 
				--end of reduced problem
				--solve big problem	
				MinusLu_BdeltaLambda:set(0.0)
				LargeProblemRHS_ElemDisc:set_multiplier_vol(DeltaLambda[1][1])
				LargeProblemRHS_ElemDisc:set_multiplier_bx(DeltaLambda[2][1])
				LargeProblemRHS_ElemDisc:set_multiplier_by(DeltaLambda[3][1])
				LargeProblemRHS_ElemDisc:set_multiplier_bz(DeltaLambda[4][1])
		
				LargeProblem_DomainDisc:adjust_solution(delta_u_p)
				print("LARGE PROBLEM SOLVER AT OPTIMIZATION STEP "..step.." AND NS STEP#   "..ns_i.." ")
                LargeProblem_DomainDisc:assemble_jacobian(A_Large, u_p)--pLaplaceLine
                LargeProblem_DomainDisc:assemble_defect(MinusLu_BdeltaLambda,u_p)
                VecScaleAdd2(MinusLu_BdeltaLambda,1.0,MinusLu_BdeltaLambda,1.0,Sensitivity)--Adds J' to vector
				VecScaleAssign(MinusLu_BdeltaLambda, math.pow(step_control, p_current-1.0), MinusLu_BdeltaLambda)--scales whole vector by step size
				LargeProblem_Solver:init(A_Large, delta_u_p)
				if LargeProblem_Solver:apply_return_defect(delta_u_p, MinusLu_BdeltaLambda) == false then print("Large problem solver failed at step "..ns_i);p_solver_failure = true; break; end		
				MinusLu_BdeltaLambda:change_storage_type_to_consistent()--this ruins the vector
				if bDebugOutput then
					print("-Lu-B.delta_lambda is saved to RHSBigProb_"..step.."_p_"..p_as_string)
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2,u3","up")
					vtkWriter:print("RHSBigProb_"..step.."_p_"..p_as_string, Lu,ns_i,ns_i, false)
				end
				--cleanup of the multipliers
				LargeProblemRHS_ElemDisc:set_multiplier_vol(0.0)
				LargeProblemRHS_ElemDisc:set_multiplier_bx(0.0)
				LargeProblemRHS_ElemDisc:set_multiplier_by(0.0)
				LargeProblemRHS_ElemDisc:set_multiplier_bz(0.0)
				--UPDATE DEFORMATION ITERATE u_new=u_old+delta_u
				VecScaleAdd2(u_p,1.0,u_p,-1.0,delta_u_p)
				--UPDATE LAGR.MULTS. ITERATE lambda_new=lambda_old+lambda_u
				Lambda[1][1]= Lambda[1][1] + DeltaLambda[1][1]
				Lambda[2][1]= Lambda[2][1] + DeltaLambda[2][1] 
				Lambda[3][1]= Lambda[3][1] + DeltaLambda[3][1]
				Lambda[4][1]= Lambda[4][1] + DeltaLambda[4][1]	

				--counter increase and convergence check
				ns_i = ns_i+1
				if ns_i > nsMaxIts then
					print("++++++++++++++++ NEWTON METHOD DID NOT CONVERGE ++++++++++++++++++++++++++")
					print("++++++++++++++++ AT ITERATION   #: "..ns_i.." ++++++++++++++++++++++++++++")
					print("++++++++++++++++ AT P "..p_as_string.." AND OPT. STEP # : "..step.."++++++")
					p_solver_failure=true;break;--breaks from newton loop into p loop
				end

				--norm of the residual Lu
				lu_norm1=L2Norm(Lu,"u1",1,"outer")
				lu_norm2=L2Norm(Lu,"u2",1,"outer")
				lu_norm3=L2Norm(Lu,"u3",1,"outer")
				lu_norm_sum = math.sqrt(lu_norm1*lu_norm1+lu_norm2*lu_norm2+lu_norm3*lu_norm3)	
				lu_euc_norm=VecNorm(Lu)		
				--norm of the 0th delta_up iterate, starts at 0, could also be overriden
				delta_u_p_1=L2Norm(delta_u_p,"u1",1,"outer")
				delta_u_p_2=L2Norm(delta_u_p,"u2",1,"outer")
				delta_u_p_3=L2Norm(delta_u_p,"u3",1,"outer")
				delta_u_p_norm_sum=math.sqrt(delta_u_p_1*delta_u_p_1+delta_u_p_2*delta_u_p_2+delta_u_p_3*delta_u_p_3)
				delta_u_p_euc_norm=VecNorm(delta_u_p)				
				delta_lambda_norm=0.0
				for j=1,m do
					delta_lambda_norm= delta_lambda_norm+DeltaLambda[j][1]*DeltaLambda[j][1]
				end
				delta_lambda_norm = math.sqrt(delta_lambda_norm)
				delta_norms_sum= delta_lambda_norm+delta_u_p_norm_sum
				
				vNS_NormDeltaUp[ns_i-1] = delta_u_p_norm_sum
				vNS_NormDeltaLambda[ns_i-1] = delta_lambda_norm
				vNS_NormSum[ns_i-1] = delta_norms_sum
				vNS_Lu2Norm[ns_i-1]=lu_norm_sum;
				vNS_EucNormDeltaUp[ns_i-1]=delta_u_p_euc_norm
				vNS_LuEuclideanNorm[ns_i-1]=lu_euc_norm
				vNS_EucNormLargeRHS[ns_i-1]=VecNorm(MinusLu_BdeltaLambda)
				--for iteration counts for every p in every step
				vNS_RHS_Solver_Iterations[ns_i-1]= PLaplaceDerivative_Solver:step();
				vNS_LargeSolver_Iterations[ns_i-1]= LargeProblem_Solver:step()
				vNS_BxSolver_Iterations[ns_i-1] = Bx_Solver:step()
				vNS_BySolver_Iterations[ns_i-1] = By_Solver:step()
				vNS_BzSolver_Iterations[ns_i-1] = Bz_Solver:step()
				vNS_BvolSolver_Iterations[ns_i-1] = Bvol_Solver:step()
				--compute the geometrical constraints, assign to vector, negative is used above
				L_lambda[1][1]=VolumeDefect(u_p, ReferenceVolume, "outer", "u1,u2,u3",6,false,1,false)
				vBarycenterDefect = BarycenterDefect(u_p,"u1,u2,u3","outer",6)
				L_lambda[2][1]= vBarycenterDefect[1]
				L_lambda[3][1]= vBarycenterDefect[2]
				L_lambda[4][1]= vBarycenterDefect[3]

				L_lambda_norm=0.0
				for j=1,m do
					L_lambda_norm= L_lambda_norm+L_lambda[j][1]*L_lambda[j][1]
				end
				L_lambda_norm = math.sqrt(L_lambda_norm)
				vNS_L_lambdaNorm[ns_i-1]=L_lambda_norm
				print("GEOMETRICAL CONSTRAINTS DEFECTS (L_LAMBDA) AT NS STEP#   "..ns_i.."  ARE:")	
				L_lambda:print()
				
				--prepare for next iteration
				--Update Lagrange multipliers of the RHS of A.sigma=Lu, all A matrices, and RHS of Large problem
				PLaplaceDerivative_ElemDisc:set_lambda_vol(Lambda[1][1])
				PLaplaceDerivative_ElemDisc:set_lambda_barycenter(Lambda[2][1],Lambda[3][1],Lambda[4][1])
				
				PLaplaceDerivativeRHS_ElemDisc:set_lambda_vol(Lambda[1][1])
				PLaplaceDerivativeRHS_ElemDisc:set_lambda_barycenter(Lambda[2][1],Lambda[3][1],Lambda[4][1])
				
				LargeProblemRHS_ElemDisc:set_lambda_vol(Lambda[1][1])
				LargeProblemRHS_ElemDisc:set_lambda_barycenter(Lambda[2][1],Lambda[3][1],Lambda[4][1])
				
				--iteration counts update
				sum_rhssolver_iterations = sum_rhssolver_iterations + PLaplaceDerivative_Solver:step()
				sum_largesolver_iterations = sum_largesolver_iterations + LargeProblem_Solver:step()
				sum_bvolsolver_iterations = sum_bvolsolver_iterations + Bvol_Solver:step()
				sum_bxsolver_iterations = sum_bxsolver_iterations + Bx_Solver:step()
				sum_bysolver_iterations = sum_bysolver_iterations + By_Solver:step()
				sum_bzsolver_iterations = sum_bzsolver_iterations + Bz_Solver:step()

				print("#   "..ns_i.." INCREMENT NORM IS:    "..delta_norms_sum)
				print("#   "..ns_i.." DELTA_U INCREMENT NORM IS:    "..delta_u_p_norm_sum.."   DELTA_LAMBDA INCREMENT NORM IS:   "..delta_lambda_norm)
				if delta_lambda_norm <= nsTol then
					break;--break from newton into p loop
				end

			end--nonlinear solver of p_current problem
			if(p_solver_failure == true) then break; end
			sum_newtonsteps = sum_newtonsteps + ns_i


			if bNewtonOutput then 
				gnuplot.write_data("__NewtonStats_step_"..step.."_p_"..p_as_string.."_.txt", {vNS_Step,vNS_NormSum,vNS_NormDeltaUp,
																													vNS_NormDeltaLambda,vNS_Lu2Norm,
																													vNS_EucNormDeltaUp,vNS_LuEuclideanNorm,
																													vNS_EucNormLargeRHS,vNS_L_lambdaNorm}) 
				gnuplot.write_data("__NewtonIterations_step_"..step.."_p_"..p_as_string.."_.txt", {vNS_Step,vNS_RHS_Solver_Iterations,
									 			vNS_BvolSolver_Iterations,vNS_BxSolver_Iterations,vNS_BySolver_Iterations,vNS_BzSolver_Iterations,vNS_LargeSolver_Iterations}) 
			end  
            
			print("		++++++++++++++++ NEWTON METHOD CONVERGED ++++++++++++++++++++++")
			print("		++++++++++++++++ FOR p = "..p_as_string.." ++++++++++++++++++++")
			print("		++++++++++++++++ AT OPTIMIZATION STEP "..step.." ++++++++++++++")
			print("		++++++++++++++++ AFTER "..ns_i.." ITERATIONS ++++++++++++++++++")
			
			if(p_current == p_max)  then break; end
			
			p_current = p_current + p_increase;

			if(p_current > p_max) then p_current = p_max; end
				
			if bOutputIntermediateUp then
				print("u_p is saved to "..p_as_string)
				vtkWriter:clear_selection()
				vtkWriter:select_nodal("u1,u2,u3","up")
				--vtkWriter:select_all(false)	-- write all variables
				vtkWriter:print("up_intermediate_"..p_as_string, u_p,step,step, false)
			end
			--Update p-value on the rhs vectors and matrices
			PLaplaceDerivative_ElemDisc:set_p(p_current);
			PLaplaceDerivativeRHS_ElemDisc:set_p(p_current);LargeProblemRHS_ElemDisc:set_p(p_current);

		end--end p-loop
		--immediately check for solver failure
		if p_solver_failure == true then
			steps_without_line_search = 0

			Lambda[1][1]=0.0
			Lambda[2][1]=0.0
			Lambda[3][1]=0.0
			Lambda[4][1]=0.0
			step_control = 0.5*step_control

			--discard current iteration count
			sum_total_iterations = 0
			sum_largesolver_iterations = 0
			sum_bvolsolver_iterations = 0
			sum_bxsolver_iterations = 0
			sum_bysolver_iterations = 0
			sum_bzsolver_iterations = 0
			sum_rhssolver_iterations = 0
			sum_newtonsteps = 0			
		else
			--Update geometry
			TransformDomainByDisplacement(u_p, "u1,u2,u3")

			--Solve primal problem on temporary vector and check if descent direction
			NavierStokes_Solver:init(navier_Op)----TODO:navierLine2
			if NavierStokes_Solver:prepare(v_temp) == false then print("NavierStokes:");exit();end
			if NavierStokes_Solver:apply(v_temp) == false then print("NavierStokes:");exit();end
			
			--Compute objective function on new geometry
			local drag_current = 0.5 * visc * Drag(u_zeros,v_temp,"v1,v2,v3","u1,u2,u3","outer",3)
			ShapeDerivative = VecProd(Sensitivity,u_p)--here J_prime is additive and u_p is consistent
			vShapeDerivative[step]=ShapeDerivative
			--CALCULATION OF: |J_prime(Omega)nu_u+lambda.g_u(Omega)nu_u|
			--We make J_prime, Bvol, Bx, By consistent (and they become useless)
			Sensitivity:change_storage_type_to_consistent();
			Bvol:change_storage_type_to_consistent();Bx:change_storage_type_to_consistent();By:change_storage_type_to_consistent();Bz:change_storage_type_to_consistent()
			--Then we use VecScaleAdd3 to save to a GF operation lambda.g_u(Omega)nu_u
			G_u= AdvancedGridFunction(PLaplaceDerivative_ApproxSpace);G_u:set(0.0);
			--We store the result in J_prime, either ways we will set it to zero
			VecScaleAdd3(G_u,Lambda[1][1],Bvol,Lambda[2][1],Bx,Lambda[3][1],By);
			VecScaleAdd3(Sensitivity, 1.0, Sensitivity, Lambda[4][1],Bz,1.0, G_u);
			--SupNorm is calculated using the VecMaxNorm function
			SupNorm = VecMaxNorm(Sensitivity);
			vSupNorm[step]=SupNorm;
			gnuplot.write_data("__ParametersPerStep.txt", {vStep,vShapeDerivative,vSupNorm},false)
			Sensitivity:set(0.0);
			print("SHAPE DERIVATIVE VALUE:"..ShapeDerivative)
			print("SUPREMUM NORM VALUE:"..SupNorm)
			--check if drag is descent direction
			print("NEW DRAG: "..drag_current)
			print("OLD DRAG: "..drag_old)
			print("NEW -OLD DRAG: "..drag_current - drag_old)
			if ((drag_current-drag_old) > 0.0
				or (drag_current - drag_old > lineSearchParam*ShapeDerivative))
				and false
				then
				print("STEP: "..step.." NOT A DESCENT DIRECTION OR SOLVER FAILURE")
				VecScaleAssign(u_p_negative,-1.0,u_p);
				TransformDomainByDisplacement(u_p_negative, "u1,u2,u3")
				steps_without_line_search = 0
				step_control = 0.5*step_control			
				Lambda[1][1]=0.0
				Lambda[2][1]=0.0
				Lambda[3][1]=0.0
				Lambda[4][1]=0.0
				--discard current iteration count
				sum_total_iterations = 0
				sum_largesolver_iterations = 0
				sum_bvolsolver_iterations = 0
				sum_bxsolver_iterations = 0
				sum_bysolver_iterations = 0
				sum_bzsolver_iterations = 0
				sum_rhssolver_iterations = 0
				sum_newtonsteps = 0
			else
				if steps_without_line_search > 10 then
					step_control = math.min(2.0*step_control, step_control_init)
				end

				print("STEP: "..step.." IS A DESCENT DIRECTION") 		
				VecScaleAssign(v,1.0,v_temp)--store temporary Navier-Stokes solutions			
				vDrag[step+1] = drag_current;--save new reference			
				gnuplot.write_data("__Drag.txt", {vStep,vDrag},false)			
				if bOutputFlows then
					print("Flows are saved to '" .. flowOutputFile .. "'...")
					vtkWriter:clear_selection()
					vtkWriter:select_nodal(flowNames, nodalFlows)
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print(flowOutputFile, v, step+1,step+1, false)
				end
				if bOutputPressure then
					print("Pressure field is saved")
					--PRESSURE
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("p", "pressure")
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print("pressureOutputFile", v,step,step, false)
				end
				grid_positions:set(0.0);SaveNodalPositions2GridFunction(grid_positions,"u1,u2,u3")
				SaveToFile(grid_positions,"nodal_positions_"..step+1)
				if bDebugNodalPositions then
					print("Nodal positions saved to vtk")
					vtkWriter:clear_selection()
					vtkWriter:select_nodal("u1,u2,u3","up")
					--vtkWriter:select_all(false)	-- write all variables
					vtkWriter:print("grid_positions", grid_positions,step+1,step+1, false)
				end
				if bOutputPProblem then
				print("u_p final is saved to ..."..p_as_string)
				vtkWriter:clear_selection()
				vtkWriter:select_nodal("u1,u2,u3","up")
				--vtkWriter:select_all(false)	-- write all variables
				vtkWriter:print("up_final_"..p_as_string, u_p,step+1,step+1, false)
				end
				vTotalLinearIterations[step] = sum_largesolver_iterations + sum_bvolsolver_iterations + sum_bxsolver_iterations + sum_bysolver_iterations 
											 + sum_bzsolver_iterations + sum_rhssolver_iterations
				vBvolSolverIterations[step] = sum_bvolsolver_iterations
				vBxSolverIterations[step] = sum_bxsolver_iterations
				vBySolverIterations[step] = sum_bysolver_iterations
				vBzSolverIterations[step] = sum_bzsolver_iterations
				vRHSSolver[step] = sum_rhssolver_iterations
				vNonLinearIterationsPerP[step] = sum_newtonsteps
				gnuplot.write_data("__Iterations_per_step.txt", {vStep, vNonLinearIterationsPerP, vTotalLinearIterations, vRHSSolver, 
											vBvolSolverIterations, vBxSolverIterations, vBySolverIterations, vBzSolverIterations, vLargeSolverIterations})
				drag_old = drag_current
				vStepSize[step] = step_control
				vStepLength[step] = step_length
				gnuplot.write_data("__ControlSize.txt", {vStep, vStepSize, vStepLength},false)
				break;
			end
		end--if solver failure then simply to next iteration
	end--end of step size control loop
	if p_solver_failure == false then

		step = step + 1
		u1_gradient_global=GlobalGridFunctionGradientData(u_p,"u1")
		u2_gradient_global=GlobalGridFunctionGradientData(u_p,"u2")
		u3_gradient_global=GlobalGridFunctionGradientData(u_p,"u3")
		u1_value_global=GlobalGridFunctionNumberData(u_p,"u1")
		u2_value_global=GlobalGridFunctionNumberData(u_p,"u2")
		u3_value_global=GlobalGridFunctionNumberData(u_p,"u3")
		--Reset the global data from Navier Stokes
		v1_gradient_global=GlobalGridFunctionGradientData(v,"v1")
		v2_gradient_global=GlobalGridFunctionGradientData(v,"v2")
		v3_gradient_global=GlobalGridFunctionGradientData(v,"v3")
		v1_value_global=GlobalGridFunctionNumberData(v,"v1")
		v2_value_global=GlobalGridFunctionNumberData(v,"v2")
		v3_value_global=GlobalGridFunctionNumberData(v,"v3")
		---Pressure
		p_value_global=GlobalGridFunctionNumberData(v,"p")
		--Reset the global data from Adjoint Navier Stokes
		q1_gradient_global=GlobalGridFunctionGradientData(q,"q1")
		q2_gradient_global=GlobalGridFunctionGradientData(q,"q2")
		q3_gradient_global=GlobalGridFunctionGradientData(q,"q3")
		q1_value_global=GlobalGridFunctionNumberData(q,"q1")
		q2_value_global=GlobalGridFunctionNumberData(q,"q2")
		q3_value_global=GlobalGridFunctionNumberData(q,"q3")
		---Pressure
		h_value_global=GlobalGridFunctionNumberData(q,"h")

		--Set anew in AdjointNavierStokes
		AdjointFlow_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
		AdjointFlow_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
		AdjointFlow_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
		AdjointFlow_ElemDisc:set_velocity_d1(v1_value_global)
		AdjointFlow_ElemDisc:set_velocity_d2(v2_value_global)
		AdjointFlow_ElemDisc:set_velocity_d3(v3_value_global)

		--Set anew in PLaplace Matrix for Large Problem
		PLaplaceDerivative_ElemDisc:set_deformation_d1(u1_value_global)
		PLaplaceDerivative_ElemDisc:set_deformation_d2(u2_value_global)
		PLaplaceDerivative_ElemDisc:set_deformation_d3(u3_value_global)
		PLaplaceDerivative_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		PLaplaceDerivative_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		PLaplaceDerivative_ElemDisc:set_deformation_vector_d3(u3_gradient_global)

		--Set anew in PLaplace RHS for Large Problem
		LargeProblemRHS_ElemDisc:set_deformation_d1(u1_value_global)
		LargeProblemRHS_ElemDisc:set_deformation_d2(u2_value_global)
		LargeProblemRHS_ElemDisc:set_deformation_d3(u3_value_global)
		LargeProblemRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		LargeProblemRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		LargeProblemRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
		LargeProblemRHS_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
		LargeProblemRHS_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
		LargeProblemRHS_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
		LargeProblemRHS_ElemDisc:set_velocity_d1(v1_value_global)
		LargeProblemRHS_ElemDisc:set_velocity_d2(v2_value_global)
		LargeProblemRHS_ElemDisc:set_velocity_d3(v3_value_global)
		LargeProblemRHS_ElemDisc:set_pressure(p_value_global)
		LargeProblemRHS_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
		LargeProblemRHS_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
		LargeProblemRHS_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
		LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
		LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
		LargeProblemRHS_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
		LargeProblemRHS_ElemDisc:set_adjoint_pressure(h_value_global)

		--Set anew in PLaplace RHS
		PLaplaceDerivativeRHS_ElemDisc:set_deformation_d1(u1_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_deformation_d2(u2_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_deformation_d3(u3_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d1(v1_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d2(v2_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_velocity_vector_d3(v3_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_velocity_d1(v1_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_velocity_d2(v2_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_velocity_d3(v3_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_pressure(p_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
		PLaplaceDerivativeRHS_ElemDisc:set_adjoint_pressure(h_value_global)

		--Volume linear system
		BVolume_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		BVolume_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		BVolume_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
		BVolume_ElemDisc:set_deformation_d1(u1_value_global)
		BVolume_ElemDisc:set_deformation_d2(u2_value_global)
		BVolume_ElemDisc:set_deformation_d3(u3_value_global)

		--Bx linear system
		XBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		XBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		XBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
		XBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
		XBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
		XBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
		--By linear system
		YBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		YBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		YBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
		YBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
		YBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
		YBarycenter_ElemDisc:set_deformation_d3(u3_value_global)
		
		--Bz linear system
		ZBarycenter_ElemDisc:set_deformation_vector_d1(u1_gradient_global)
		ZBarycenter_ElemDisc:set_deformation_vector_d2(u2_gradient_global)
		ZBarycenter_ElemDisc:set_deformation_vector_d3(u3_gradient_global)
		ZBarycenter_ElemDisc:set_deformation_d1(u1_value_global)
		ZBarycenter_ElemDisc:set_deformation_d2(u2_value_global)
		ZBarycenter_ElemDisc:set_deformation_d3(u3_value_global)

    	--Jprime discretization    
        Jprime_ElemDisc:set_kinematic_viscosity(visc)
        --VELOCITY GRADIENT
        Jprime_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
        Jprime_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
        Jprime_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
        --VELOCITY VECTOR
        Jprime_ElemDisc:set_velocity_d1(v1_value_global)
        Jprime_ElemDisc:set_velocity_d2(v2_value_global)
        Jprime_ElemDisc:set_velocity_d3(v3_value_global)
        --PRESSURE
        Jprime_ElemDisc:set_pressure(p_value_global)
        --ADJOINT VELOCITY VECTOR
        Jprime_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
        Jprime_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
        Jprime_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
        --ADJOINT VELOCITY GRADIENT
        Jprime_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
        Jprime_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
        Jprime_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
        --ADJOINT PRESSURE
        Jprime_ElemDisc:set_adjoint_pressure(h_value_global)
	end--performed only if solver was succesful else goes back to beginning
end--end optimization loop







