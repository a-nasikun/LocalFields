#include "VectorFields.h"
#include <igl/unproject_onto_mesh.h>

#include "TestSolver.h"

int eigToShow = 0, basisId=0, numSample=1000, selectedVertex;

int main(int argc, char *argv[])
{
	bool selectFace = false;
	Eigen::MatrixXd C;


	VectorFields vectorFields;
	igl::opengl::glfw::Viewer		viewer;
	Eigen::MatrixXd					V;
	Eigen::MatrixXi					F, E;

	/* READING DATA */
	//string meshFile = "../Models/AIM177_Vase-Lion/177_vase-lion.off";
	string meshFile = "../Models/AIM894_Chinese Dragon/894_dragon_tris.obj";
	//string meshFile = "../Models/AIM894_Chinese Dragon/dragon_2000.obj";
	//string meshFile = "../Models/Cube/Cube_1400.obj";
	//string meshFile = "../Models/Cube/Cube_488.obj";
	//string meshFile = "../Models/Armadillo/Armadillo_1083.obj";
	//string meshFile = "../Models/Armadillo/Armadillo_10812.obj";
	//string meshFile = "../Models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "../Models/Plane/squarePlane_regular_784.obj";
	//string meshFile = "../Models/Plane/squared_small.obj";
	//string meshFile = "../Models/Sphere/sphere_small.obj";
	//string meshFile = "../Models/Sphere/sphere_half_1321.obj";
	//string meshFile = "../Models/Sphere/round_sphere_small.obj";
	//string meshFile = "../Models/Sphere/round_sphere_1500.obj";
	//string meshFile = "../Models/Sphere/round_sphere_10242.obj";
	//string meshFile = "../Models/Sphere/sphere_tiny.obj";
	//string meshFile = "../Models/Thorus/Thorus_small.obj";
	//string meshFile = "../Models/Thorus/Thorus_2304.obj";
	//string meshFile = "../Models/Monkey/Monkey_967.obj";

	/* ========================= PRE-PROCESS ==============================*/
	cout << "========================= PRE-PROCESS ==============================\n"; 
	vectorFields.readMesh(meshFile);
	vectorFields.computeEdges();
	vectorFields.computeAverageEdgeLength();
	vectorFields.computeFaceCenter();
	vectorFields.computeFaceNormal();
	vectorFields.constructVertexAdjacencyMatrix();
	vectorFields.constructVFNeighbors();
	vectorFields.constructVFNeighborsFull();
	vectorFields.constructFaceAdjacency3NMatrix();
	vectorFields.constructFaceAdjacency2RingMatrix();
	//vectorFields.constructFaceAdjacencyMatrix_IGL();
	//vectorFields.testAdjacencyAndEdges();

	//vectorFields.testAdjMV();
	
	vectorFields.getVF(V, F);
	viewer.data().set_mesh(V, F);
	C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);

	// ID SELECTION
	// plane, regular => 1130
	// Armadillo, 1083 => 542, 347, (467), (569)

	/* MATRIX CONSTRUCTIONS */
	srand(time(NULL));
	//int randVId = rand() % V.rows(); 
	int randFId = rand() % F.rows();
	vectorFields.constructNeigborRings(0);
	vectorFields.constructMassMatrices();

	//vectorFields.constructRotationMatrix();
	vectorFields.constructMappingMatrix();
	//vectorFields.constructGradient3D();
	//vectorFields.constructGradient2D();
	//vectorFields.constructArbitraryField();

	/*========== Test for GRADIENT and CURL-related functions ==========*/
	//vectorFields.computeGradArbField3D();
	//vectorFields.computeGradArbField2D();
	//vectorFields.computeCoGradArbField2D();
	//vectorFields.computeCoGradArbField3D();
	//vectorFields.computeCurl2D();
	//vectorFields.computeCurl3D();
	//vectorFields.computeDivergent3D();
	//vectorFields.computeDivergent2D();
	//vectorFields.computeCurlGradArbField3D();
	//vectorFields.computeCurlGradArbField2D();
	//vectorFields.computeCurlCoGradArbField3D();
	//vectorFields.computeCurlCoGradArbField2D();
	//vectorFields.computeDivGradArbField3D();
	//vectorFields.computeDivGradArbField2D();
	//vectorFields.computeDivCoGradArbField3D();
	//vectorFields.computeDivCoGradArbField2D();

	//vectorFields.constructSingularities();
	
	/* =========== Test on PROBLEM SOLVING-related functionalities ================*/
	vectorFields.constructStiffnessMatrices();
	vectorFields.constructMatrixB();
	//vectorFields.checkB2DStructure();

	/* ====================== GLOBAL PROBLEM ====================*/
	cout << "\n========================= GLOBAL PROBLEM =============================\n";
	vectorFields.setupGlobalProblem();


	/* ====================== LOCAL ELEMENTS ====================*/
	cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
	vectorFields.constructSamples(numSample);
	vectorFields.constructBasis();
	vectorFields.setAndSolveUserSystem();

	/* =========== Test on EIGENFIELDS-related funcionalities ============= */
	//vectorFields.testCurlEnergy();
	//vectorFields.constructLaplace2D();
	//vectorFields.computeEigenLaplace2D();
	//vectorFields.computeEigenLaplace3D();
	//vectorFields.computeEigenstructureGradient3D();
	//vectorFields.constructLaplace3D();

	//vectorFields.constructFaceAdjacencyMatrix();
	//vectorFields.constructNeigborRings(randFId);
	//vectorFields.constructRotationMatrix();
	//vectorFields.testDijkstraFace();

	/* ==================== VISUALIZATION ======================== */
	/* GLOBAL  */
	//vectorFields.visualizeApproximatedFields(viewer);
	//vectorFields.visualizeGlobalConstraints(viewer);
	//vectorFields.visualizeSingularitiesConstraints(viewer);

	/* LOCAL  */
	vectorFields.visualizeApproxResult(viewer, 0);	
	vectorFields.visualizeUserConstraints(viewer);
	//vectorFields.visualizeSamples(viewer);
	vectorFields.visualizeSingularitiesConstraints(viewer);
	

	/* FOR TESTING PURPOSE */
	//vectorFields.visualizeNeighboringRings(viewer);
	//vectorFields.visualizeDijkstraFace(viewer);
	//vectorFields.visualizeArbField(viewer);
	//vectorFields.visualizeEigenfields(viewer);
	//vectorFields.visualizeVertexFacesNeighbors(viewer, 0);

	/* MEASURE ACCURACY */
	//vectorFields.measureApproxAccuracyL2Norm();

	/* TEST SOLVERS */
	//testLAPACKEdsysv();
	//testCUDA_LULinearSolver();
	//testLDLTLinearSolver();
	//testMKL_Pardiso();

	const auto &key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		//int selectedFace;
		srand(time(NULL));
		selectedVertex = rand() % V.rows();

		switch (key)
		{
		case '1':
			//vectorFields.visualizeApproxResult(viewer, 0);
			//eigToShow = max(eigToShow - 1, 0);
			//vectorFields.visualizeEigFieldsDiv(viewer, eigToShow);
			break;
		case '2':
			//vectorFields.visualizeApproxResult(viewer, 1);
			//eigToShow = min(eigToShow + 1, 100);
			//vectorFields.visualizeEigFieldsDiv(viewer, eigToShow);
			break;
		case '3':
			basisId = max(basisId - 1, 0);
			vectorFields.visualizeBasis(viewer, basisId);
			break;
		case '4':
			basisId = min(basisId + 1, 2*numSample-1);
			vectorFields.visualizeBasis(viewer, basisId);
			break;
		case '5':
			basisId = max(basisId - 1, 0);
			vectorFields.visualizeBasisNormalized(viewer, basisId);
			break;
		case '6':
			basisId = min(basisId + 1, 2 * numSample - 1);
			vectorFields.visualizeBasisNormalized(viewer, basisId);
			break;
		case '7':
			vectorFields.visualizeBasisSum(viewer, 0);
			break;
		case '8':
			vectorFields.visualizeBasisSum(viewer, 1);
			break;
		case '9':
			vectorFields.visualizeApproximatedFields(viewer);
			vectorFields.visualizeGlobalConstraints(viewer);
			vectorFields.visualizeSingularitiesConstraints(viewer);
			break;
		case '0':
			vectorFields.visualizeApproxResult(viewer, 0);
			vectorFields.visualizeUserConstraints(viewer);
			//vectorFields.visualizeSamples(viewer);
			vectorFields.visualizeSingularitiesConstraints(viewer);
			break;
		case 'x':
		case 'X':
			selectFace = !selectFace;
			cout << "Select face: " << selectFace << endl; 
			//if (selectFace) cout << "Face is selected" << endl; 
			break;
		case 'c':
		case 'C':
			C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
			break; 
		default:
			break;
			return false;
		}
		//viewer.data.set_vertices(V);
		viewer.data().compute_normals();
		//viewer.core.align_camera_center(V, F);
		return true;
	};

	
	viewer.callback_mouse_move =
		[&V, &F, &C, &selectFace](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	{
		int fid;
		Eigen::Vector3f bc;
		// Cast a ray in the view direction starting from the mouse position
		double x = viewer.current_mouse_x;
		double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		if (selectFace) {
			if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view /** viewer.core.model*/,
				viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
			{

				// paint hit red
				C.row(fid) << 1, 0, 0;
				viewer.data().set_colors(C);
				return true;
			}
		}
		return false;
	};
	viewer.callback_init = [&](igl::opengl::glfw::Viewer &viewer) {return false; };
	//viewer.data().set_mesh(V, F);
	viewer.callback_key_down = key_down;

	Eigen::Vector4f bgCol(1.0, 1.0, 1.0, 1.0);
	viewer.core.background_color = bgCol;
	viewer.data().point_size = 10.0f;
	viewer.data().line_width = 1.0f; 

	return viewer.launch();
}