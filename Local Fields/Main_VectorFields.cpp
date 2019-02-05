#include "VectorFields.h"
#include <igl/unproject_onto_mesh.h>

#include "TestSolver.h"

int eigToShow = 0, basisId=0, numSample=5000, selectedVertex;

int main(int argc, char *argv[])
{
	bool selectFace = false;
	Eigen::MatrixXd C;

	VectorFields vectorFields;
	igl::opengl::glfw::Viewer		viewer;
	Eigen::MatrixXd					V;
	Eigen::MatrixXi					F, E;

	// Hell there this is main function.

	/* READING DATA */
	
	//string meshFile = "../LocalFields/Models/Cube/Cube_1400.obj";
	//string meshFile = "../LocalFields/Models/Plane/square_plane.obj";

	//string meshFile = "../LocalFields/Models/Sphere/round_sphere_small.obj";
	//string meshFile = "../LocalFields/Models/Sphere/round_sphere_1500.obj";
	//string meshFile = "../LocalFields/Models/Sphere/round_sphere_10242.obj";
	//string meshFile = "../LocalFields/Models/Thorus/Thorus_2304.obj";
	//string meshFile = "../LocalFields/Models/Thorus/torus.obj";

	string meshFile = "../LocalFields/Models/Armadillo/Armadillo_1083.obj";
	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_10812.obj";
	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/894_dragon_tris.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/dragon_2000.obj";
	//string meshFile = "../LocalFields/Models/AIM_fertility_watertight/fertility.obj";
	//string meshFile = "../LocalFields/Models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off";
	//string meshFile = "../LocalFields/Models/Bunny/Bunny.obj";

	/* MODEL FOR TESTING, LARGE ONES */
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Neptune_clean__watertight_4M triangles/803_neptune_4Mtriangles_manifold.off";

	/* ========================= PRE-PROCESS ==============================*/
	cout << "========================= PRE-PROCESS ==============================\n"; 
	vectorFields.readMesh(meshFile);
	//vectorFields.computeEdges();
	vectorFields.computeAverageEdgeLength();
	vectorFields.computeFaceCenter();
	vectorFields.computeFaceNormal();
	vectorFields.constructVFNeighbors();
	//vectorFields.constructVFNeighborsFull();
	//vectorFields.constructVFAdjacency();
	//vectorFields.testAdjacency();
	vectorFields.constructFaceAdjacency3NMatrix();
	vectorFields.constructFaceAdjacency2RingMatrix();
	vectorFields.selectFaceToDraw(2500);
	
	vectorFields.getVF(V, F);
	viewer.data().set_mesh(V, F);

	/* MATRIX CONSTRUCTIONS */
	vectorFields.constructMassMatrices();
	vectorFields.constructRotationMatrix();
	vectorFields.constructMappingMatrix();
	
	/* =========== Test on PROBLEM SOLVING-related functionalities ================*/
	vectorFields.constructGradient3D();
	vectorFields.constructStiffnessMatrices();
	vectorFields.constructMatrixB();
	vectorFields.constructConstraints();
	//vectorFields.checkB2DStructure();

	/* ====================== GLOBAL PROBLEM ====================*/
	//cout << "\n========================= GLOBAL PROBLEM =============================\n";
	//vectorFields.setupGlobalProblem();
	
	/* ====================== LOCAL ELEMENTS ====================*/
	//cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
	//vectorFields.constructSamples(numSample);
	//vectorFields.constructBasis();
	//vectorFields.setAndSolveUserSystem();
	//vectorFields.measureApproxAccuracyL2Norm();

	/* ====================== TESTING BASIS ====================*/
	//vectorFields.constructArbitraryField();
	//vectorFields.constructArbitraryField2D();
	//vectorFields.testBasis();
	//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.1, 0.1, 0.8), 5000);
	//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 5000);

	//vectorFields.measureDirichletEnergy();


	/* ====================== APP: SMOOTHING VECTOR FIELDS  ====================*/
	Eigen::VectorXd v_in = vectorFields.arbField2D;
	Eigen::VectorXd v_out;
	const double mu = 4; 

	/* ====================== APP: SMOOTHING TENSOR FIELDS (CURVATURE) ====================*/
	vectorFields.ConstructCurvatureTensor();


	/* ==================== VISUALIZATION ======================== */
	/* GLOBAL  */
	//vectorFields.visualizeApproximatedFields(viewer);
	//vectorFields.visualizeGlobalConstraints(viewer);
	//vectorFields.visualizeSingularitiesConstraints(viewer);
	//vectorFields.visualizeSharedEdges(viewer);

	/* LOCAL  */
	//vectorFields.visualizeApproxResult(viewer);	
	//vectorFields.visualizeUserConstraints(viewer);
	//vectorFields.visualizeSamples(viewer);
	//vectorFields.visualizeSingularitiesConstraints(viewer);
	
	/* VISUALIZATION FOR TESTING PURPOSE */
	//vectorFields.visualizeNeighboringRings(viewer);
	//vectorFields.visualizeDijkstraFace(viewer);
	//vectorFields.visualizeArbField(viewer);
	//vectorFields.visualizeVertexFacesNeighbors(viewer, 0);
	
	//vectorFields.visualizeCurveConstraints(viewer);
	//vectorFields.visualizeSoftConstraints(viewer);

	/* MEASURE ACCURACY */
	//vectorFields.measureApproxAccuracyL2Norm();

	/* TEST SOLVERS */
	//testLAPACKEdsysv();
	//testCUDA_LULinearSolver();
	//testLDLTLinearSolver();
	//testMKL_Pardiso();

	/* VISUALIZATION OF APPLICATIONS */
	vectorFields.visualizeCurvatureTensor(viewer);

	const auto &key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		int selectedFace;
		srand(time(NULL));
		selectedVertex = rand() % V.rows();

		switch (key)
		{
		case '1':
			//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.1, 0.1, 0.8), 5000);			
			vectorFields.visualizeApproximatedFields(viewer);
			vectorFields.visualizeGlobalConstraints(viewer);
			break;
		case '2':
			vectorFields.visualizeApproxResult(viewer);
			vectorFields.visualizeUserConstraints(viewer);
			break;
		case '3':
			vectorFields.visualizeGlobalConstraints(viewer);
			//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 5000);
			//vectorFields.visualizeGlobalConstraints(viewer);
			//basisId = max(basisId - 1, 0);
			//vectorFields.visualizeBasis(viewer, basisId);
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
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			break;
		//case 'x':
		//case 'X':
			//C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
			//selectFace = !selectFace;
			//cout << "Select face: " << selectFace << endl; 
			//if (selectFace) cout << "Face is selected" << endl; 
			//selectedFace = rand() % F.rows();
			//cout << "Face: " << selectedFace << endl;
			//vectorFields.visualizeRandomFace(viewer, selectedFace);
			break;
		case 'c':
		case 'C':
			vectorFields.computeSmoothingApprox(mu, v_in, v_out);			
			v_in = v_out; 
			break; 
		case 'x':
		case 'X':
			vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.6, 0.2, 0.3), 100);
			break;
		case 'v':
		case 'V':
			vectorFields.visualize2DfieldsScaled(viewer, v_in, Eigen::RowVector3d(0.6, 0.2, 0.3), 100);
			break; 
		case 'b':
		case 'B':
			vectorFields.visualize2DfieldsScaled(viewer, v_out, Eigen::RowVector3d(0.6, 0.2, 0.3), 100);
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
				cout << "Face " << fid << " is selected." << endl; 
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