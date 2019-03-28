#include "VectorFields.h"
#include <igl/unproject_onto_mesh.h>

#include "TestSolver.h"


int eigToShow = 0, basisId = 0, selectedVertex;
int numSample = 500;
int eigToShow2 = 0;
int eigsToCompute = 20; 

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
	string meshFile = "../LocalFields/Models/Sphere/round_sphere_10242.obj";
	//string meshFile = "../LocalFields/Models/Thorus/Thorus_2304.obj";
	//string meshFile = "../LocalFields/Models/Thorus/torus.obj";

	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_1083.obj";
	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_10812.obj";
	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/894_dragon_tris.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/dragon_2000.obj";
	//string meshFile = "../LocalFields/Models/AIM_fertility_watertight/fertility.obj";
	//string meshFile = "../LocalFields/Models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off";
	//string meshFile = "../LocalFields/Models/Bunny/Bunny.obj";

	/* MODEL FOR TESTING, LARGE ONES */
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_4k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Rocker-arm/38_rocker-arm.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus5_long_36k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus5_33k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus2_60k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Brezel/Brezel_1920.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Neptune_clean__watertight_4M triangles/803_neptune_4Mtriangles_manifold.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_100.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Raptor/178_raptor.off";

	/* ========================= PRE-PROCESS ==============================*/
	cout << "========================= PRE-PROCESS ==============================\n"; 
	vectorFields.readMesh(meshFile);
	//Eigen::SparseMatrix<double> ChrisSparseMat;
	//ReadChristopherStiffnessMatrix("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Local Fields/Models/hodgeLaplace.txt", ChrisSparseMat);
	//WriteSparseMatrixToMatlab(ChrisSparseMat, "hello");
	//vectorFields.readArrowMesh("../LocalFields/Models/arrow.obj");
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
	vectorFields.selectFaceToDraw(15000);
	
	vectorFields.getVF(V, F);
	viewer.data().set_mesh(V, F);
	viewer.append_mesh();
	viewer.data().set_mesh(V, F);
	viewer.data().show_lines = false; 
	viewer.selected_data_index = 0; 

	/* MATRIX CONSTRUCTIONS */
	vectorFields.constructMassMatrices();
	vectorFields.constructRotationMatrix();
	vectorFields.constructMappingMatrix();
	
	/* =========== Test on PROBLEM SOLVING-related functionalities ================*/
	vectorFields.constructGradient3D();
	vectorFields.constructStiffnessMatrices();
	vectorFields.constructMatrixB();
	//vectorFields.constructConstraints();
	//vectorFields.checkB2DStructure();

	/* ====================== GLOBAL PROBLEM ====================*/
	//cout << "\n========================= GLOBAL PROBLEM =============================\n";
	//vectorFields.setupGlobalProblem();
	
	/* ====================== LOCAL ELEMENTS ====================*/
	cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
	vectorFields.constructSamples(numSample);
	vectorFields.constructBasis();
	//vectorFields.constructBasisEigenVects();
	//vectorFields.setupReducedBiLaplacian();
	//vectorFields.setAndSolveUserSystem();
	//vectorFields.measureApproxAccuracyL2Norm();
	//vectorFields.measureDirichletEnergy();

	//vectorFields.writeBasisToFile();
	//vectorFields.writeField3DToFile();
	//vectorFields.measureL2NormEigVectors();

	//vectorFields.computeEigenFields(eigsToCompute);
	//vectorFields.retrieveEigenFields();
	//vectorFields.computeApproxEigenFields(eigsToCompute);
	//vectorFields.retrieveApproxEigenFields();

	//vectorFields.testEnergyOfLocalPatch(viewer);

	//vectorFields.printDataForVTK();
	//vectorFields.writeEigenFieldsForVTK();

	/* ====================== TESTING BASIS ====================*/
	//vectorFields.constructArbitraryField();
	//vectorFields.constructArbitraryField2D();
	//double error; 
	//vectorFields.testBasis_NoRegularizer(error);
	//vectorFields.testBasis_WithRegularizer();
	//vectorFields.projectionTest();
	//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.1, 0.1, 0.8), 1.0);
	//vectorFields.visualizeApproximatedFields(viewer);	
	//vectorFields.visualize2Dfields(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 2.0, true);
	//vectorFields.visualize2Dfields(viewer, vectorFields.projRef, Eigen::RowVector3d(0.0, 0.9, 0.1), 2.0, true);
	//vectorFields.visualize2Dfields(viewer, vectorFields.projApprox, Eigen::RowVector3d(0.8, 0.1, 0.1), 2.0, true);
	//vectorFields.visualizeGlobalConstraints(viewer);
	//vectorFields.measureDirichletEnergy();

	/* ====================== PARALLEL TRANSPORT ====================*/
	//vectorFields.computeDijkstraForParallelTransport(200, 5000);
	//vectorFields.constructParallelTransport();
	//vectorFields.visualizeParallelTransportPath(viewer);
	//vectorFields.visualizeParallelTransport(viewer);

	/* ====================== APP: SMOOTHING VECTOR FIELDS  ====================*/
	//Eigen::VectorXd v_in = vectorFields.arbField2D;
	Eigen::VectorXd v_in = vectorFields.getRefFields();
	Eigen::VectorXd v_out;
	const double mu = 0.04; 

	/* ====================== APP: SMOOTHING TENSOR FIELDS (CURVATURE) ====================*/
	//vectorFields.ConstructCurvatureTensor(viewer);
	//vectorFields.ComputeCurvatureFields();

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
	//vectorFields.testEdgesAddition(viewer);
	//vectorFields.visualizePatchDijkstra(viewer);
	
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
	//vectorFields.visualizeCurvatureTensor(viewer);


	/* FOR GENERATING IMAGES on PAPER */
	vectorFields.visualizeSubdomain(viewer);
	bool evenSpaceField = true; 
	//vectorFields.visualizeSamples(viewer);

	/* Variables for faces of face selection */
	vector<int> ChosenFaces;
	Eigen::RowVector3d constraintDir;

	const auto &key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		int selectedFace;
		srand(time(NULL));
		selectedVertex = rand() % V.rows();

		switch (key)
		{
		case '1':
			//vectorFields.visualizeSubdomain(viewer);
			//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.1, 0.1, 0.8), 1.0);			
			vectorFields.visualizeApproximatedFields(viewer);
			vectorFields.visualizeGlobalConstraints(viewer);
			break;
		case '2':
			vectorFields.visualizeApproxResult(viewer);
			vectorFields.visualizeGlobalConstraints(viewer);
			//evenSpaceField = !evenSpaceField; 
			//vectorFields.visualize1FieldOnCenter(viewer, evenSpaceField);
			break;
		case '3':
			//vectorFields.visualizeGlobalConstraints(viewer);
			//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 1.0);
			//vectorFields.visualizeGlobalConstraints(viewer);
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
			eigToShow = max(eigToShow - 1, 0);
			vectorFields.visualizeEigenfields(viewer, eigToShow);
			printf("[Full] Eigen vector: %d (eigval=%.3f)\n", eigToShow, vectorFields.eigValuesFull(eigToShow));
			break;
		case '8':
			eigToShow = min(eigToShow + 1, 2 * 100 - 1);
			vectorFields.visualizeEigenfields(viewer, eigToShow);
			printf("[Full] Eigen vector: %d (eigval=%.3f)\n", eigToShow, vectorFields.eigValuesFull(eigToShow));
			break;
		case '9':
			//vectorFields.visualizeApproximatedFields(viewer);
			//vectorFields.visualizeGlobalConstraints(viewer);
			//vectorFields.visualizeSingularitiesConstraints(viewer);
			vectorFields.resetInteractiveConstraints();
			break;
		case '0':			
			viewer.selected_data_index = 0;
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			viewer.selected_data_index = 1;
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			viewer.selected_data_index = 0;			
			break;
		case 'y':
		case 'Y':
			eigToShow2 = max(eigToShow2 - 1, 0);
			vectorFields.visualizeApproxEigenfields(viewer, eigToShow2);
			printf("[Approx] Eigen vector: %d (eigval=%.3f)\n", eigToShow2, vectorFields.eigValuesReduced(eigToShow2));
			break;
		case 'u':
		case 'U':
			eigToShow2 = min(eigToShow2 + 1, eigsToCompute - 1);
			vectorFields.visualizeApproxEigenfields(viewer, eigToShow2);
			printf("[Approx] Eigen vector: %d (eigval=%.3f)\n", eigToShow2, vectorFields.eigValuesReduced(eigToShow2));
			break;

		//case 'b':
		//case 'B':
		//	vectorFields.visualizeLocalSubdomain(viewer);
		//	break;
		//case 'x':
		//case 'X':
			//C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
			//selectFace = !selectFace;
			//vectorFields.visualizeRandomFace(viewer, selectedFace);

			//cout << "Select face: " << selectFace << endl; 
			//if (selectFace) cout << "Face is selected" << endl; 
			//selectedFace = rand() % F.rows();
			//cout << "Face: " << selectedFace << endl;
		case 'x':
		case 'X':
			C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
			selectFace = !selectFace;			
			//vectorFields.visualizeRandomFace(viewer, selectedFace);
			break;
		case 'c':
		case 'C':
			printf("Computing smoothing on Reduced basis\n");
			vectorFields.computeSmoothingApprox(mu, v_in, v_out);			
			v_in = v_out; 
			break; 
		//case 'x':
		//case 'X':
		//	printf("Computing smoothing on full res\n");
		//	vectorFields.computeSmoothing(mu, v_in, v_out);
		//	vectorFields.visualize2Dfields(viewer, v_out, Eigen::RowVector3d(0.9, 0.1, 0.0), 2.0, true);
		//	//vectorFields.visualize2Dfields(viewer, v_in, Eigen::RowVector3d(0.0, 0.1, 0.9), 2.0, true);
		//	v_in = v_out;
		//	break; 
		case 'v':
		case 'V':
			vectorFields.visualize2DfieldsScaled(viewer, v_out, Eigen::RowVector3d(0.6, 0.2, 0.3), 1.0);
			break; 
		case 'b':
		case 'B':
			v_in = vectorFields.arbField2D; 
			break;
		case 'r':
		case 'R':
			cout << "\n========================= GLOBAL PROBLEM =============================\n";
			vectorFields.setupGlobalProblem();
			vectorFields.visualizeApproximatedFields(viewer);
			vectorFields.visualizeGlobalConstraints(viewer);
			break;
		case ' ':
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			//cout << "\n========================= GLOBAL PROBLEM =============================\n";
			//vectorFields.setupGlobalProblem();			
			//vectorFields.visualizeApproximatedFields(viewer);
			cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
			vectorFields.setAndSolveUserSystem();
			vectorFields.visualizeApproxResult(viewer);
			vectorFields.visualizeGlobalConstraints(viewer);
			break; 
		case 'Z':
		case 'z':
			vectorFields.write2DFieldsForVTK(vectorFields.XFullDim, "VectorFields_2", "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/VTK and ParaView/Test Data/VFieldsDes_SmallDragon_vBased.vtk");
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
		[&V, &F, &C, &selectFace, &ChosenFaces, &constraintDir, &vectorFields](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	{
		int fid;
		Eigen::Vector3f bc;
		/* Collect the selected faces */
		
		// Cast a ray in the view direction starting from the mouse position
		double x = viewer.current_mouse_x;
		double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		if (selectFace) {
			if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view /** viewer.core.model*/,
				viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
			{				
				// paint hit red
				C.row(fid) << 1, 0, 0;
				//viewer.data().set_colors(C);
				ChosenFaces.push_back(fid);
				if (ChosenFaces.size() == 1)
				{
					viewer.data().add_points(vectorFields.FC.row(fid), Eigen::RowVector3d(0.0, 0.1, 0.0));
				}

				return true;
			}
		}
		else {
			/* Getting the length + direction of the constraints*/
			const int constraintSize = ChosenFaces.size();
			if (constraintSize > 0)
			{
				//printf("Size of the contraints vector is %d\n", constraintSize);
				constraintDir = vectorFields.FC.row(ChosenFaces[constraintSize-1]) - vectorFields.FC.row(ChosenFaces[0]);
				viewer.data().add_edges(vectorFields.FC.row(ChosenFaces[0]), vectorFields.FC.row(ChosenFaces[0]) + constraintDir, Eigen::RowVector3d(1.0, 0.0, 0.1));
				vectorFields.pushNewUserConstraints(ChosenFaces[0], ChosenFaces[constraintSize - 1]);
				printf("Pair [%d]->[%d] is inserted\n", ChosenFaces[0], ChosenFaces[constraintSize - 1]);

				/* UPdating the mesh information */
				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				//cout << "\n========================= GLOBAL PROBLEM =============================\n";
				//vectorFields.setupGlobalProblem();			
				//vectorFields.visualizeApproximatedFields(viewer);
				cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
				vectorFields.setAndSolveUserSystem();
				vectorFields.visualizeApproxResult(viewer);
				vectorFields.visualizeGlobalConstraints(viewer);

				ChosenFaces.clear();
				//cout << "Hello there\n";
			}
			return false; 
		}
		return false;
	};
	viewer.callback_init = [&](igl::opengl::glfw::Viewer &viewer) {return false; };
	//viewer.data().set_mesh(V, F);
	viewer.callback_key_down = key_down;

	Eigen::Vector4f bgCol(1.0, 1.0, 1.0, 1.0);
	viewer.core.background_color = bgCol;
	viewer.data().point_size = 10.0f;
	viewer.data().line_width = 1.5f; 

	return viewer.launch();
}