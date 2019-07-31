#include "VectorFields.h"
#include "NRoSyFields.h"
#include "TensorFields.h"

#include <igl/unproject_onto_mesh.h>

#include "TestSolver.h"

int eigToShow = 0, basisId = 0, selectedVertex;
int numSample = 50;
int eigToShow2 = 0;
int eigsToCompute = 500; 
int vfSaveId = 0;

enum class FieldsType {VECTOR, NROSY, TENSOR};
FieldsType fieldsType = FieldsType::NROSY;


int main(int argc, char *argv[])
{
	/* TEST MATLAB DATA C++ */
	//evalSurfaceGraph();

	bool selectFace = false;
	Eigen::MatrixXd C;

	
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

	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_1083.obj";
	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_10812.obj";	
	string meshFile = "../LocalFields/Models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/894_dragon_tris.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/dragon_2000.obj";
	//string meshFile = "../LocalFields/Models/AIM_fertility_watertight/fertility.obj";
	//string meshFile = "../LocalFields/Models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off";
	//string meshFile = "../LocalFields/Models/Bunny/Bunny.obj";

	/* MODEL FOR TESTING, LARGE ONES */
	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_4k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Kitten-watertight/366_kitten_5000.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Kitten-watertight/366_kitten_final.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Bimba_1M faces_clean_watertight/272_bimba_clean_1Mf.obj";	
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Rocker-arm/38_rocker-arm.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus5_long_36k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus5_33k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus2_60k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Brezel/Brezel_1920.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Armadillo/Armadillo_2525.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Neptune_clean__watertight_4M triangles/803_neptune_4Mtriangles_manifold.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Raptor/178_raptor.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Cube/Cube_round_50k_2.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Cube/Cube_sharp_50k_2.obj";
	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";

	/* ====================== VECTOR FIELDS  ====================*/
	

	/* ====================== APP: SMOOTHING VECTOR FIELDS  ====================*/
	//Eigen::VectorXd v_in = vectorFields.arbField2D;
	//Eigen::VectorXd v_in = vectorFields.getRefFields();
	Eigen::VectorXd v_in;
	Eigen::VectorXd v_out;
	const double mu = 0.04; 

	/* ====================== APP: SMOOTHING TENSOR FIELDS (CURVATURE) ====================*/
	cout << "\n========================= TENSOR FIELDS =============================\n";
	
	TensorFields tensorFields;
	VectorFields vectorFields;
	NRoSyFields nRoSyFields;

	switch (fieldsType)
	{
	case FieldsType::VECTOR:
		cout << "\n========================= VECTOR FIELDS =============================\n";
		vectorFields.TEST_VECTOR(viewer, meshFile);
		vectorFields.getVF(V, F);
		break;
	case FieldsType::NROSY:
		cout << "\n========================= N-ROSY FIELDS =============================\n";
		nRoSyFields.TEST_NROSY(viewer, meshFile);
		nRoSyFields.getVF(V, F);
		break;
	case FieldsType::TENSOR:
		cout << "\n========================= TENSOR FIELDS =============================\n";
		tensorFields.TEST_TENSOR(viewer, meshFile);
		tensorFields.getVF(V, F);
		break;
	default:
		break;
	}
	
	Eigen::MatrixXd inputTensor = tensorFields.Tensor;
	Eigen::MatrixXd outputTensor;
	Eigen::MatrixXd inputTensorRed = tensorFields.Tensor;
	Eigen::MatrixXd outputTensorRed;

	//viewer.data().add_points(V.row(0), Eigen::RowVector3d(1, 0, 0));
		

	/* TEST SOLVERS */
	//testLAPACKEdsysv();
	//testCUDA_LULinearSolver();
	//testLDLTLinearSolver();
	//testMKL_Pardiso();


	/* FOR GENERATING IMAGES on PAPER */
	//vectorFields.visualizeSubdomain(viewer);
	bool evenSpaceField = true; 
	//vectorFields.visualizeSamples(viewer);

	/* Variables for faces of face selection */
	vector<int> ChosenFaces;
	Eigen::RowVector3d constraintDir;

	/* ====================== N-RoSy Fields  ====================*/
	///NRoSyFields nRoSyFields;
	///nRoSyFields.readMesh(V, F);
	///nRoSyFields.constructNRoSyFields(4, vectorFields.eigFieldFull2D);
	///nRoSyFields.createNRoSyFromVectors(vectorFields.eigFieldFull2D.col(0));
	/////nRoSyFields.visualizeNRoSyFields(viewer);
	///nRoSyFields.visualizeRepVectorFields(viewer);

	string filename_vfields;
	Eigen::Vector3d lambda(1, 1e-4, 0.9);

	bool showSmoothed = false; 
	bool represent_tensor_in_voigt = false; 

	/* N-RoSy stuff */
	NRoSy nRoSy; 

	const auto &key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		int selectedFace;
		srand(time(NULL));
		selectedVertex = rand() % V.rows();

		/* For selecting faces */
		double xMouse;
		double yMouse;
		int fid;
		Eigen::Vector3f bc;



		switch (key)
		{
		case '-':
			//vectorFields.projectionTest();
			//vectorFields.visualize2Dfields(viewer, vectorFields.projRef, Eigen::RowVector3d(0.0, 0.9, 0.1), 2.0, true);
			//vectorFields.visualize2Dfields(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 3.0, false);	
			vectorFields.visualize2Dfields(viewer, vectorFields.wbEigen, Eigen::RowVector3d(0.8, 0.1, 0.8), 300.0, false);
			//vectorFields.visualizeGlobalConstraints(viewer); 	
			break;
		case '`':			
			//vectorFields.visualizeSubdomain(viewer);
			///vectorFields.visualize2Dfields(viewer, vectorFields.projRef, Eigen::RowVector3d(0.1, 0.1, 0.8), 3.0, false);			
			vectorFields.visualizeSoftConstraints(viewer);
			break;
		case '1':
			if (fieldsType == FieldsType::VECTOR)
			{				
				//vectorFields.visualizeSubdomain(viewer);
				//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.1, 0.1, 0.8), 1.0);			
				vectorFields.visualizeApproximatedFields(viewer);
				//vectorFields.visualizeGlobalConstraints(viewer);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				nRoSyFields.visualizeSoftConstraints(viewer);
				nRoSyFields.visualizeConstrainedFields(viewer);
				//nRoSyFields.visualizeConstraints(viewer);
				
				//nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.Xf, nRoSy);
				//nRoSyFields.visualizeNRoSyFields(viewer, nRoSy, Eigen::RowVector3d(0.1, 0.1, 0.8));

			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
			}			
			break;
		case '2':
			if (fieldsType == FieldsType::VECTOR)
			{
				///vectorFields.visualize2Dfields(viewer, vectorFields.pertFields, Eigen::RowVector3d(0.0, 0.9, 0.1), 3.0, false);
				//vectorFields.visualizeApproxResult(viewer);
				
				vectorFields.visualize2Dfields(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 3.0, false);
				//vectorFields.visualizeGlobalConstraints(viewer);
				//evenSpaceField = !evenSpaceField; 
				//vectorFields.visualize1FieldOnCenter(viewer, evenSpaceField);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				//nRoSyFields.visualizeSoftConstraints(viewer);
				nRoSyFields.visualizeConstrainedFields_Reduced(viewer);
				//nRoSyFields.visualizeConstraints(viewer);

				//nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.wb, nRoSy);
				//nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.XfBar, nRoSy);
				//nRoSyFields.visualizeNRoSyFields(viewer, nRoSy, Eigen::RowVector3d(0.8, 0.1, 0.1));
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				tensorFields.visualizeReducedTensorFields(viewer);
			}
			
			break;
		case '3':
			//vectorFields.visualizeGlobalConstraints(viewer);
			//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 1.0);
			//vectorFields.visualizeGlobalConstraints(viewer);
			basisId = max(basisId - 1, 0);		

			if (fieldsType == FieldsType::VECTOR)
			{
				vectorFields.visualizeBasis(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				nRoSyFields.visualizeBasis(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				
			}
			break;
		case '4':
			basisId = min(basisId + 1, 2*numSample-1);
			vectorFields.visualizeBasis(viewer, basisId);
			break;
		case '5':
			basisId = max(basisId - 1, 0);

			if (fieldsType == FieldsType::VECTOR)
			{
				vectorFields.visualizeBasisNormalized(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				nRoSyFields.visualizeBasis(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				tensorFields.visualizeBasis(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			
			break;
		case '6':
			basisId = min(basisId + 1, 2 * numSample - 1);
			if (fieldsType == FieldsType::VECTOR)
			{
				vectorFields.visualizeBasisNormalized(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				nRoSyFields.visualizeBasis(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				tensorFields.visualizeBasis(viewer, basisId);
				printf("Showing the %d-th basis \n", basisId);
			}
			break;
		case '7':
			eigToShow = max(eigToShow - 1, 0);
			if (fieldsType == FieldsType::VECTOR) 
			{				
				vectorFields.visualizeEigenfields(viewer, eigToShow);
				printf("[Full] Eigen vector: %d (eigval=%.3f)\n", eigToShow, vectorFields.eigValuesFull(eigToShow));
			}
			else if (fieldsType == FieldsType::NROSY) 
			{
				//nRoSyFields.visualizeRepVectorFields(viewer, nRoSyFields.eigFieldsNRoSyRef.col(eigToShow));
				nRoSyFields.visualizeEigenFields(viewer, eigToShow);
				printf("[Full] Eigen vector: %d (eigval=%.3f)\n", eigToShow, nRoSyFields.eigValuesNRoSyRef(eigToShow));
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				tensorFields.visualizeEigenTensorFields(viewer, eigToShow);
				printf("[Ref] Eigen tensor fields: %d (eigval=%.3f)\n", eigToShow, tensorFields.eigValuesTensorRef(eigToShow));
			}
			break;
		case '8':
			eigToShow = min(eigToShow + 1, eigsToCompute - 1);
			if (fieldsType == FieldsType::VECTOR)
			{
				vectorFields.visualizeEigenfields(viewer, eigToShow);
				printf("[Full] Eigen vector: %d (eigval=%.3f)\n", eigToShow, vectorFields.eigValuesFull(eigToShow));
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				//nRoSyFields.visualizeRepVectorFields(viewer, nRoSyFields.eigFieldsNRoSyRef.col(eigToShow));
				nRoSyFields.visualizeEigenFields(viewer, eigToShow);
				printf("[Full] Eigen vector: %d (eigval=%.3f)\n", eigToShow, nRoSyFields.eigValuesNRoSyRef(eigToShow));
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				tensorFields.visualizeEigenTensorFields(viewer, eigToShow);
				printf("[Ref] Eigen tensor fields: %d (eigval=%.3f)\n", eigToShow, tensorFields.eigValuesTensorRef(eigToShow));
			}
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
			vectorFields.visualizeApproxEigenfields(viewer, eigToShow2, eigToShow);
			printf("[Approx] Eigen vector: %d (eigval=%.3f)\n", eigToShow2, vectorFields.eigValuesReduced(eigToShow2));
			break;
		case 'u':
		case 'U':
			eigToShow2 = min(eigToShow2 + 1, eigsToCompute - 1);
			vectorFields.visualizeApproxEigenfields(viewer, eigToShow2, eigToShow);
			printf("[Approx] Eigen vector: %d (eigval=%.3f)\n", eigToShow2, vectorFields.eigValuesReduced(eigToShow2));
			break;

		//case 'b':
		//case 'B':
		//	vectorFields.visualizeLocalSubdomain(viewer);
		//	break;
		case 'n':
		case 'N':
			//selectFace = true; 
			xMouse = viewer.current_mouse_x;
			yMouse = viewer.core.viewport(3) - viewer.current_mouse_y;
			
			//if (selectFace) {
				if (igl::unproject_onto_mesh(Eigen::Vector2f(xMouse, yMouse), viewer.core.view /** viewer.core.model*/,
					viewer.core.proj, viewer.core.viewport, V, F, selectedFace, bc))
				{
					// paint hit red
					printf("Selected face is %d \n", selectedFace);
					C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
					C.row(selectedFace) << 1, 0, 0;
					viewer.data().set_colors(C);					
					//viewer.data().add_points(vectorFields.FC.row(selectedFace), Eigen::RowVector3d(0.0, 0.1, 0.0));

					return true;
				}
			//}
			break; 
		case 'M':
		case 'm':
			showSmoothed = !showSmoothed;
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			if (!showSmoothed) {				
				tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
			}
			else
			{
				tensorFields.smoothingRef(viewer, inputTensor, outputTensor);
				tensorFields.visualizeSmoothedTensorFields(viewer);
				inputTensor = outputTensor;
			}
			break;
		// Reduced smoothing
		case ',':
			showSmoothed = !showSmoothed;
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			if (!showSmoothed) {
				tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
			}
			else
			{
				tensorFields.smoothingRed(viewer, inputTensorRed, outputTensorRed);
				tensorFields.visualizeSmoothedAppTensorFields(viewer);
				inputTensorRed = outputTensorRed;
			}
			break;

		/* Case x to activate for user inputted constraints */	
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
			represent_tensor_in_voigt = !represent_tensor_in_voigt;
			if (fieldsType == FieldsType::TENSOR)
			{
				if (represent_tensor_in_voigt)
				{
					cout << "Representing in double conversion \n";
					tensorFields.tensorConvertNConvert(viewer);
				}
				else
				{
					cout << "Representing the original data \n";
					tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
				}
			}
			break;
		case 'r':
		case 'R':
			cout << "\n========================= GLOBAL PROBLEM =============================\n";
			vectorFields.setupGlobalProblem(lambda);
			vectorFields.visualizeApproximatedFields(viewer);
			//vectorFields.visualizeGlobalConstraints(viewer);
			break;
		case 's':
		case 'S':
			filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_" + std::to_string(vfSaveId) + ".txt";
			vfSaveId++;
			WriteEigenVectorToTxtFile(vectorFields.arbField2D, filename_vfields);
			break;
		case ' ':
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			//cout << "\n========================= GLOBAL PROBLEM =============================\n";
			//vectorFields.setupGlobalProblem();			
			//vectorFields.visualizeApproximatedFields(viewer);
			cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
			vectorFields.setAndSolveUserSystem(lambda);
			vectorFields.visualizeApproxResult(viewer);
			//vectorFields.visualizeGlobalConstraints(viewer);
			
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
		[&V, &F, &C, &lambda, &selectFace, &ChosenFaces, &constraintDir, &vectorFields](igl::opengl::glfw::Viewer& viewer, int, int)->bool
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
				//vectorFields.setupGlobalProblem(lambda);			
				//vectorFields.visualizeApproximatedFields(viewer);
				//vectorFields.visualizeGlobalConstraints(viewer);
				//vectorFields.visualizeAreaOfLaplaceConstraint(viewer);
				///cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
				///vectorFields.setAndSolveUserSystem(lambda);
				///vectorFields.visualizeApproxResult(viewer);
				///vectorFields.visualizeUserConstraints(viewer);
	
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
	viewer.data().line_width = 1.0f; 

	return viewer.launch();

	/* Trick for remote desktop */
	//getchar();
	//return 1;
}