#include "VectorFields.h"
#include "NRoSyFields.h"
#include "TensorFields.h"

#include <igl/unproject_onto_mesh.h>

#include "TestSolver.h"

int eigToShow = 0, basisId = 0, selectedVertex;
int numSample = 10000;
int eigToShow2 = 0;
int eigsToCompute = 500; 
int saveId = 0;
bool singularConstraint = false; 

enum class FieldsType {VECTOR, NROSY, TENSOR};
FieldsType fieldsType = FieldsType::VECTOR;



int main(int argc, char *argv[])
{
	/* TEST MATLAB DATA C++ */
	//evalSurfaceGraph();

	bool selectFace = false;
	Eigen::MatrixXd C;

	HANDLE msHandle; 
	LPTSTR slotName = TEXT("\\\\.\\mailslot\\sandy2");

	msHandle = CreateMailslot(slotName, 0, MAILSLOT_WAIT_FOREVER, (LPSECURITY_ATTRIBUTES)NULL);
	if (msHandle == INVALID_HANDLE_VALUE)  cout << "Cannto create mailslot\n"; else cout << "Creation of mailslot " << slotName << " is successful\n";

	
	igl::opengl::glfw::Viewer		viewer;
	Eigen::MatrixXd					V;
	Eigen::MatrixXi					F, E;

	// Hell there this is main function.

	/* READING DATA */
	const string model = "Fertility_";
	
	//string meshFile = "../LocalFields/Models/Cube/Cube_1400.obj";
	//string meshFile = "../LocalFields/Models/Plane/square_plane.obj";

	//string meshFile = "../LocalFields/Models/Sphere/round_sphere_small.obj";
	//string meshFile = "../LocalFields/Models/Sphere/round_sphere_1500.obj";
	//string meshFile = "../LocalFields/Models/Sphere/round_sphere_10242.obj";
	//string meshFile = "../LocalFields/Models/Thorus/Thorus_2304.obj";	
	//string meshFile = "../LocalFields/Models/Thorus/torus.obj";

	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_1083.obj";
	string meshFile = "../LocalFields/Models/Armadillo/Armadillo_10812.obj";	
	//string meshFile = "../LocalFields/Models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/894_dragon_tris.obj";
	//string meshFile = "../LocalFields/Models/AIM894_Chinese Dragon/dragon_2000.obj";
	//string meshFile = "../LocalFields/Models/AIM_fertility_watertight/fertility.obj";
	//string meshFile = "../LocalFields/Models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off";
	//string meshFile = "../LocalFields/Models/Bunny/Bunny.obj";

	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/heart/Heart3.obj";
	/* MODEL FOR TESTING, LARGE ONES */
	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_4k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";
	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Kitten-watertight/366_kitten_5000.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Kitten-watertight/366_kitten_final.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Vase-lion/177_vase-lion.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Bimba_1M faces_clean_watertight/bimba.obj";	
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/TOSCA_hires-mat/wolf_500k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Raptor/raptor.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/octopus_large/octopus_large.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/TOSCA_hires-mat/centaur1_425k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/TOSCA_hires-mat/cat4_750k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Rocker-arm/38_rocker-arm.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Pulley_full/pulley_40k.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Rocker-arm/38_rocker-arm_800k.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/blade_smooth/blade_smooth.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus5_long_36k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus5_33k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/HighGenus/Genus2_60k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Brezel/Brezel_1920.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Armadillo/Armadillo_43243.obj";
	//string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Armadillo/Armadillo_2525.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Neptune_clean__watertight_4M triangles/803_neptune_4Mtriangles_manifold.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Isidore_horse/424_Isidore_horse.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Raptor/178_raptor.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Cube/Cube_round_50k_2.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Cube/Cube_sharp_50k_2.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Thorus_73k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Thorus/Torus_3k_jv.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Plane/squarePlane_16k.obj";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Sphere/sphere10k.off";
	//string meshFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/Dragon/Dragon_150k.off";


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
		//vectorFields.computeApproxEigenFields_Mult();
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
		//viewer.data().set_colors(Eigen::RowVector3d(186.0/255.0, 225.0/255.0, 255.0/255.0));
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

	double tensor_lambda = 0.1;

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
			//vectorFields.visualize2Dfields(viewer, vectorFields.wbEigen, Eigen::RowVector3d(0.8, 0.1, 0.8), 3.0, false);
			//vectorFields.visualize2Dfields(viewer, vectorFields.wbEigen, Eigen::RowVector3d(180.0/255.0, 173.0/255.0, 234.0 / 255.0), 3.0, false);
			vectorFields.visualize2Dfields(viewer, vectorFields.wbEigen, Eigen::RowVector3d(0.0 / 255.0, 200.0 / 255.0, 0.0 / 255.0), 3.0, false);
			vectorFields.selectAdaptiveRegions(viewer);
			//vectorFields.visualizeGlobalConstraints(viewer); 	
			break;
		case '`':			
			//vectorFields.visualizeSubdomain(viewer);
			///vectorFields.visualize2Dfields(viewer, vectorFields.projRef, Eigen::RowVector3d(0.1, 0.1, 0.8), 3.0, false);			
			//vectorFields.visualizeSoftConstraints(viewer);
			vectorFields.visualizeSamples(viewer);
			break;
		case '1':
			if (fieldsType == FieldsType::VECTOR)
			{				
				//vectorFields.visualizeSubdomain(viewer);
				//vectorFields.visualize2DfieldsScaled(viewer, vectorFields.arbField2D, Eigen::RowVector3d(0.1, 0.1, 0.8), 1.0);			
				vectorFields.visualizeApproximatedFields(viewer);
				vectorFields.visualizeGlobalConstraints(viewer);
				//vectorFields.visualizeSingularitiesConstraints(viewer);				
				
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				//nRoSyFields.visualizeSoftConstraints(viewer);
				nRoSyFields.visualizeConstraints(viewer);
				nRoSyFields.visualizeConstrainedFields(viewer);				
				
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
				vectorFields.visualizeApproxResult(viewer);
				
				//vectorFields.visualize2Dfields(viewer, vectorFields.wb, Eigen::RowVector3d(0.8, 0.1, 0.1), 3.0, false);
				vectorFields.visualizeGlobalConstraints(viewer);
				//evenSpaceField = !evenSpaceField; 
				//vectorFields.visualize1FieldOnCenter(viewer, evenSpaceField);
				//vectorFields.selectAdaptiveRegions(viewer);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				//nRoSyFields.visualizeSoftConstraints(viewer);
				nRoSyFields.visualizeConstraints(viewer);
				//nRoSyFields.visualizeConstrainedFields_Reduced(viewer);

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
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
			///tensorFields.visualizeSubdomain(viewer);
			if(basisId%2==0)
				basisId = max(basisId - 1, 0);
			else 
				basisId = max(basisId - 1, 0);

			if (fieldsType == FieldsType::VECTOR)
			{
				//vectorFields.visualizeSubdomain(viewer);
				printf("Showing the %d-th basis \n", basisId);
				vectorFields.visualizeBasisNormalized(viewer, basisId);
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
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
			///tensorFields.visualizeSubdomain(viewer);
			if(basisId%2==1)
				basisId = min(basisId + 1, 2 * numSample - 1);
			else 
				basisId = min(basisId + 1, 2 * numSample - 1);
			if (fieldsType == FieldsType::VECTOR)
			{
				printf("Showing the %d-th basis \n", basisId);
				//vectorFields.visualizeSubdomain(viewer);
				vectorFields.visualizeBasisNormalized(viewer, basisId);
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
			nRoSyFields.resetInteractiveConstraints();
			break;
		case '0':			
			viewer.selected_data_index = 0;
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
			//viewer.selected_data_index = 1;
			//viewer.data().clear();
			//viewer.data().set_mesh(V, F);
			//viewer.selected_data_index = 0;			
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
			viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
			//viewer.data().set_colors(Eigen::RowVector3d(186.0 / 255.0, 225.0 / 255.0, 255.0 / 255.0));
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
			viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
			//viewer.data().set_colors(Eigen::RowVector3d(186.0 / 255.0, 225.0 / 255.0, 255.0 / 255.0));
			if (!showSmoothed) {
				tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
			}
			else
			{
				tensorFields.smoothingRed(viewer, inputTensorRed, tensor_lambda, outputTensorRed);
				tensorFields.visualizeSmoothedAppTensorFields(viewer);
				inputTensorRed = outputTensorRed;
			}
			break;

		case 'g':
		case 'G':
			if (fieldsType == FieldsType::VECTOR)
			{
				singularConstraint = !singularConstraint;
			}
			break;
		/* Case x to activate for user inputted constraints */	
		case 'x':
		case 'X':
			if (fieldsType == FieldsType::VECTOR)
			{
				C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
				selectFace = !selectFace;			
				//vectorFields.visualizeRandomFace(viewer, selectedFace);
			} 
			else if (fieldsType == FieldsType::TENSOR)
			{
				showSmoothed = !showSmoothed;
				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
				//viewer.data().set_colors(Eigen::RowVector3d(186.0 / 255.0, 225.0 / 255.0, 255.0 / 255.0));
				if (!showSmoothed) {
					tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
				}
				else
				{
					tensor_lambda /= 5.0;
					printf("lambda=%.3f \n", tensor_lambda);
					tensorFields.smoothingRed(viewer, inputTensorRed, tensor_lambda, outputTensorRed);
					tensorFields.visualizeSmoothedAppTensorFields(viewer);
					//inputTensorRed = outputTensorRed;
				}
			}
			break;
		case 'c':
		case 'C':			
			if (fieldsType == FieldsType::VECTOR)
			{
				printf("Computing smoothing on Reduced basis\n");
				vectorFields.computeSmoothingApprox(mu, v_in, v_out);
				v_in = v_out;
			}
			else if(fieldsType==FieldsType::TENSOR)
			{
				showSmoothed = !showSmoothed;
				viewer.data().clear();
				viewer.data().set_mesh(V, F);
				viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
				//viewer.data().set_colors(Eigen::RowVector3d(186.0 / 255.0, 225.0 / 255.0, 255.0 / 255.0));
				if (!showSmoothed) {
					tensorFields.visualizeTensorFields(viewer, tensorFields.tensorFields);
				}
				else
				{
					tensor_lambda *= 5.0;
					printf("lambda=%.3f \n", tensor_lambda);
					tensorFields.smoothingRed(viewer, inputTensorRed, tensor_lambda, outputTensorRed);
					tensorFields.visualizeSmoothedAppTensorFields(viewer);
					//inputTensorRed = outputTensorRed;
				}
			}
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
			if (fieldsType == FieldsType::NROSY) {				
				nRoSyFields.readFieldsFromMailSlot(msHandle);
			}
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
			if (fieldsType == FieldsType::VECTOR)
			{
				cout << "\n========================= GLOBAL PROBLEM =============================\n";
				vectorFields.constructInteractiveConstraints();
				//vectorFields.constructHardConstraintsWithSingularitiesWithGauss(viewer);
				vectorFields.setupGlobalProblem(lambda);
				vectorFields.visualizeApproximatedFields(viewer);
				vectorFields.visualizeGlobalConstraints(viewer);

				//nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.XfBar, nRoSy);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				nRoSyFields.nRoSyFieldsDesignRef();
				viewer.data().clear(); viewer.data().set_mesh(V, F);
				nRoSyFields.visualizeConstraints(viewer);
				nRoSyFields.visualizeConstrainedFields(viewer);

				nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.alignFields, nRoSy);
				nRoSyFields.visualizeNRoSyFields(viewer, nRoSy, Eigen::RowVector3d(0.0, 0.7, 0.0));
			}
			else if (fieldsType == FieldsType::TENSOR)
			{

			}
			break;
		case 's':
		case 'S':
			if (fieldsType == FieldsType::VECTOR)
			{
				filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/"+model+"constraintFields_user_" + std::to_string(saveId);
				//filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/" + model + "constraintFields_Local_user_" + std::to_string(saveId);
				saveId++;
				vectorFields.writeVectorFieldsToFile(vectorFields.XFullDim, filename_vfields + ".txt");
				vectorFields.writeConstraintsToFile(filename_vfields + "_constraints.txt");

				//vectorFields.writeVectorFieldsToFile_Local(vectorFields.XFullDim, filename_vfields + "_red.txt");
				//vectorFields.writeVectorFieldsToFile_Local(vectorFields.Xf, filename_vfields + "_ref.txt");
				//vectorFields.writeConstraintsToFile_Local(filename_vfields + "_constraints.txt");
				//WriteEigenVectorToTxtFile(vectorFields.arbField2D, filename_vfields);

				//vectorFields.loadVectorFieldsFromFile(filename_vfields + ".txt", vectorFields.XFullDim);
				//vectorFields.visualize2Dfields(viewer, vectorFields.XFullDim, Eigen::RowVector3d(0.9, 0.1, 0.1), 4.0, false);
				//vectorFields.loadConstraintsFromFile(filename_vfields + "_constraints.txt");
				//vectorFields.visualizeGlobalConstraints(viewer);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				/* Saving the file */
				//filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/" + model +"2fields_" + to_string(nRoSyFields.Basis.cols()) + "_" + to_string((int)nRoSyFields.numSupport) + "_Local" + to_string(saveId++);
				filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/" + model + "2fields_" + to_string(nRoSyFields.Basis.cols()) + "_" + to_string((int)nRoSyFields.numSupport) + "_approx_" + to_string(saveId++);
				nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.XfBar, nRoSy);
				//nRoSyFields.writeNRoSyFieldsToFile_Local(nRoSy, filename_vfields+"_red.txt");
				nRoSyFields.writeNRoSyFieldsToFile(nRoSy, filename_vfields + ".txt");				
				nRoSyFields.writeConstraintsToFile(filename_vfields + "_constraints.txt");

				//nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.Xf, nRoSy);
				//nRoSyFields.writeNRoSyFieldsToFile_Local(nRoSy, filename_vfields + "_ref.txt");

				/* Load and display the fields */
				//cout << "========================================== \n";
				//cout << "Loading the " << nCounter-1 << " fields \n";
				//nRoSyFields.loadConstraintsFromFile(filename_vfields + "_constraints.txt");
				//nRoSyFields.loadNRoSyFieldsFromFile(filename_vfields + ".txt", nRoSy);
				//
				////nRoSyFields.visualizeConstraints(viewer);
				//nRoSyFields.visualizeNRoSyFields(viewer, nRoSy, Eigen::RowVector3d(0.9, 0.1, 0.1));
			}
			else if (fieldsType == FieldsType::TENSOR)
			{

			}
			break;
		case 'k':
		case 'K':
			if (fieldsType == FieldsType::VECTOR)
			{
				//filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_constraintFields_user_" + std::to_string(saveId);
				filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_constraintFields_Local_user_" + std::to_string(saveId);
				saveId++;
				//vectorFields.writeVectorFieldsToFile(vectorFields.XFullDim, filename_vfields + ".txt");
				//vectorFields.writeConstraintsToFile(filename_vfields + "_constraints.txt");
				//WriteEigenVectorToTxtFile(vectorFields.arbField2D, filename_vfields);

				//vectorFields.loadVectorFieldsFromFile(filename_vfields + ".txt", vectorFields.XFullDim);
				//vectorFields.visualize2Dfields(viewer, vectorFields.XFullDim, Eigen::RowVector3d(0.9, 0.1, 0.1), 4.0, false);
				vectorFields.loadConstraintsFromFile(filename_vfields + "_constraints.txt");
				vectorFields.visualizeGlobalConstraints(viewer);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				//filename_vfields = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/" + model + "_2fields_" + to_string(nRoSyFields.Basis.cols()) + "_" + to_string((int)nRoSyFields.numSupport) + "_approx" + to_string(saveId++);
				filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_constraintFields_Local_user_10_constraints.txt";
				nRoSyFields.loadConstraintsFromFile(filename_vfields);
				nRoSyFields.visualizeConstraints(viewer);
			}
			else if (fieldsType == FieldsType::TENSOR)
			{

			}
			break;
		case ' ':
			//viewer.data().clear();
			//viewer.data().set_mesh(V, F);

			if (fieldsType == FieldsType::VECTOR)
			{
				//cout << "\n========================= GLOBAL PROBLEM =============================\n";
				//vectorFields.setupGlobalProblem();			
				//vectorFields.visualizeApproximatedFields(viewer);
				cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
				//vectorFields.constructInteractiveConstraints();
				//vectorFields.setAndSolveInteractiveSystem(lambda);
				vectorFields.setAndSolveUserSystem(lambda);
				vectorFields.visualizeApproxResult(viewer);
				//vectorFields.visualizeGlobalConstraints(viewer);
			}
			else if (fieldsType == FieldsType::NROSY)
			{
				cout << "\n========================= N-ROSY FIELDS =============================\n";
				/* Full res*/
				//if (F.rows() < 50000)
				//{
				//	nRoSyFields.nRoSyFieldsDesignRef();
				//	nRoSyFields.visualizeConstraints(viewer);
				//	nRoSyFields.visualizeConstrainedFields(viewer);
				//}


				/* Reduced */
				//nRoSyFields.constructInteractiveConstraints();
				//nRoSyFields.nRoSyFieldsDesign_Reduced_HardConstraints();
				//nRoSyFields.nRoSyFieldsDesign_Reduced_Splines();
				nRoSyFields.setAndSolveInteractiveSystem();
				
				//viewer.data().clear(); viewer.data().set_mesh(V, F);
				nRoSyFields.visualizeConstraints(viewer);
				nRoSyFields.visualizeConstrainedFields_Reduced(viewer);

				//nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.alignFields, nRoSy);
				//nRoSyFields.visualizeNRoSyFields(viewer, nRoSy, Eigen::RowVector3d(0.0, 0.9, 0.1));				
			}
			else if (fieldsType == FieldsType::TENSOR)
			{
				
			}
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
		[&V, &F, &C, &lambda, &selectFace, &ChosenFaces, &constraintDir, &vectorFields, &nRoSyFields](igl::opengl::glfw::Viewer& viewer, int, int)->bool
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
				if (singularConstraint)
				{
					ChosenFaces.push_back(fid);
					selectFace = !selectFace;
					for (int i : ChosenFaces)
						viewer.data().add_points(vectorFields.V.row(F(fid,0)), Eigen::RowVector3d(0.0, 0.8, 0.0));
						//viewer.data().add_points(vectorFields.FC.row(fid), Eigen::RowVector3d(0.0, 0.8, 0.0));
				}
				else {
					// paint hit red
					C.row(fid) << 1, 0, 0;
					//viewer.data().set_colors(C);
					ChosenFaces.push_back(fid);
					if (ChosenFaces.size() == 1)
					{
						viewer.data().add_points(vectorFields.FC.row(fid), Eigen::RowVector3d(0.0, 0.1, 0.0));
					}
				}

				/* Select vertices closest to the center*/
	

				return true;
			}
		}
		else {
			/* Getting the length + direction of the constraints*/
			const int constraintSize = ChosenFaces.size();
			if (constraintSize > 0)
			{				
	
				/* UPdating the mesh information */
				///viewer.data().clear();
				///viewer.data().set_mesh(V, F);

				//printf("Size of the contraints vector is %d\n", constraintSize);
				if (fieldsType == FieldsType::VECTOR)
				{
					if (singularConstraint)
					{
						for (int i : ChosenFaces) {
							vectorFields.userSingularConstraints.push_back(F(i, 0));
							printf("Inserting %d as singularity constraints \n", F(i,0));
							vectorFields.addSingularityConstraints();
						}
					}
					else {
						constraintDir = vectorFields.FC.row(ChosenFaces[constraintSize - 1]) - vectorFields.FC.row(ChosenFaces[0]);
						viewer.data().add_edges(vectorFields.FC.row(ChosenFaces[0]), vectorFields.FC.row(ChosenFaces[0]) + constraintDir, Eigen::RowVector3d(1.0, 0.0, 0.1));
						vectorFields.pushNewUserConstraints(ChosenFaces[0], ChosenFaces[constraintSize - 1]);
						vectorFields.addHardConstraints();
					}

					//srand(time(NULL));
					//int vConst = rand() % V.rows();
					//vectorFields.userSingularConstraints.push_back(vConst);
					//printf("Inserting %d as singularity constraints \n", vConst);

					//printf("Pair [%d]->[%d] is inserted\n", ChosenFaces[0], ChosenFaces[constraintSize - 1]);

					/* Interactivity in reduced space */
					/////viewer.data().clear();
					/////viewer.data().set_mesh(V, F);
					///vectorFields.constructInteractiveConstraints();
					///vectorFields.setAndSolveInteractiveSystem(lambda);
					

					///vectorFields.constructInteractiveSingularities();
					///vectorFields.constructInteractiveConstraintsWithSingularities(viewer);
					
					//vectorFields.setupGlobalProblem(lambda);
					

					//vectorFields.setAndSolveUserSystem(lambda);
					vectorFields.setAndSolveInteractiveSystem(lambda);

					vectorFields.visualizeApproxResult(viewer);
					vectorFields.visualizeGlobalConstraints(viewer);
				}
				else if (fieldsType == FieldsType::NROSY)
				{
					constraintDir = nRoSyFields.FC.row(ChosenFaces[constraintSize - 1]) - nRoSyFields.FC.row(ChosenFaces[0]);
					viewer.data().add_edges(nRoSyFields.FC.row(ChosenFaces[0]), nRoSyFields.FC.row(ChosenFaces[0]) + constraintDir, Eigen::RowVector3d(1.0, 0.0, 0.1));
					nRoSyFields.pushNewUserConstraints(ChosenFaces[0], ChosenFaces[constraintSize - 1]);
					printf("Pair [%d]->[%d] is inserted\n", ChosenFaces[0], ChosenFaces[constraintSize - 1]);

					nRoSyFields.constructInteractiveConstraints();
					nRoSyFields.setAndSolveInteractiveSystem();

					nRoSyFields.visualizeConstrainedFields_Reduced(viewer);
					nRoSyFields.visualizeConstraints(viewer);
					
					NRoSy nRoSy_;
					nRoSyFields.convertRepVectorsToNRoSy(nRoSyFields.XfBar, nRoSy_);
					nRoSyFields.sendFieldsToMailSlot(nRoSy_);

					/* Full res*/
					//if (F.rows() < 50000)
					//{
					//	nRoSyFields.nRoSyFieldsDesignRef();
					//	nRoSyFields.visualizeConstraints(viewer);
					//	nRoSyFields.visualizeConstrainedFields(viewer);
					//}

					
					/* Reduced */
					///nRoSyFields.constructInteractiveConstraints();
					///nRoSyFields.nRoSyFieldsDesign_Reduced_Splines();
					///
					///nRoSyFields.visualizeConstraints(viewer);
					///nRoSyFields.visualizeConstrainedFields_Reduced(viewer);

				}
				else if (fieldsType == FieldsType::TENSOR)
				{
					
				}
				
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
	viewer.data().line_width = 1.5f; 
	//viewer.data().set_colors(Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333));
	return viewer.launch();

	/* Trick for remote desktop */
	//getchar();
	//return 1;
}