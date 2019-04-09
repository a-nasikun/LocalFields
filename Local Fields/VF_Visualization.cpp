#include "VectorFields.h"

/* ====================== VISUALIZATION of IMPORTANT ELEMENTS ============================*/

void VectorFields::selectFaceToDraw(const int& numFaces)
{
	/*Getting faces to draw, using farthest point sampling (could take some time, but still faster than drawing everything for huge mesh) */

	if (numFaces < F.rows())
	{
		FaceToDraw.resize(numFaces);
		Eigen::VectorXd D(F.rows());

		/* Initialize the value of D */
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		//srand(time(NULL));
		//FaceToDraw[0] = rand() % F.rows();
		FaceToDraw[0] = 0;

		for (int i = 1; i < numFaces; i++) {
			Eigen::VectorXi::Index maxIndex;
			computeDijkstraDistanceFaceForSampling(FaceToDraw[i - 1], D);
			D.maxCoeff(&maxIndex);
			FaceToDraw[i] = maxIndex;
		}
	}
	else
	{
		FaceToDraw.resize(F.rows());
		for (int i = 0; i < F.rows(); i++)
		{
			FaceToDraw[i] = i;
		}
	}
}

void VectorFields::visualizeMassMatrix(igl::opengl::glfw::Viewer &viewer, const MassMatrixToShow &type)
{
	Eigen::VectorXd z;
	Eigen::MatrixXd vColor;

	switch (type)
	{
	case MASS_MV:
		z.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			z(i) = MV.coeff(i, i);
		}
		break;

	case MASS_MVinv:
		z.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			z(i) = MVinv.coeff(i, i);
		}
		break;

	case MASS_MF2D:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF2D.coeff(2 * i + 1, 2 * i + 1);
		}
		break;

	case MASS_MF2Dinv:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF2Dinv.coeff(2 * i + 1, 2 * i + 1);
		}
		break;

	case MASS_MF3D:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF3D.coeff(3 * i + 2, 3 * i + 2);
		}
		break;

	case MASS_MF3Dinv:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF3Dinv.coeff(3 * i + 2, 3 * i + 2);
		}
		break;

	default:
		break;
	}

	igl::jet(z, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeGradient(igl::opengl::glfw::Viewer &viewer, const GradientToShow &type)
{
	switch (type)
	{
	case GRAD_2D:

		break;

	case GRAD_3D:

		break;

	default:
		break;
	}
}

void VectorFields::visualizeLocalFrames(igl::opengl::glfw::Viewer &viewer)
{
	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, e, f;
		e << A.coeffRef(3 * i, 2 * i), A.coeffRef(3 * i + 1, 2 * i), A.coeffRef(3 * i + 2, 2 * i);
		f << A.coeffRef(3 * i, 2 * i + 1), A.coeffRef(3 * i + 1, 2 * i + 1), A.coeffRef(3 * i + 2, 2 * i + 1);
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;

		// First basis
		viewer.data().add_edges(c, c + avgEdgeLength*e, Eigen::RowVector3d(1.0, 0.1, 0.2));
		// Second basis
		viewer.data().add_edges(c, c + avgEdgeLength*f, Eigen::RowVector3d(0.0, 0.1, 0.9));
	}
}

void VectorFields::visualizeApproximatedFields(igl::opengl::glfw::Viewer &viewer)
{
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);
	Eigen::RowVector3d color = Eigen::RowVector3d(0.1, 0.1, 0.9);	

	visualize2Dfields(viewer, Xf, color, 2.0, false);
	//visualize2Dfields(viewer, Xf, color, 2.0, true);
}

void VectorFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized)
{
	/* For Timing*/
	chrono::high_resolution_clock::time_point	t1, t2, te1, te2, ta1, ta2;
	chrono::duration<double>					duration, da, de;
	t1 = chrono::high_resolution_clock::now();
	//cout << "> Adding edges... ";

//<<<<<<< HEAD
//	for (int i = 0; i < F.rows(); i++)
//	{
//		Eigen::RowVector3d c;
//		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
//		viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, color);
//	}
//}
//void VectorFields::visualize2DfieldsNormalized(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const int &numFaces)
//{
//	/* Some constants for arrow drawing */
//	const double HEAD_RATIO = 5.0;
//	const double EDGE_RATIO = 10.0;
//
//=======
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 3.0;
	const double EDGE_RATIO = scale;
	double lengthScale = EDGE_RATIO*avgEdgeLength;
//>>>>>>> master

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	Eigen::SparseMatrix<double> MRot1(2 * FaceToDraw.size(), 2 * FaceToDraw.size()), MRot2(2 * FaceToDraw.size(), 2 * FaceToDraw.size());
	vector<Eigen::Triplet<double>> R1Triplet, R2Triplet;
	R1Triplet.reserve(2 * 2 * FaceToDraw.size());
	R2Triplet.reserve(2 * 2 * FaceToDraw.size());
	Eigen::MatrixXd FCLoc(FaceToDraw.size(), 3);

	/* Defining the rotation matrix (2-by-2) on the local frame */
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	for (int i = 0; i < FaceToDraw.size(); i++)
	{
		/* Rotation matrix for the first head */
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, rotMat1(0, 0)));
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, rotMat1(1, 0)));
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, rotMat1(0, 1)));
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, rotMat1(1, 1)));

		/* Rotation matrix for the second head */
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, rotMat2(0, 0)));
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, rotMat2(1, 0)));
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, rotMat2(0, 1)));
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, rotMat2(1, 1)));

		/* Getting the face center of selected faces */
		FCLoc.row(i) = FC.row(FaceToDraw[i]);
	}
	MRot1.setFromTriplets(R1Triplet.begin(), R1Triplet.end());
	MRot2.setFromTriplets(R2Triplet.begin(), R2Triplet.end());

	/* Getting the local data from the population of data */
	Eigen::SparseMatrix<double> ALoc(3 * FaceToDraw.size(), 2 * FaceToDraw.size());
	vector<Eigen::Triplet<double>> ATriplet;
	ATriplet.reserve(6 * FaceToDraw.size());
	Eigen::VectorXd fieldLoc(2 * FaceToDraw.size()), fields3D(3 * FaceToDraw.size()), rot1Field, rot2Field;
	Eigen::MatrixXd TFields(FaceToDraw.size(), F.cols()), Head1Fields(FaceToDraw.size(), F.cols()), Head2Fields(FaceToDraw.size(), F.cols());

	for (int i = 0; i < FaceToDraw.size(); i++)
	{
		/* Getting the selected ALoc from A */
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 0, A.coeff(3 * FaceToDraw[i] + 0, 2 * FaceToDraw[i] + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 0, A.coeff(3 * FaceToDraw[i] + 1, 2 * FaceToDraw[i] + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 0, A.coeff(3 * FaceToDraw[i] + 2, 2 * FaceToDraw[i] + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 1, A.coeff(3 * FaceToDraw[i] + 0, 2 * FaceToDraw[i] + 1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 1, A.coeff(3 * FaceToDraw[i] + 1, 2 * FaceToDraw[i] + 1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 1, A.coeff(3 * FaceToDraw[i] + 2, 2 * FaceToDraw[i] + 1)));

		/* Getting the selected face */
		fieldLoc.block(2 * i, 0, 2, 1) = field2D.block(2 * FaceToDraw[i], 0, 2, 1);
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());
	fields3D = ALoc * fieldLoc;

	/* The head of the arrows */
	rot1Field = MRot1*fieldLoc;
	rot1Field = ALoc * rot1Field;
	rot2Field = MRot2*fieldLoc;
	rot2Field = ALoc * rot2Field;

	/* Transform field to Matrix format */
	for (int i = 0; i < FaceToDraw.size(); i++)
	{
		TFields.row(i) = (fields3D.block(3 * i, 0, 3, 1)).transpose();
		Head1Fields.row(i) = (rot1Field.block(3 * i, 0, 3, 1)).transpose();
		Head2Fields.row(i) = (rot2Field.block(3 * i, 0, 3, 1)).transpose();
	}

	/* If user wants normalized fields, then so do it */
	if (normalized)
	{
		TFields.rowwise().normalize();
		Head1Fields.rowwise().normalize();
		Head2Fields.rowwise().normalize();
	}

	/* Draw the fields */
	viewer.data().add_edges(FCLoc, FCLoc + TFields*lengthScale, color);
	viewer.data().add_edges(FCLoc + TFields*lengthScale, FCLoc + TFields*lengthScale + Head1Fields*lengthScale / HEAD_RATIO, color);
	viewer.data().add_edges(FCLoc + TFields*lengthScale, FCLoc + TFields*lengthScale + Head2Fields*lengthScale / HEAD_RATIO, color);

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::visualize2DfieldsNormalized(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const int &numFaces)
{
	visualize2Dfields(viewer, field2D, color, 2.0, true);	
}

void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double &scale)
{
	visualize2Dfields(viewer, field2D, color, scale);	
}

void VectorFields::visualize2DfieldsRegular(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	visualize2Dfields(viewer, field2D, color, 1.0);
}

/* Efficient visualization for sparse matrix */
void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::SparseMatrix<double> &Field2D, const int &idx, const Eigen::RowVector3d &color)
{	
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 5.0;
	//const double EDGE_RATIO = 0.075;			// for the eigenfields patch
	//const double EDGE_RATIO = 0.125;
	const double EDGE_RATIO = 2;

	/*Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	Eigen::RowVector3d c, g; 
	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d f, h1, h2, e;
	Eigen::Vector2d v;
	Eigen::MatrixXd ALoc(3, 2);

	for (Eigen::SparseMatrix<double>::InnerIterator it(Field2D, idx); it; ++it) {
		if (it.row() % 2 == 1) continue;
		v << it.value(), Field2D.coeff(it.row() + 1, idx);
		ALoc = A.block(3 * it.row()/2, 2 * it.row()/2, 3, 2);

		c = FC.row(it.row()/2);
		f = (ALoc * v).transpose();
		h1 = (ALoc* (rotMat1*v)).transpose();
		h2 = (ALoc* (rotMat2*v)).transpose();
		e = c + f*lengthScale;
		viewer.data().add_edges(c, e, color);
		viewer.data().add_edges(e, e + h1*lengthScale / HEAD_RATIO, color);
		viewer.data().add_edges(e, e + h2*lengthScale / HEAD_RATIO, color);
	}

	
}

void VectorFields::visualize3Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field3D, const Eigen::RowVector3d &color)
{
	// Reset
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);

	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = field3D.block(3 * i, 0, 3, 1).transpose();
		VectorBlock.row(i) = g;
		blockLength(i) = g.norm();
		totalGrad += blockLength(i);
	}

	//averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 0.2*avgEdgeLength; // / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, Eigen::RowVector3d(1.0, 0.1, 0.2));
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
		//viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, Eigen::RowVector3d(1.0, 0.1, 0.2));
	}
}

void VectorFields::visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);

	int bId = id;
	Eigen::RowVector3d color;
	if (id % 2 == 0) {
		color = Eigen::RowVector3d(1.0, 0.1, 0.2);
		//color = Eigen::RowVector3d(0.6, 0.2, 0.2);
	}
	else {
		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	}

	if (id >= 2 * Sample.size()) {
		bId = 2 * Sample.size() - 1;
	}

	printf("Showing the %d BasisTemp field (Sample=%d) \n", bId, Sample[id/2]);
	visualize2DfieldsScaled(viewer, BasisTemp, bId, color);

	Eigen::RowVector3d const c1 = (V.row(F(Sample[bId / 2], 0)) + V.row(F(Sample[bId / 2], 1)) + V.row(F(Sample[bId / 2], 2))) / 3.0;
	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
}

void VectorFields::visualizeBasisNormalized(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);

	int bId = id;
	Eigen::RowVector3d color, color2;
	if (id % 2 == 0) {
		//color = Eigen::RowVector3d(1.0, 0.1, 0.2);
		color  = Eigen::RowVector3d(0.6, 0.2, 0.2);
		color2 = Eigen::RowVector3d(0.4, 0.8, 0.8);
	}
	else {
		//color = Eigen::RowVector3d(0.0, 0.1, 1.0);
		color  = Eigen::RowVector3d(0.2, 0.2, 0.6);
		color2 = Eigen::RowVector3d(0.8, 0.8, 0.4);
	}

	//if (id >= 2 * Sample.size()) {
	//	bId = 2 * Sample.size() - 1;
	//}

	

	//Eigen::RowVector3d const c1 = (V.row(F(Sample[bId / 2], 0)) + V.row(F(Sample[bId / 2], 1)) + V.row(F(Sample[bId / 2], 2))) / 3.0;
	//viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));

	//if (id < 2)
	//{
	//	printf(">  Showing local eigenfields [%d] \n", id);
	//	visualize2Dfields(viewer, eigFieldsLocal.col(id), color2, 0.6, false);
	//} 
	//else
	//{
		printf("Showing the %d Basis field\n", bId);
		visualize2DfieldsScaled(viewer, Basis, bId, color);
	//}
}

void VectorFields::visualizeBasisSum(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color;
	if (id == 0)
		color = Eigen::RowVector3d(1.0, 0.4, 0.4);
	else
		color = Eigen::RowVector3d(0.0, 0.4, 0.9);

	visualize2DfieldsScaled(viewer, BasisSum.col(id), color,1.0);
}

void VectorFields::visualizeApproxResult(igl::opengl::glfw::Viewer &viewer)
{
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color;	
	color = Eigen::RowVector3d(0.9, 0.1, 0.1);

	//cout << "Size of X_Lifted " << XFullDim.rows() << "x" << XFullDim.cols() << "." << endl;
	//visualize2Dfields(viewer, XFullDim, color, 2, false);
	visualize2Dfields(viewer, XFullDim, color, 2, true);
	//cout << "XFULL approx. \n " << XFullDim.block(0, 0, 100, 1) << endl; 
}

void VectorFields::visualizeUserConstraints(igl::opengl::glfw::Viewer &viewer)
{
	for (int i = 0; i < userConstraints.size(); i++) {
		Eigen::RowVector3d c, g, v1, v2, v3;
		//c = (V.row(F(userConstraints[i], 0)) + V.row(F(userConstraints[i], 1)) + V.row(F(userConstraints[i], 2))) / 3.0;														// first vertex of each face [NOT GOOD]	
		c = FC.row(userConstraints[i]);
		g = (A.block(3 * userConstraints[i], 2 * userConstraints[i], 3, 2) * cBar.block(2 * i, 0, 2, 1)).transpose();
		viewer.data().add_edges(c, c + g / g.norm()*avgEdgeLength, Eigen::RowVector3d(0.1, 0.1, 0.2));
		viewer.data().add_points(c, Eigen::RowVector3d(0.1, 0.1, 0.2));
	}
}

void  VectorFields::visualizeGlobalConstraints(igl::opengl::glfw::Viewer &viewer)
{
	/* ORIGINAL + OVERLAY on 2nd Mesh */
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 3.0;
	const double ARRAW_RATIO = 4.0; 
	const double EDGE_RATIO = 4.0;
	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d color(0.0, 0.0, 0.2);
	
	viewer.selected_data_index = 1; 
	viewer.data().line_width = 5.0;
	viewer.data().point_size = 5.0;
	viewer.data().show_lines = false;

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);
	
	Eigen::MatrixXd ALoc(3, 2);
	for (int i = 0; i < globalConstraints.size(); i++) {
		Eigen::RowVector3d cc, g, h1, h2, v1, v2, v3, e, f, n;
		Eigen::Vector2d v; 
		cc = FC.row(globalConstraints[i]);
		n = NF.row(globalConstraints[i]);
		n *= (avgEdgeLength/10.0); 
		ALoc = A.block(3 * globalConstraints[i], 2 * globalConstraints[i], 3, 2);
		v = c.block(2 * i, 0, 2, 1);
		g = (ALoc * v).transpose();
		cc += n; 
		f = cc + g*lengthScale;
		viewer.data().add_edges(cc, f, color);
		viewer.data().add_points(cc, color);

		h1 = (ALoc * (rotMat1*c.block(2 * i, 0, 2, 1))).transpose();
		h2 = (ALoc * (rotMat2*c.block(2 * i, 0, 2, 1))).transpose();

		viewer.data().add_edges(f, f + h1*lengthScale/HEAD_RATIO, color);
		viewer.data().add_edges(f, f + h2*lengthScale/HEAD_RATIO, color);
	}
	 
	viewer.selected_data_index = 0; 	

}

void VectorFields::visualizeSingularitiesConstraints(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;

	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}

	for (int id = 0; id < SingNeighCC.size(); id++) {
		const double diff = 0.6 / (double)SingNeighCC[id].size();
		for (int i = 0; i < SingNeighCC[id].size(); i++) {
			z(SingNeighCC[id][i]) = 0.3 + i*diff;

			Eigen::Vector3d basis = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i], 3, 1);
			basis *= avgEdgeLength; 
			Eigen::RowVector3d c = FC.row(SingNeighCC[id][i]);
			viewer.data().add_edges(c, c + basis.transpose(), Eigen::RowVector3d(0.5, 0.1, 0.6));
		}
	}

	for (int id = 0; id < mappedBasis.size(); id++) {
		for (int i = 0; i < mappedBasis[id].size(); i++) {
			Eigen::MatrixXd Map2D(3, 2);
			Map2D = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i], 3, 2);
			Eigen::Vector3d oldBasis = Map2D.col(0);
			Eigen::Vector3d newBasis = Map2D * mappedBasis[id][i];

			Eigen::RowVector3d c = FC.row(SingNeighCC[id][i]);
			viewer.data().add_edges(c, c + avgEdgeLength * oldBasis.transpose(), Eigen::RowVector3d(0.0, 0.1, 0.3));
			viewer.data().add_edges(c, c + avgEdgeLength * Map2D.col(1).transpose(), Eigen::RowVector3d(0.0, 0.3, 0.1));

			viewer.data().add_edges(c, c + avgEdgeLength * newBasis.transpose(), Eigen::RowVector3d(0.0, 0.0, 0.9));
			newBasis = Map2D * mappedBasis2[id][i];
			viewer.data().add_edges(c, c + avgEdgeLength * newBasis.transpose(), Eigen::RowVector3d(0.0, 0.9, 0.0));
		}
	}

	igl::jet(z, false, vColor);
	//viewer.data().set_colors(vColor);

	for (int i : singularities) {
		viewer.data().add_points(V.row(i), Eigen::RowVector3d(0.1, 0.9, 0.3));
	}
}
void VectorFields::write2DFieldsForVTK(const Eigen::VectorXd &field2D, const string& dataName, const string& filename)
{
	ofstream file(filename);

	cout << "Trying to write files \n";

	if (file.is_open())
	{
		/* Print headers */
		file << "# vtk DataFile Version 2.0\n";
		file << dataName << "\n";
		file << "ASCII\n";
		file << "DATASET POLYDATA\n";

		/* Print vertices */
		file << "POINTS " << V.rows() << " double\n";
		for (int i = 0; i < V.rows(); i++)
		{
			file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
		}
		
		/* Print faces */
		file << "POLYGONS " << F.rows() << " " << F.rows() * (F.cols() + 1) << "\n";
		for (int i = 0; i < F.rows(); i++)
		{
			file << F.cols() << " " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\n";
		}

		cout << "Converting to world space \n";
		/* Print eigenfields*/
		Eigen::VectorXd field3D;
		Eigen::MatrixXd Fields3D;
		double const scale = 1.0;
		//field3D = A * field2D;
		field3D = A * Xf;

		//cout << "FACE-based\n";
		///* FACE-base DATA */
		//file << "CELL_DATA " << F.rows() << "\n";
		//file << "VECTORS " << dataName << "_faceBased double\n";
		//for (int i = 0; i < F.rows(); i++)
		//{
		//	Eigen::Vector3d v;
		//	v << field3D(3 * i + 0), field3D(3 * i + 1), field3D(3 * i + 2);
		//	v.normalize();
		//	//file << scale*field3D(3 * i + 0) << " " << scale*field3D(3 * i + 1) << " " << scale*field3D(3 * i + 2) << "\n";
		//	file << v(0) << " " << v(1) << " " << v(2) << "\n";
		//}
		 
		//file << "POINTS_DATA" << V.rows() << "\n";
		//file << "NORMALS vertex_normals float \n";
		//Eigen::MatrixXd NV;
		//igl::per_vertex_normals(V, F, NV);
		//for (int i = 0; i < V.rows(); i++)
		//{
		//	file << NV(i, 0) << " " << NV(i, 1) << " " << NV(i, 2) << "\n";
		//}

		//file << "NORMALS face_normals float\n";
		//for (int i = 0; i < F.rows(); i++)
		//{
		//	file << NF(i, 0) << " " << NF(i, 1) << " " << NF(i, 2) << "\n";
		//}

		cout << "VERTEX-based Mode\n"; 
		/* VERTEX-based DATA */
		Eigen::VectorXd vFields3D;
		vFields3D.setZero(3 * V.rows());
		Eigen::VectorXi VNumNeighFaces;
		VNumNeighFaces.setZero(V.rows());
		Eigen::VectorXd VSumAngles;
		VSumAngles.setZero(V.rows());

		for (int i = 0; i < F.rows(); i++)
		{
			//vFields3D.block(3 * F(i, 0), 0, 3, 1) += doubleArea(i)*field3D.block(3 + i, 0, 3, 1);
			//vFields3D.block(3 * F(i, 1), 0, 3, 1) += doubleArea(i)*field3D.block(3 + i, 0, 3, 1);
			//vFields3D.block(3 * F(i, 2), 0, 3, 1) += doubleArea(i)*field3D.block(3 + i, 0, 3, 1);
			vFields3D.block(3 * F(i, 0), 0, 3, 1) += field3D.block(3 + i, 0, 3, 1);
			vFields3D.block(3 * F(i, 1), 0, 3, 1) += field3D.block(3 + i, 0, 3, 1);
			vFields3D.block(3 * F(i, 2), 0, 3, 1) += field3D.block(3 + i, 0, 3, 1);
			VNumNeighFaces(F(i, 0)) += 1;
			VNumNeighFaces(F(i, 1)) += 1;
			VNumNeighFaces(F(i, 2)) += 1;
			VSumAngles(F(i, 0)) += doubleArea(i);
			VSumAngles(F(i, 1)) += doubleArea(i);
			VSumAngles(F(i, 2)) += doubleArea(i);
		}
		
		for (int i = 0; i < V.rows(); i++)
		{
			vFields3D.block(3 * i, 0, 3, 1) /= (double)VNumNeighFaces(i);
			//vFields3D.block(3 * i, 0, 3, 1) /= (double)VSumAngles(i);
			vFields3D.block(3 * i, 0, 3, 1) = vFields3D.block(3 * i, 0, 3, 1).normalized();
		}
		
		file << "POINT_DATA " << V.rows() << "\n";
		file << "VECTORS " << dataName << "_VertexBased double\n";
		for (int i = 0; i < V.rows(); i++)
		{
			file << scale*vFields3D(3 * i + 0) << " " << scale*vFields3D(3 * i + 1) << " " << scale*vFields3D(3 * i + 2) << "\n";
		}

		file << "NORMALS vertex_normals float \n";
		Eigen::MatrixXd NV;
		igl::per_vertex_normals(V, F, NV);
		for (int i = 0; i < V.rows(); i++)
		{
			file << NV(i, 0) << " " << NV(i, 1) << " " << NV(i, 2) << "\n";
		}
		cout << "Closing files\n";
		file.close();
	}
}

void VectorFields::visualizeSmoothing(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& v)
{
	visualize2DfieldsScaled(viewer, v, Eigen::RowVector3d(0.1, 0.9, 0.3), 1.0);
}

void VectorFields::visualizeCurvatureTensor(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Visualizing the curvature fields\n";
	// Average edge length for sizing
	const Eigen::RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8), green(0.2, 0.8, 0.2);

	/* From tensor */
	Eigen::VectorXd maxCurvField = CurvatureTensorField2D.col(0);
	Eigen::VectorXd minCurvField = CurvatureTensorField2D.col(1);

	/* FACE BASED */	
	visualize2DfieldsScaled(viewer,  maxCurvField, red,  0.3);
	visualize2DfieldsScaled(viewer, -maxCurvField, red,  0.3);
	visualize2DfieldsScaled(viewer,  minCurvField, blue, 0.3);
	visualize2DfieldsScaled(viewer, -minCurvField, blue, 0.3);
}

/* ====================== VISUALIZATION for TESTING ELEMENTS ============================*/
void VectorFields::visualizeFaceNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx) {
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;

	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}
	z(idx) = 1.0;

	//for (int i = 0; i < AdjMF3N.cols(); i++) {
	//	if(i==0)
	//	z(AdjMF3N(idx, i)) = 0.5; 
	//}

	for (std::set<int, double>::iterator it = AdjMF2Ring[idx].begin(); it != AdjMF2Ring[idx].end(); ++it) {
		z(*it) = 0.5f;
	}
	
	// Visualizing the COLOR
	igl::jet(z, true, vColor);
	printf("Size of color: %dx%d\n", vColor.rows(), vColor.cols());
	viewer.data().set_colors(vColor);

	// Visualilzing the EDGE
	for (int i = 0; i<F.cols(); i++) {
		if (i == 0)
			viewer.data().add_edges(V.row(EdgePairMatrix(idx, 2 * i)), V.row(EdgePairMatrix(idx, 2 * i + 1)), Eigen::RowVector3d(0.2, 0.1, 0.2));
	}
}

void VectorFields::visualizeVertexFacesNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx)
{
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;

	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}

	for (std::set<VtoFPair>::iterator it = VFNeighFull[idx].begin(); it != VFNeighFull[idx].end(); ++it) {
		z(it->fId) = 0.5;
	}

	igl::jet(z, true, vColor);
	viewer.data().set_colors(vColor);

	viewer.data().add_points(V.row(idx), Eigen::RowVector3d(0.3, 0.3, 0.9));
}


void VectorFields::visualizeDijkstra(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(arbField, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeEigenfields(igl::opengl::glfw::Viewer &viewer, int i)
{
	Eigen::VectorXd eigfields = eigFieldFull2D.col(i);
	double l2norm = eigfields.transpose()*MF2D*eigfields;

	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	/* Define some colors */
	Eigen::RowVector3d purple(136.0 / 255.0, 86.0 / 255.0, 167.0 / 255.0);
	Eigen::RowVector3d blue(0.1, 0.0, 0.9);

	/* Visualizing the fields*/
	//visualize2Dfields(viewer, eigfields, blue, 3, true);
	visualize2Dfields(viewer, eigfields, purple, 3, false);
}

void VectorFields::visualizeApproxEigenfields(igl::opengl::glfw::Viewer &viewer, int i, int iRef)
{
	Eigen::VectorXd eigfields = Basis*eigFieldReduced2D.col(i);
	double l2norm = eigfields.transpose()*MF2D*eigfields;
	eigfields /= sqrt(l2norm);
	l2norm = eigfields.transpose()*MF2D*eigfields;

	/* Inverse the sign to match the reference */
	if (eigFieldFull2D.size() > 0)
	{
		if ((eigfields.transpose()*MF2D*eigFieldFull2D.col(iRef)) < 0) eigfields *= -1; 

		/* Computing the L2norm error */	
		Eigen::VectorXd diff = eigFieldFull2D.col(iRef) - eigfields;
		double norm1 = diff.transpose()*MF2D*diff;
		double norm2 = eigFieldFull2D.col(iRef).transpose()*MF2D*eigFieldFull2D.col(iRef);
		l2norm = sqrt(norm1 / norm2);
		printf("L2 norm of %d approx is: %.5f\n", i, l2norm);
	}


	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	/* Define some colors */
	Eigen::RowVector3d purple(136.0 / 255.0, 86.0 / 255.0, 167.0 / 255.0);
	Eigen::RowVector3d red(0.9, 0.1, 0.1);	
	
	//visualize2Dfields(viewer, eigfields, red, 3, false);
	visualize2Dfields(viewer, eigfields, red, 3, true);
}

void VectorFields::visualizeArbField(igl::opengl::glfw::Viewer &viewer)
{
	/* as a Function */
	//Eigen::MatrixXd vColor;
	//igl::jet(arbField, true, vColor);
	//viewer.data().set_colors(vColor);

	/* As a Vector field */
	Eigen::RowVector3d color(0.8, 0.4, 0.2);
	visualize2DfieldsScaled(viewer, arbField2D, color, 1.0);
}

void VectorFields::visualizeRandomFace(igl::opengl::glfw::Viewer &viewer, const int &faceID)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d c;
	c = FC.row(faceID);
	viewer.data().add_points(c, Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeDijkstraFace(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(dijkstraFace, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeSubdomain(igl::opengl::glfw::Viewer &viewer)
{
	/* Color map => JET */
	//Eigen::VectorXd dom(F.rows());
	//for (int i = 0; i < F.rows(); i++) dom(i) = 0.0;
	//
	//for (std::set<int>::iterator it = SubDomain.begin(); it != SubDomain.end(); ++it) {
	//	dom(*it) = 0.3;
	//	if (*it == Sample[0]) dom(*it) = 1.0;
	//}
	//
	//for (std::set<int>::iterator it = Boundary.begin(); it != Boundary.end(); ++it) {
	//	dom(*it) = 0.7;
	//}


	Eigen::MatrixXd vColor;
	//igl::jet(dom, true, vColor);
	//igl::jet(localSystem, false, vColor);

	/* My Own */
	vColor.resize(F.rows(), 3);
	// 0:background => eeeeee; 0.3:selected region; 0.7:boundary
	for (int i = 0; i < F.rows(); i++)
	{
		if (localSystem(i) < 0.1)
		{
			vColor.row(i) = Eigen::RowVector3d(0.93333333, 0.93333333, 0.9333333);
		}
		else if (localSystem(i) > 0.6)
		{
			vColor.row(i) = Eigen::RowVector3d(0.96078431372, 0.36470588235, 0.2431372549);
		}
		else
		{
			//vColor.row(i) = Eigen::RowVector3d(1, 0.88235294117, 0.77647058823);
			//vColor.row(i) = Eigen::RowVector3d(0.89803921568, 0.94901960784, 0.78823529411);
			vColor.row(i) = Eigen::RowVector3d(186.0/255.0, 212.0/255.0, 170.0/255.0);
			
		}
	}

	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeSamples(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::RowVector3d color(0.0, 0.8, 0.0);
	Eigen::RowVector3d c;

	/* Point based */
	//for (int i : Sample) {
	//	c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	//	viewer.data().add_points(c, color);
	//}

	/* Color based */
	Eigen::MatrixXd FColor;
	igl::jet(-sampleDistance, true, FColor);
	viewer.data().set_colors(FColor);
}

void VectorFields::visualizeSharedEdges(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::RowVector3d color(0.8, 0.1, 0.0);
	Eigen::RowVector3d color2(0.8, 0.1, 0.1);

	for (int i = 0; i < sharedEdgesVect.size(); i++) {
		for (int j = 0; j < sharedEdgesVect[i].size(); j++) {
			if (j % 2 == 1) continue;

			double k = 1.0 / (double)sharedEdgesVect[i].size() * (double) j; 
			Eigen::RowVector3d color3(0.1, 0.1, k);
			viewer.data().add_edges(V.row(sharedEdgesVect[i][j]), V.row(sharedEdgesVect[i][j + 1]), color3);
			viewer.data().add_points(V.row(sharedEdgesVect[i][j]), color3);
		}
	}

	// Visualizing the basis
	Eigen::MatrixXd ALoc(3, 2);
	for (int i = 0; i < SingNeighCC.size(); i++)
	{
		for (int j = 0; j < SingNeighCC[i].size(); j++)
		{
			ALoc = A.block(3 * SingNeighCC[i][j], 2 * SingNeighCC[i][j], 3, 2);
			viewer.data().add_edges(FC.row(SingNeighCC[i][j]), FC.row(SingNeighCC[i][j]) + avgEdgeLength*ALoc.col(0).transpose(), Eigen::RowVector3d(0.0, 0.9, 0.3));
		}
	}

}

void VectorFields::visualizeLocalSubdomain(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd fColor; 
	igl::jet(localSystem, false, fColor);
	viewer.data().set_colors(fColor);
}

void VectorFields::visualizeParallelTransportPath(igl::opengl::glfw::Viewer &viewer)
{
	const double colorScale = 0.5 / (double) PTpath.size();
	Eigen::MatrixXd fColor;
	Eigen::VectorXd z(F.rows());

	for (int i = 0; i < F.rows(); i++) z(i) = 0.0;

	for (int i = 0; i < PTpath.size(); i++)
	{
		z(PTpath[i]) = 0.3 + i*colorScale;
		//viewer.data().add_edges()
	}

	igl::jet(z, false, fColor);
	viewer.data().set_colors(fColor);
}

void VectorFields::visualizeParallelTransport(igl::opengl::glfw::Viewer &viewer)
{
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 3.0;
	const double EDGE_RATIO = 1.0;

	/*Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	Eigen::RowVector3d c, g;
	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d f, h1, h2, e;
	Eigen::Vector2d v;
	Eigen::MatrixXd ALoc(3, 2);
	Eigen::RowVector3d color(0.8, 0.0, 0.3);
	for (int i = 0; i < parallelTransport.size(); i++) {
		v = parallelTransport[i];
		//v << it.value(), Field2D.coeff(it.row() + 1, idx);
		ALoc = A.block(3 * PTpath[i], 2 * PTpath[i], 3, 2);

		c = FC.row(PTpath[i]);
		f = (ALoc * v).transpose();
		h1 = (ALoc* (rotMat1*v)).transpose();
		h2 = (ALoc* (rotMat2*v)).transpose();
		e = c + f*lengthScale;
		viewer.data().add_edges(c, e, color);
		viewer.data().add_edges(e, e + h1*lengthScale / HEAD_RATIO, color);
		viewer.data().add_edges(e, e + h2*lengthScale / HEAD_RATIO, color);



		//f = (A.block(3 * PTpath[i], 2 * PTpath[i], 3, 2)*parallelTransport[i]).transpose();
		//viewer.data().add_edges(c, c + avgEdgeLength*f, color);

		//viewer.data().add_edges(V.row(PTsharedEdges[2*i+0]), V.row(PTsharedEdges[2 * i + 1]), color);
	}

	// Showing basis
	/*
	Eigen::Vector3d basis1, basis2;
	Eigen::RowVector3d color2(0.1, 0.2, 0.7);
	for (int i = 0; i < parallelTransport.size(); i++) {
	basis1 = A.block(3 * PTpath[i], 2 * PTpath[i], 3, 1);
	c = FC.row(PTpath[i]);
	viewer.data().add_edges(c, c+avgEdgeLength*basis1.transpose(), color2);
	}
	*/

	// Add points -> init of the path (target of Dijkstra)
	viewer.data().add_points(FC.row(PTpath[0]), Eigen::RowVector3d(0.0, 0.9, 0.1));
}
void VectorFields::visualizeCurveConstraints(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::VectorXd func(F.rows());
	Eigen::MatrixXd FColor(F.rows(), 3);
	const double colorDiff = 0.5 / (double) (curvesConstraints.size());

	for (int i = 0; i < F.rows(); i++) func(i) = 0.0;

	/*Draw the color*/
	for (int i = 0; i < curvesConstraints.size(); i++)
	{
		for (int j = 0; j < curvesConstraints[i].size(); j++)
		{
			func(curvesConstraints[i][j]) = (i+4)*colorDiff;
		}
	}


	igl::jet(func, false, FColor);
	viewer.data().set_colors(FColor);
}

void VectorFields::visualizeSoftConstraints(igl::opengl::glfw::Viewer &viewer)
{
	/* Color */
	Eigen::RowVector3d color(0.1, 0.1, 0.1);

	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 5.0;
	const double EDGE_RATIO = 1.0;

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d f, h1, h2, e, c;
	Eigen::Vector2d v;
	Eigen::MatrixXd ALoc(3, 2);
	int face;
	for (int i = 0; i < constraintVect2D.size(); i++)
	{
		for (int j = 0; j < constraintVect2D[i].size(); j++)
		{
			face = curvesConstraints[i][j];
			c = FC.row(face);
			//f = VectorBlock.row(i);
			v = constraintVect2D[i][j];
			ALoc = A.block(3 * face, 2 * face, 3, 2);
			f = (ALoc * v).transpose();
			h1 = (ALoc* (rotMat1*v)).transpose();
			h2 = (ALoc* (rotMat2*v)).transpose();
			e = c + f*lengthScale;
			//cout << "c: " << c << "e: " << e << endl; 
			viewer.data().add_edges(c, e, color);
			viewer.data().add_edges(e, e + h1*lengthScale / HEAD_RATIO, color);
			viewer.data().add_edges(e, e + h2*lengthScale / HEAD_RATIO, color);
		}
	}
}

void VectorFields::visualize1FieldOnCenter(igl::opengl::glfw::Viewer &viewer, const bool& even)
{
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 5.0;
	const double EDGE_RATIO = 1.0;

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	/* Drawing faces */
	Eigen::RowVector3d c, g;	

	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d f, h1, h2, e, color;
	Eigen::Vector2d v;
	Eigen::MatrixXd ALoc(3, 2);

	if (even)	color = Eigen::RowVector3d(1.0, 0.1, 0.0);
	else		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	
	c = FC.row(Sample[0]);
	if(even)	v << 1.0, 0.0;
	else 		v << 0.0, 1.0;
	ALoc = A.block(3 * Sample[0], 2 * Sample[0], 3, 2);
	f = (ALoc * v).transpose();
	h1 = (ALoc* (rotMat1*v)).transpose();
	h2 = (ALoc* (rotMat2*v)).transpose();
	e = c + f*lengthScale;
	viewer.data().add_edges(c, e, color);
	viewer.data().add_edges(e, e + h1*lengthScale / HEAD_RATIO, color);
	viewer.data().add_edges(e, e + h2*lengthScale / HEAD_RATIO, color);
	
}


void VectorFields::visualizePatchDijkstra(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::VectorXd distColor;
	distColor.setZero(F.rows());
	
	for (int i = 0; i < F.rows(); i++)
	{
		if (patchDijkstraDist(i) < numeric_limits<double>::infinity())
		{
			distColor(i) = patchDijkstraDist(i);
		}
	}

	Eigen::MatrixXd FColor;
	igl::jet(distColor, false, FColor);
	viewer.data().set_colors(FColor);
}


void VectorFields::visualizeAreaOfLaplaceConstraint(igl::opengl::glfw::Viewer &viewer)
{
	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			viewer.data().add_points(FC.row(floor(it.col()/2)), Eigen::RowVector3d(0.1, 0.1, 0.4));
			
		}
	}
}