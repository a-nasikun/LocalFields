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
	visualize2DfieldsScaled(viewer, Xf, color, 5000);
}

void VectorFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, avgField;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
		blockLength(i) = g.norm();
		totalGrad += blockLength(i);
	}


	avgField = totalGrad / (double)F.rows();
	//double lengthScale = 1.0*avgEdgeLength / avgField;
	double lengthScale = 1.0*avgEdgeLength;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;

		//if (i == *(NeighRing[0].begin())){
		//	viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale*2.0, Eigen::RowVector3d(0.1, 0.1, 0.2));
		//	//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale*2.0, Eigen::RowVector3d(1.1, 0.1, 0.2));
		//}else 
		{
			viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
			//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, Eigen::RowVector3d(1.0, 0.1, 0.2));
			//viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, Eigen::RowVector3d(1.0, 0.1, 0.2));
		}

	}
}

void VectorFields::visualize2DfieldsNormalized(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, avgField;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
	}

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, color);
	}
}

void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());

	for (int i = 0; i < F.rows(); i += 1)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		//c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
		c = FC.row(i);
		//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
	}

	double lengthScale = 5.0*avgEdgeLength;
	for (int i = 0; i < F.rows(); i += 1)
	{
		Eigen::RowVector3d c;
		//c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		c = FC.row(i);
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
	}
}

void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const int &numFaces)
{
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);

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
	Eigen::MatrixXd VectorBlock(FaceToDraw.size(), F.cols());
	for (int i = 0; i < FaceToDraw.size(); i += 1)
	{
		c = FC.row(FaceToDraw[i]);												
		g = (A.block(3 * FaceToDraw[i], 2 * FaceToDraw[i], 3, 2) * field2D.block(2 * FaceToDraw[i], 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
	}
	//cout << "picking face to draw: done \n" << endl;

	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d f, h1, h2, e;
	Eigen::Vector2d v;
	Eigen::MatrixXd ALoc(3,2);
	for (int i = 0; i<FaceToDraw.size(); i += 1)
	{		
		c = FC.row(FaceToDraw[i]);
		//f = VectorBlock.row(i);
		v = field2D.block(2 * FaceToDraw[i], 0, 2, 1);
		ALoc = A.block(3 * FaceToDraw[i], 2 * FaceToDraw[i], 3, 2);
		f = (ALoc * v).transpose();
		h1 = (ALoc* (rotMat1*v)).transpose();
		h2 = (ALoc* (rotMat2*v)).transpose(); 
		e = c + f*lengthScale;
		viewer.data().add_edges(c, e, color);
		viewer.data().add_edges(e, e + h1*lengthScale / HEAD_RATIO, color);
		viewer.data().add_edges(e, e + h2*lengthScale / HEAD_RATIO, color);
	}
	//cout << "Drawing: done \n" << endl;
}

void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double &percent)
{

}

void VectorFields::visualize2DfieldsRegular(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, avgField;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
		blockLength(i) = g.norm();
		totalGrad += blockLength(i);
	}

	avgField = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
	}
}

/* Efficient visualization for sparse matrix */
void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::SparseMatrix<double> &Field2D, const int &idx, const Eigen::RowVector3d &color)
{	
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 5.0;
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
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

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

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

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
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	int bId = id;
	Eigen::RowVector3d color;
	if (id % 2 == 0) {
		color = Eigen::RowVector3d(1.0, 0.1, 0.2);
	}
	else {
		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	}

	if (id >= 2 * Sample.size()) {
		bId = 2 * Sample.size() - 1;
	}

	printf("Showing the %d BasisTemp field\n", bId);
	//visualize2DfieldsScaled(viewer, BasisTemp.col(bId), color);
	visualize2DfieldsScaled(viewer, BasisTemp, bId, color);

	Eigen::RowVector3d const c1 = (V.row(F(Sample[bId / 2], 0)) + V.row(F(Sample[bId / 2], 1)) + V.row(F(Sample[bId / 2], 2))) / 3.0;
	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
}

void VectorFields::visualizeBasisNormalized(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	int bId = id;
	Eigen::RowVector3d color;
	if (id % 2 == 0) {
		color = Eigen::RowVector3d(1.0, 0.1, 0.2);
	}
	else {
		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	}

	if (id >= 2 * Sample.size()) {
		bId = 2 * Sample.size() - 1;
	}

	printf("Showing the %d Basis field\n", bId);
	//visualize2DfieldsScaled(viewer, Basis.col(bId), color);
	visualize2DfieldsScaled(viewer, Basis, bId, color);

	Eigen::RowVector3d const c1 = (V.row(F(Sample[bId / 2], 0)) + V.row(F(Sample[bId / 2], 1)) + V.row(F(Sample[bId / 2], 2))) / 3.0;
	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
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

	visualize2DfieldsScaled(viewer, BasisSum.col(id), color);
	//visualize2DfieldsScaled(viewer, BasisSumN.col(id), color);
	//for (int i = 0; i < Sample.size(); i++) {
	//	Eigen::RowVector3d const c1 = (V.row(F(Sample[i], 0)) + V.row(F(Sample[i], 1)) + V.row(F(Sample[i], 2))) / 3.0;
	//	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
	//}
}

void VectorFields::visualizeApproxResult(igl::opengl::glfw::Viewer &viewer)
{
	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color;	
	color = Eigen::RowVector3d(0.9, 0.1, 0.1);

	//cout << "Size of X_Lifted " << XFullDim.rows() << "x" << XFullDim.cols() << "." << endl; 
	//visualize2DfieldsNormalized(viewer, XFullDim, color);
	visualize2DfieldsScaled(viewer, XFullDim, color, 5000);
	//visualize2DfieldsRegular(viewer, XFullDim, color);
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
	const double ARRAW_RATIO = 1.0; 
	for (int i = 0; i < globalConstraints.size(); i++) {
		Eigen::RowVector3d cc, g, v1, v2, v3;
		cc = FC.row(globalConstraints[i]);
		g = (A.block(3 * globalConstraints[i], 2 * globalConstraints[i], 3, 2) * c.block(2 * i, 0, 2, 1)).transpose();
		viewer.data().add_edges(cc, cc + ARRAW_RATIO*g*avgEdgeLength, Eigen::RowVector3d(0.1, 0.1, 0.2));
		viewer.data().add_points(cc, Eigen::RowVector3d(0.1, 0.1, 0.2));
	}
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
			//viewer.data().add_edges(c, c + basis.transpose(), Eigen::RowVector3d(0.5, 0.1, 0.6));
		}
	}
	igl::jet(z, false, vColor);
	//viewer.data().set_colors(vColor);

	for (int i : singularities) {
		viewer.data().add_points(V.row(i), Eigen::RowVector3d(0.1, 0.9, 0.3));
	}
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

void VectorFields::visualizeArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(arbField, true, vColor);
	viewer.data().set_colors(vColor);
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
	Eigen::VectorXd dom(F.rows());
	for (int i = 0; i < F.rows(); i++) dom(i) = 0.0;

	for (std::set<int>::iterator it = SubDomain.begin(); it != SubDomain.end(); ++it) {
		dom(*it) = 0.5;
		if (*it == 542) dom(*it) = 1.0;
	}

	for (std::set<int>::iterator it = Boundary.begin(); it != Boundary.end(); ++it) {
		dom(*it) = 0.25;
	}


	Eigen::MatrixXd vColor;
	igl::jet(dom, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeSamples(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::RowVector3d color(0.1, 0.1, 0.8);
	Eigen::RowVector3d c;

	for (int i : Sample) {
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_points(c, color);
	}
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