#include "VectorFields.h"

/* ====================== VISUALIZATION of IMPORTANT ELEMENTS ============================*/
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
	viewer.data().clear();
	viewer.data().set_mesh(V, F);
	Eigen::RowVector3d color1 = Eigen::RowVector3d(1.0, 0.1, 0.2);
	Eigen::RowVector3d color2 = Eigen::RowVector3d(0.0, 0.1, 1.0);
	visualize2DfieldsScaled(viewer, Xf.col(0), color1);
	//visualize2Dfields(viewer, Xf.col(0), color1);
	//visualize2Dfields(viewer, Xf.col(1), color2);	
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

	double lengthScale = 1.0*avgEdgeLength;
	for (int i = 0; i < F.rows(); i += 1)
	{
		Eigen::RowVector3d c;
		//c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		c = FC.row(i);
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
	}
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
	visualize2DfieldsScaled(viewer, BasisTemp.col(bId), color);

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

	printf("Showing the %d BasisTemp field\n", bId);
	visualize2DfieldsScaled(viewer, Basis.col(bId), color);

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

void VectorFields::visualizeApproxResult(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color;
	if (id == 0)
		color = Eigen::RowVector3d(1.0, 0.4, 0.4);
	else
		color = Eigen::RowVector3d(0.0, 0.4, 0.9);

	//cout << "Size of X_Lifted " << XFullDim.rows() << "x" << XFullDim.cols() << "." << endl; 
	//visualize2DfieldsNormalized(viewer, XFullDim.col(id), color);
	visualize2DfieldsScaled(viewer, XFullDim.col(id), color);
	//visualize2DfieldsRegular(viewer, XFullDim.col(id), color);
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
	for (int i = 0; i < globalConstraints.size(); i++) {
		Eigen::RowVector3d cc, g, v1, v2, v3;
		cc = FC.row(globalConstraints[i]);
		g = (A.block(3 * globalConstraints[i], 2 * globalConstraints[i], 3, 2) * c.block(2 * i, 0, 2, 1)).transpose();
		viewer.data().add_edges(cc, cc + g / g.norm()*avgEdgeLength, Eigen::RowVector3d(0.1, 0.1, 0.2));
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
	c = (V.row(F(faceID, 0)) + V.row(F(faceID, 1)) + V.row(F(faceID, 2))) / 3.0;
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

}