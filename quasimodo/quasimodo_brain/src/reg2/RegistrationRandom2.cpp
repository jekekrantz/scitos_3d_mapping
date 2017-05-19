#include "RegistrationRandom2.h"
#include <iostream>
#include <fstream>
#include <omp.h>
#include <algorithm>

namespace reglib
{

RegistrationRandom2::RegistrationRandom2(unsigned int steps){
	only_initial_guess		= false;
	visualizationLvl		= 0;

	initTransform = Eigen::Matrix4d::Identity();

	Eigen::Affine3d startrot = Eigen::Affine3d::Identity();
	startrot = Eigen::AngleAxisd(-30*2*M_PI/360.0, Eigen::Vector3d::UnitX());
	initTransform = startrot.matrix();

	nr_src_target = 300000;
	nr_dst_target = 300000;

	src_nr_arraypoints = 0;
	src_arraypoints = 0;
	src_trees3d = 0;
	src_a3d = 0;

	dst_nr_arraypoints = 0;
	dst_arraypoints = 0;
	dst_trees3d = 0;
	dst_a3d = 0;

//    dfunc_ = new DistanceWeightFunction2();
//    dfunc_->f = THRESHOLD;
//    dfunc_->p = 0.005;

//    nfunc_ = new DistanceWeightFunction2();
//    nfunc_->f = THRESHOLD;
//    nfunc_->p = 0.50;

	steprx		= stepry	= steprz	= steps;
	start_rx	= start_ry	= start_rz	= 0;
	stop_rx		= stop_ry	= stop_rz	= 2.0 * M_PI;// * double(steps)/double(steps+1);

	steptx		= stepty	= steptz	= 1;
	start_tx	= start_ty	= start_tz	= 0;
	stop_tx		= stop_ty	= stop_tz	= 0;

		GeneralizedGaussianDistribution * ggddfunc	= new GeneralizedGaussianDistribution(true,false);
		ggddfunc->nr_refineiters		= 1;
		ggddfunc->debugg_print		= false;

		DistanceWeightFunction2PPR3 * dfuncTMP = new DistanceWeightFunction2PPR3(ggddfunc);
		dfunc_ = dfuncTMP;
		dfuncTMP->bidir					= false;
		dfuncTMP->zeromean				= true;
		dfuncTMP->maxp					= 0.9999;
////		dfuncTMP->maxd					= 1.0;//dstdval*10;
////		dfuncTMP->histogram_size		= 1000;
////		dfuncTMP->fixed_histogram_size	= false;
////		dfuncTMP->startreg				= 0.000;
////		dfuncTMP->max_under_mean		= false;
////		dfuncTMP->debugg_print			= false;
////		dfuncTMP->maxp					= 0.9999;
////		dfuncTMP->maxd					= 1.0;//dstdval*10;
////		dfuncTMP->histogram_size		= 1000;
////		dfuncTMP->fixed_histogram_size	= false;
////		dfuncTMP->startmaxd				= dfuncTMP->maxd;
////		dfuncTMP->starthistogram_size	= dfuncTMP->histogram_size;
////		dfuncTMP->blurval				= 0.0;
////		dfuncTMP->maxnoise				= dstdval;
////		dfuncTMP->compute_infront		= true;
////		dfuncTMP->ggd					= true;
////		dfuncTMP->reset();
////		dfunc->computeModel(dvec);

		GeneralizedGaussianDistribution * ggdnfunc	= new GeneralizedGaussianDistribution(true,true);
		ggdnfunc->nr_refineiters		= 10;
////		ggdnfunc->power					= 1.0;
////		ggdnfunc->costpen				= -1;
////		ggdnfunc->debugg_print			= false;
		DistanceWeightFunction2PPR3 * nfuncTMP		= new DistanceWeightFunction2PPR3(ggdnfunc);
		nfunc_ = nfuncTMP;
		nfuncTMP->bidir					= false;
		nfuncTMP->zeromean				= true;
		nfuncTMP->maxp					= 0.9999;
		nfuncTMP->maxd					= 2.0;
		nfuncTMP->histogram_size		= 200;
		nfuncTMP->fixed_histogram_size	= true;
		nfuncTMP->startmaxd				= nfuncTMP->maxd;
		nfuncTMP->starthistogram_size	= nfuncTMP->histogram_size;
		nfuncTMP->blurval				= 2.0;


//		nfuncTMP->startreg				= 0.05;
//		nfuncTMP->debugg_print			= debugg;
//		nfuncTMP->stdval2				= 1;
//		nfuncTMP->maxnoise				= 1;
//		nfuncTMP->ggd					= true;
//		nfuncTMP->reset();
//		nfunc->computeModel(nvec);


		dfunc_->debugg_print = false;
		nfunc_->debugg_print = false;

}

RegistrationRandom2::~RegistrationRandom2(){
	delete dfunc_;
	delete nfunc_;

	src_nr_arraypoints = 0;
	if(src_arraypoints != 0){delete src_arraypoints; src_arraypoints = 0;}
	if(src_trees3d != 0){delete src_trees3d; src_trees3d = 0;}
	if(src_a3d != 0){delete src_a3d; src_a3d = 0;}

	dst_nr_arraypoints = 0;
	if(dst_arraypoints != 0){delete dst_arraypoints; dst_arraypoints = 0;}
	if(dst_trees3d != 0){delete dst_trees3d; dst_trees3d = 0;}
	if(dst_a3d != 0){delete dst_a3d; dst_a3d = 0;}
}

void RegistrationRandom2::setSrc(std::vector<superpoint> & src_){
    src = src_;
	src_random = getRandom(src);

	Tsrc = src_random;
	transformSuperPoints(Tsrc,initTransform);
	deMean(Tsrc);

	rotationweight = getMeanDist(Tsrc);

	src_nr_arraypoints = std::min(int(Tsrc.size()),int(nr_src_target));
	src_arraypoints = new double[3*src_nr_arraypoints];

	for(unsigned int i = 0; i < src_nr_arraypoints; i++){
		const superpoint & p	= Tsrc[i];
		src_arraypoints[3*i+0]	= p.x;
		src_arraypoints[3*i+1]	= p.y;
		src_arraypoints[3*i+2]	= p.z;
	}

	src_a3d = new ArrayData3D<double>;
	src_a3d->data	= src_arraypoints;
	src_a3d->rows	= src_nr_arraypoints;
	src_trees3d	= new Tree3d(3, *src_a3d, nanoflann::KDTreeSingleIndexAdaptorParams(10));
	src_trees3d->buildIndex();
}

void RegistrationRandom2::setDst(std::vector<superpoint> & dst_){
    dst = dst_;
	dst_random = getRandom(dst);

	Tdst = dst_random;
	transformSuperPoints(Tdst,initTransform);
	deMean(Tdst);

	dst_nr_arraypoints = std::min(int(Tdst.size()),int(nr_dst_target));
	dst_arraypoints = new double[3*dst_nr_arraypoints];

	for(unsigned int i = 0; i < dst_nr_arraypoints; i++){
		const superpoint & p	= Tdst[i];
		dst_arraypoints[3*i+0]	= p.x;
		dst_arraypoints[3*i+1]	= p.y;
		dst_arraypoints[3*i+2]	= p.z;
	}

	dst_a3d = new ArrayData3D<double>;
	dst_a3d->data	= dst_arraypoints;
	dst_a3d->rows	= dst_nr_arraypoints;
	dst_trees3d	= new Tree3d(3, *dst_a3d, nanoflann::KDTreeSingleIndexAdaptorParams(10));
	dst_trees3d->buildIndex();
}

double getTime(){
	struct timeval start1;
	gettimeofday(&start1, NULL);
	return double(start1.tv_sec+(start1.tv_usec/1000000.0));
}

double transformationdiff(Eigen::Matrix4d A, Eigen::Matrix4d B, double rotationweight){
	Eigen::Matrix4d C = A.inverse()*B;
	double r = fabs(1-C(0,0))+fabs(C(0,1))+fabs(C(0,2))  +  fabs(C(1,0))+fabs(1-C(1,1))+fabs(C(1,2))  +  fabs(C(2,0))+fabs(C(2,1))+fabs(1-C(2,2));
	double t = sqrt(C(0,3)*C(0,3)+C(1,3)*C(1,3)+C(2,3)*C(2,3));
	return r*rotationweight+t;
}

bool compareFusionResults (FusionResults i,FusionResults j) { return i.score > j.score; }

std::vector<superpoint> RegistrationRandom2::getRandom(std::vector<superpoint> points){
	unsigned int nr = points.size();
	std::vector<unsigned int> inds;
	inds.resize(nr);
	for(unsigned int i = 0; i < nr; i++){
		inds[i] = i;
	}

	for(unsigned int i = 0; i < nr; i++){
		unsigned int rind = rand()%nr;
		unsigned int tmp = inds[i];
		inds[i] = inds[rind];
		inds[rind] = tmp;
	}

	std::vector<superpoint> ret;
	ret.resize(nr);
	for(unsigned int i = 0; i < nr; i++){
		ret[i] = points[inds[i]];
	}
	return ret;
}

double RegistrationRandom2::getMeanDist(std::vector<superpoint> & data){
	unsigned int nr_data = data.size();

	double sum = 0;
	for(unsigned int i = 0; i < nr_data; i++){
		double x = data[i].x;
		double y = data[i].y;
		double z = data[i].z;
		sum += sqrt(x*x+y*y+z*z);
	}
	return sum/double(nr_data);
}

void RegistrationRandom2::deMean(std::vector<superpoint> & data){
	unsigned int nr_data = data.size();

	double mean_x = 0;
	double mean_y = 0;
	double mean_z = 0;

	for(unsigned int i = 0; i < nr_data; i++){
		mean_x += data[i].x;
		mean_y += data[i].y;
		mean_z += data[i].z;
	}
	mean_x /= double(nr_data);
	mean_y /= double(nr_data);
	mean_z /= double(nr_data);

	for(unsigned int i = 0; i < nr_data; i++){
		data[i].x -= mean_x;
		//data[i].y -= mean_y;
		data[i].z -= mean_z;
	}
}

void pruneDuplicates(std::vector<FusionResults> & fr_X, double threshold, double rotationWeight){
	unsigned int nr = fr_X.size();
	std::vector<bool> to_delete;
	to_delete.resize(nr);
	for(unsigned int i = 0; i < nr; i++){
		to_delete[i] = false;
	}

	for(unsigned int i = 0; i < nr; i++){
		if(!to_delete[i]){
			for(unsigned int j = i+1; j < nr; j++){
				if(!to_delete[j]){
					to_delete[j] = transformationdiff(fr_X[i].guess,fr_X[j].guess,rotationWeight) < threshold;
				}
			}
		}
	}

	std::vector<FusionResults> fr_X2;
	for(unsigned int i = 0; i < nr; i++){
		if(!to_delete[i]){
			fr_X2.push_back(fr_X[i]);
		}
	}
	fr_X = fr_X2;
}

void RegistrationRandom2::transformSuperPoints(std::vector<superpoint> & spvec, Eigen::Matrix4d cp){
	float m00 = cp(0,0); float m01 = cp(0,1); float m02 = cp(0,2); float m03 = cp(0,3);
	float m10 = cp(1,0); float m11 = cp(1,1); float m12 = cp(1,2); float m13 = cp(1,3);
	float m20 = cp(2,0); float m21 = cp(2,1); float m22 = cp(2,2); float m23 = cp(2,3);

	unsigned int nrp = spvec.size();
	for(unsigned long i = 0; i < nrp; i++){
		superpoint & p = spvec[i];

		float x		= p.x;
		float y		= p.y;
		float z		= p.z;
		float nx	= p.nx;
		float ny	= p.ny;
		float nz	= p.nz;

		float tx	= m00*x + m01*y + m02*z + m03;
		float ty	= m10*x + m11*y + m12*z + m13;
		float tz	= m20*x + m21*y + m22*z + m23;

		float tnx	= m00*nx + m01*ny + m02*nz;
		float tny	= m10*nx + m11*ny + m12*nz;
		float tnz	= m20*nx + m21*ny + m22*nz;

		p.x	= tx;
		p.y	= ty;
		p.z	= tz;
		p.nx = tnx;
		p.ny = tny;
		p.nz = tnz;
	}
}

std::vector<reglib::superpoint> RegistrationRandom2::getMatches(const std::vector<reglib::superpoint> & X, const std::vector<reglib::superpoint> & Y, Tree3d * Y_tree){
	unsigned int nrp = X.size();
	std::vector<reglib::superpoint> Q;
	Q.resize(nrp);

	//#pragma omp parallel for
	for(unsigned int i=0; i < nrp; ++i) {
		std::vector<size_t>   ret_indexes(1);
		std::vector<double> out_dists_sqr(1);
		nanoflann::KNNResultSet<double> resultSet(1);
		resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
		double qp [3];
		const superpoint & p = X[i];

		qp[0] = p.x;
		qp[1] = p.y;
		qp[2] = p.z;

		Y_tree->findNeighbors(resultSet, qp, nanoflann::SearchParams(10));
		Q[i] = Y[ret_indexes[0]];
	}
	return Q;
}

std::vector<double> RegistrationRandom2::getStdDiv(std::vector<reglib::superpoint> & A, std::vector<reglib::superpoint> & Q){
	unsigned int nr = Q.size();
	std::vector<double> ret;
	ret.resize(nr);
	for(unsigned long i = 0; i < nr; i++){
		ret[i] = sqrt(1.0/A[i].point_information+1.0/Q[i].point_information);//information to covariance then compute stddiv
	}
	return ret;
}


double RegistrationRandom2::getStdDivVec(const std::vector<double> & V){
	unsigned int nr = V.size();
	double sum = 0;
	for(unsigned long i = 0; i < nr; i++){
		double v = V[i];
		sum += v*v;
	}
	return sqrt(sum/double(nr));
}


void RegistrationRandom2::getResiduals(std::vector<double> & residuals_d, std::vector<double> & residuals_a, std::vector<double> & standard_div, std::vector<reglib::superpoint> & A, std::vector<reglib::superpoint> & Q){
	unsigned int nr_matches = A.size();
	residuals_d.resize(nr_matches);
	residuals_a.resize(nr_matches);

	for(unsigned long i = 0; i < nr_matches; i++){
		const reglib::superpoint & a = A[i];
		const reglib::superpoint & q = Q[i];

		double dx = a.x-q.x;
		double dy = a.y-q.y;
		double dz = a.z-q.z;
		float qnx = q.nx;
		float qny = q.ny;
		float qnz = q.nz;
		float di = qnx*dx + qny*dy + qnz*dz;
		residuals_d[i] = di/standard_div[i];
		residuals_a[i] = 1.0-(a.nx*qnx + a.ny*qny + a.nz*qnz);
	}
}

RegistrationOptimization RegistrationRandom2::getOptimization(std::vector<double> & standard_div, DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, std::vector<reglib::superpoint> & X, std::vector<reglib::superpoint> & Q){
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X_cld;
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y_cld;
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X_cld2;
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y_cld2;
	if(visualizationLvl == 5){
		X_cld = getPointCloudFromVector(X , 3, 0,   255,	0);
		Y_cld = getPointCloudFromVector(Q , 3, 255, 0,		0);
		X_cld2 = getPointCloudFromVector(X ,3, 0,   255,    0);
		Y_cld2 = getPointCloudFromVector(Q ,3, 0,   0,		255);
	}

    Eigen::Matrix<double, 6, 1> ATb = Eigen::Matrix<double, 6, 1>::Zero();
    Eigen::Matrix<double, 6, 6> ATA = Eigen::Matrix<double, 6, 6>::Identity();

	double info = pow(dfunc->getNoise()*nfunc->getNoise(),-2);

    unsigned int nr_matches = X.size();
    for(unsigned int i = 0; i < nr_matches; i++){
        const reglib::superpoint & x = X[i];
        const reglib::superpoint & q = Q[i];
        const double & stddiv = standard_div[i];

        const double & sx = x.x;
        const double & sy = x.y;
        const double & sz = x.z;
        const double & nx = x.nx;
        const double & ny = x.ny;
        const double & nz = x.nz;

        const double & dx = q.x;
        const double & dy = q.y;
        const double & dz = q.z;
        const double & dnx = q.nx;
        const double & dny = q.ny;
        const double & dnz = q.nz;

        double diffX = sx-dx;
        double diffY = sy-dy;
        double diffZ = sz-dz;

        const double & di    = nx*diffX + ny*diffY + nz*diffZ;
        const double & angle = nx*dnx+ny*dny+nz*dnz;

		double prob = dfunc->getProb(di/stddiv)*nfunc->getProb(1.0-angle);
		double weight = info*prob/(stddiv*stddiv);

		if(visualizationLvl == 5){
			X_cld2->points[i].r = 255*prob;
			X_cld2->points[i].g = 255*prob;
			X_cld2->points[i].b = 255*prob;
		}

        const double & a = nz*sy - ny*sz;
        const double & b = nx*sz - nz*sx;
        const double & c = ny*sx - nx*sy;

        ATA.coeffRef ( 0) += weight * a  * a;
        ATA.coeffRef ( 1) += weight * a  * b;
        ATA.coeffRef ( 2) += weight * a  * c;
        ATA.coeffRef ( 3) += weight * a  * nx;
        ATA.coeffRef ( 4) += weight * a  * ny;
        ATA.coeffRef ( 5) += weight * a  * nz;
        ATA.coeffRef ( 7) += weight * b  * b;
        ATA.coeffRef ( 8) += weight * b  * c;
        ATA.coeffRef ( 9) += weight * b  * nx;
        ATA.coeffRef (10) += weight * b  * ny;
        ATA.coeffRef (11) += weight * b  * nz;
        ATA.coeffRef (14) += weight * c  * c;
        ATA.coeffRef (15) += weight * c  * nx;
        ATA.coeffRef (16) += weight * c  * ny;
        ATA.coeffRef (17) += weight * c  * nz;
        ATA.coeffRef (21) += weight * nx * nx;
        ATA.coeffRef (22) += weight * nx * ny;
        ATA.coeffRef (23) += weight * nx * nz;
        ATA.coeffRef (28) += weight * ny * ny;
        ATA.coeffRef (29) += weight * ny * nz;
        ATA.coeffRef (35) += weight * nz * nz;

        const double & d = weight * (nx*dx + ny*dy + nz*dz - nx*sx - ny*sy - nz*sz);

        ATb.coeffRef (0) += a * d;
        ATb.coeffRef (1) += b * d;
        ATb.coeffRef (2) += c * d;
        ATb.coeffRef (3) += nx * d;
        ATb.coeffRef (4) += ny * d;
        ATb.coeffRef (5) += nz * d;
    }

	if(visualizationLvl == 5){
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (X_cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X_cld), "scloud1",1);
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y_cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y_cld), "dcloud1",1);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scloud1");
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "dcloud1");

		viewer->addPointCloud<pcl::PointXYZRGBNormal> (X_cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X_cld2), "scloud2",2);
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y_cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y_cld2), "dcloud2",2);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scloud2");
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "dcloud2");

		viewer->spin();
		viewer->removeAllPointClouds();
	}

    return RegistrationOptimization (ATb,ATA);
}

void RegistrationRandom2::show(std::vector<reglib::superpoint> & X1, std::vector<reglib::superpoint> & Y1,std::vector<reglib::superpoint> & X2, std::vector<reglib::superpoint> & Y2 ){
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X_cld = getPointCloudFromVector(X1 ,3, 0,   255,0);
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y_cld = getPointCloudFromVector(Y1      ,3, 255, 0,  0);
	viewer->addPointCloud<pcl::PointXYZRGBNormal> (X_cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X_cld), "scloud1",1);
	viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y_cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y_cld), "dcloud1",1);
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scloud1");
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "dcloud1");
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X_cld2 = getPointCloudFromVector(X2 ,3, 0,   255,    0);
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y_cld2 = getPointCloudFromVector(Y2 ,3, 0,     0,  255);
	viewer->addPointCloud<pcl::PointXYZRGBNormal> (X_cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X_cld2), "scloud2",2);
	viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y_cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y_cld2), "dcloud2",2);
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scloud2");
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "dcloud2");
	viewer->spin();
	viewer->removeAllPointClouds();
}


FusionResults RegistrationRandom2::refine(FusionResults guess, unsigned int nrp){
	unsigned int src_nrp = std::min(int(src.size()),int(nrp));
	unsigned int dst_nrp = std::min(int(dst.size()),int(nrp));

	Eigen::Matrix4d cp = guess.guess;

	std::vector<superpoint> X = Tsrc;
	std::vector<superpoint> Y = Tdst;
	X.resize(src_nrp);
	transformSuperPoints(X,cp);

    std::vector<superpoint> Xstart = X;

	double score = 0;




	//copy  functions
	DistanceWeightFunction2 * dfunc = dfunc_->clone();
	DistanceWeightFunction2 * nfunc = nfunc_->clone();

	//Compute first pair of matches
	std::vector<reglib::superpoint> Q = getMatches(X,Y,dst_trees3d);


	//Compute regularization
	std::vector<double> residuals_d;
	std::vector<double> residuals_a;
	std::vector<double> residual_standard_div = getStdDiv(X,Q);

	getResiduals(residuals_d,residuals_a,residual_standard_div,X,Q);
	dfunc->setRegularization(getStdDivVec(residuals_d));
	nfunc->setRegularization(getStdDivVec(residuals_a));

	//Reset functions
	dfunc->reset();
	nfunc->reset();

//	dfunc->debugg_print = true;
//	nfunc->debugg_print = true;

	double start = getTime();

	bool timestopped = false;
	double maxtime = 0.1;

	//printf("scores = [");

	//iterate over regularization
	for(unsigned int funcupdate=0; funcupdate < 100; ++funcupdate) {
		if( (getTime()-start) > maxtime ){timestopped = true; break;}
		if(visualizationLvl == 2){show(Xstart,Y,X,Y );}

		//Normal ICP loop
		for(unsigned int rematching=0; rematching < 20; ++rematching) {
			Eigen::Matrix4d before_rematch = cp;
			if(rematching > 0){
                Q = getMatches(X,Y,dst_trees3d);
				residual_standard_div = getStdDiv(X,Q);
			}
            //Retrain
			if(rematching % 2 == 0){
				getResiduals(residuals_d,residuals_a,residual_standard_div,X,Q);
				dfunc->computeModel(residuals_d);
				nfunc->computeModel(residuals_a);
			}

			if(visualizationLvl == 3){show(Xstart,Y,X,Y );}

            //Optimize
            for(unsigned int optimization=0; optimization < 10; ++optimization) {
				if(visualizationLvl == 4){show(Xstart,Y,X,Y );}


                RegistrationOptimization ro = getOptimization(residual_standard_div,dfunc,nfunc,X,Q);
				Eigen::Matrix4d To = ro.solve().matrix().inverse();

				score = ro.getScore();

                transformSuperPoints(X,To);
				cp = To*cp;

				if(transformationdiff(To,Eigen::Matrix4d::Identity(),rotationweight) < dfunc->getNoise()*0.01){//Check convergence
					break;
				}
            }

			if(transformationdiff(cp,before_rematch,rotationweight) < dfunc->getNoise()*0.1){//Check convergence
				break;
			}
		}

		dfunc->update();
		nfunc->update();
		if(dfunc->regularization/dfunc->getNoise() < 0.01 && nfunc->regularization/nfunc->getNoise() < 0.01){
			break;
		}
	}

	//printf("];\n");

	if(visualizationLvl > 0){show(Xstart,Y,X,Y );}


	guess.guess = cp;
	guess.score = score;

	return guess;
}

FusionResults RegistrationRandom2::getTransform(Eigen::MatrixXd guess){
	double startTime = getTime();
	std::vector< double > rxs;
	std::vector< double > rys;
	std::vector< double > rzs;

	std::vector< double > txs;
	std::vector< double > tys;
	std::vector< double > tzs;



	std::vector<FusionResults> fr_X;

	for(double rx = 0; rx < steprx; rx++){
		for(double ry = 0; ry < stepry; ry++){
			for(double rz = 0; rz < steprz; rz++){
				for(double tx = 0; tx < steptx; tx++){
					for(double ty = 0; ty < stepty; ty++){
						for(double tz = 0; tz < steptz; tz++){

							Eigen::Affine3d randomrot = Eigen::Affine3d::Identity();
							randomrot =	Eigen::AngleAxisd(start_rx + rx*(stop_rx-start_rx)/double(steprx+1), Eigen::Vector3d::UnitX()) *
										Eigen::AngleAxisd(start_ry + ry*(stop_ry-start_ry)/double(stepry+1), Eigen::Vector3d::UnitY()) *
										Eigen::AngleAxisd(start_rz + rz*(stop_rz-start_rz)/double(steprz+1), Eigen::Vector3d::UnitZ());

							//Eigen::Affine3d guess =  Ymean*randomrot*Xmean.inverse();
							Eigen::Affine3d guess =  randomrot;

							fr_X.push_back(FusionResults(guess.matrix(),0));

						}
					}
				}
			}
		}
	}

	if(visualizationLvl > 0){
		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X_cld = getPointCloudFromVector(src,3,0,255,0);
		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y_cld = getPointCloudFromVector(dst,3,255,0,0);
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (X_cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X_cld), "scloud1",1);
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y_cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y_cld), "dcloud1",1);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scloud1");
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "dcloud1");

		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X_cld2 = getPointCloudFromVector(Tsrc,3,0,255,0);
		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y_cld2 = getPointCloudFromVector(Tdst,3,255,0,0);
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (X_cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X_cld2), "scloud2",2);
		viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y_cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y_cld2), "dcloud2",2);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scloud2");
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "dcloud2");

		viewer->spin();
		viewer->removeAllPointClouds();
	}

	printf("%s::%5.5fs\n",__PRETTY_FUNCTION__,getTime()-startTime);

	double rejection_area = 0.1;
	unsigned int target_nr = fr_X.size();
	for(unsigned int tp = 500; tp <= 16000; tp *= 2){
		printf("------------ %i ------------\n",tp);
		unsigned int nr_X = fr_X.size();
		unsigned int todo = std::min(nr_X,target_nr);
		if(visualizationLvl == 0){
			#pragma omp parallel for num_threads(8) schedule(dynamic)
			for(unsigned int ax = 0; ax < todo; ax++){//Register
				fr_X[ax] = refine(fr_X[ax],tp);
			}
		}else{
			for(unsigned int ax = 0; ax < todo; ax++){//Register
				double beforescore = fr_X[ax].score;
				fr_X[ax] = refine(fr_X[ax],tp);
				printf("%i -> %i/%i before: %f after: %f\n",tp,ax+1,todo,beforescore,fr_X[ax].score);
			}
		}

		//Sort scores
		std::sort (fr_X.begin(), fr_X.end(), compareFusionResults);


		//Prune half
		target_nr = std::max(8,int(target_nr/2));

//		while(fr_X.size() > target_nr){
//			//printf("fr_X.size(): %i\n",fr_X.size());
//			//Prune duplicates
//			pruneDuplicates(fr_X, rejection_area,rotationweight);
//			if(fr_X.size() > target_nr){
//				rejection_area *= 1.1;
//			}
//		}

		pruneDuplicates(fr_X, rejection_area,rotationweight);

		//fr_X.resize(std::min(int(target_nr),int(fr_X.size())));

	}

	FusionResults fr;
	for(unsigned int ax = 0; ax < fr_X.size() && ax < 500; ax++){
		fr.candidates.push_back(fr_X[ax].guess);
		fr.counts.push_back(1);
		fr.scores.push_back(fr_X[ax].score);
	}

	return fr;
}

}
