#include "registration/MassRegistrationPPR3.h"

#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace reglib
{


//Surface
void DataNode::init(){
    surface_nrp = 0;
    surface_p = 0;
    surface_n = 0;
    surface_i = 0;
    surface_valid = 0;
    surface_trees3d = 0;
    surface_a3d = 0;
    surface_active = 0;
    surface_nr_active = 0;
    meandist = 1;

    //Randomized color
    randr = 56+rand()%200;
    randg = 56+rand()%200;
    randb = 56+rand()%200;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr DataNode::getPCLcloud(Eigen::Matrix4d p, int r, int g, int b){
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ptr (new pcl::PointCloud<pcl::PointXYZRGB>);

    const double & m00 = p(0,0); const double & m01 = p(0,1); const double & m02 = p(0,2); const double & m03 = p(0,3);
    const double & m10 = p(1,0); const double & m11 = p(1,1); const double & m12 = p(1,2); const double & m13 = p(1,3);
    const double & m20 = p(2,0); const double & m21 = p(2,1); const double & m22 = p(2,2); const double & m23 = p(2,3);

    for(unsigned int i = 0; i < surface_nrp; i++){
        double x  = surface_p[3*i+0];
        double y  = surface_p[3*i+1];
        double z  = surface_p[3*i+2];

        pcl::PointXYZRGB p;
        p.x = m00*x + m01*y + m02*z + m03;
        p.y = m10*x + m11*y + m12*z + m13;
        p.z = m20*x + m21*y + m22*z + m23;
        p.r = r;
        p.g = g;
        p.b = b;

        cloud_ptr->points.push_back(p);
    }
    return cloud_ptr;
}

DataNode::DataNode(Model * model, bool useSurfacePoints, unsigned int surface_nr_active_target){
    init();

    double meandist_sum = 0;
    double meandist_weight = 0;
    if(useSurfacePoints){
        std::vector<superpoint> & points = model->points;

        surface_nrp = points.size();

        if(surface_nrp > 0){
            surface_p       = new double          [3*surface_nrp];
            surface_n       = new double          [3*surface_nrp];
            surface_i       = new double          [  surface_nrp];
            surface_valid   = new bool            [  surface_nrp];
            surface_active  = new unsigned int    [  surface_nrp];
            for(unsigned int i = 0; i < surface_nrp; i++){
                superpoint & p      = points[i];
                surface_p[3*i+0]    = p.x;
                surface_p[3*i+1]    = p.y;
                surface_p[3*i+2]    = p.z;
                surface_n[3*i+0]    = p.nx;
                surface_n[3*i+1]    = p.ny;
                surface_n[3*i+2]    = p.nz;
                surface_i[i]        = p.point_information;
                surface_valid[i]    = !p.is_boundry;
                surface_active[i]   = i;

                meandist_sum    += p.point_information*sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
                meandist_weight += p.point_information;
            }

            //srand(0);
            for(unsigned int i = 0; i < surface_nrp; i++){
                int rind             = rand()%surface_nrp;
                int tmp              = surface_active[i];
                surface_active[i]    = surface_active[rind];
                surface_active[rind] = tmp;
            }
            surface_nr_active   = std::min(surface_nr_active_target,surface_nrp);

            surface_a3d         = new ArrayData3D<double>;
            surface_a3d->data	= surface_p;
            surface_a3d->rows	= surface_nrp;
            surface_trees3d     = new Tree3d(3, *surface_a3d, nanoflann::KDTreeSingleIndexAdaptorParams(10));
            surface_trees3d->buildIndex();
        }
    }
    meandist = meandist_sum/meandist_weight;
}

DataNode::~DataNode(){
    if(surface_nrp != 0){
        delete[] surface_p;
        delete[] surface_n;
        delete[] surface_i;
        delete[] surface_valid;
        delete[] surface_active;
        delete   surface_a3d;
        delete   surface_trees3d;
    }
}


void EdgeData::showMatches(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, Eigen::Matrix4d p){
    char buf [1024];
    viewer->removeAllShapes();
    viewer->removeAllPointClouds();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud1 = node1->getPCLcloud(p                          , 0, 255,   0);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2 = node2->getPCLcloud(Eigen::Matrix4d::Identity(), 0,   0, 255);
    for(unsigned int i = 0; i < surface_match1.size() && i < 5000; i++){
        sprintf(buf,"line%i",i); viewer->addLine<pcl::PointXYZRGB> (cloud1->points[surface_match1[i]],cloud2->points[surface_match2[i]],1,0,0,buf);
    }
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cloud1), "scloud");
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cloud2), "dcloud");
    viewer->spin();
}

void EdgeData::showMatchesW(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, DistanceWeightFunction2 * func, Eigen::Matrix4d p1, Eigen::Matrix4d p2){
    char buf [1024];
    viewer->removeAllShapes();
    viewer->removeAllPointClouds();

    const double & m001 = p1(0,0); const double & m011 = p1(0,1); const double & m021 = p1(0,2); const double & m031 = p1(0,3);
    const double & m101 = p1(1,0); const double & m111 = p1(1,1); const double & m121 = p1(1,2); const double & m131 = p1(1,3);
    const double & m201 = p1(2,0); const double & m211 = p1(2,1); const double & m221 = p1(2,2); const double & m231 = p1(2,3);

    const double & m002 = p2(0,0); const double & m012 = p2(0,1); const double & m022 = p2(0,2); const double & m032 = p2(0,3);
    const double & m102 = p2(1,0); const double & m112 = p2(1,1); const double & m122 = p2(1,2); const double & m132 = p2(1,3);
    const double & m202 = p2(2,0); const double & m212 = p2(2,1); const double & m222 = p2(2,2); const double & m232 = p2(2,3);

    double * surface_p1             = node1->surface_p;
    double * surface_n1             = node1->surface_n;

    double * surface_p2             = node2->surface_p;
    double * surface_n2             = node2->surface_n;

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud1 = node1->getPCLcloud(p1, 0, 255,   0);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2 = node2->getPCLcloud(p2, 0,   0, 255);


    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ptr (new pcl::PointCloud<pcl::PointXYZRGB>);

    double stopD2 = pow(10.0*func->getNoise(),2);

    unsigned int nr_matches = surface_match1.size();
    for(unsigned int i = 0; i < nr_matches; i++){
        unsigned int i1 = surface_match1[i];
        unsigned int i2 = surface_match2[i];
        double rw       = surface_rangew[i];

        double src_x  = surface_p1[3*i1+0];
        double src_y  = surface_p1[3*i1+1];
        double src_z  = surface_p1[3*i1+2];
        double src_nx = surface_n1[3*i1+0];
        double src_ny = surface_n1[3*i1+1];
        double src_nz = surface_n1[3*i1+2];

        const double & sx = m001*src_x + m011*src_y + m021*src_z + m031;
        const double & sy = m101*src_x + m111*src_y + m121*src_z + m131;
        const double & sz = m201*src_x + m211*src_y + m221*src_z + m231;
        const double & nx = m001*src_nx + m011*src_ny + m021*src_nz;
        const double & ny = m101*src_nx + m111*src_ny + m121*src_nz;
        const double & nz = m201*src_nx + m211*src_ny + m221*src_nz;

        double dst_x  = surface_p2[3*i2+0];
        double dst_y  = surface_p2[3*i2+1];
        double dst_z  = surface_p2[3*i2+2];
        double dst_nx = surface_n2[3*i2+0];
        double dst_ny = surface_n2[3*i2+1];
        double dst_nz = surface_n2[3*i2+2];

        const double & dx = m002*dst_x + m012*dst_y + m022*dst_z + m032;
        const double & dy = m102*dst_x + m112*dst_y + m122*dst_z + m132;
        const double & dz = m202*dst_x + m212*dst_y + m222*dst_z + m232;
        const double & dnx = m002*dst_nx + m012*dst_ny + m022*dst_nz;
        const double & dny = m102*dst_nx + m112*dst_ny + m122*dst_nz;
        const double & dnz = m202*dst_nx + m212*dst_ny + m222*dst_nz;

        const double & angle = nx*dnx+ny*dny+nz*dnz;

//        sprintf(buf,"line%i",i);
//        if(angle < 0){
//			viewer->addLine<pcl::PointXYZRGB> (cloud1->points[surface_match1[i]],cloud2->points[surface_match2[i]],1,0,0,buf);
//            continue;
//        }

        double diffX = sx-dx;
        double diffY = sy-dy;
        double diffZ = sz-dz;

//        double d2 = diffX*diffX + diffY*diffY + diffZ*diffZ;
//        if(d2 > stopD2){
//			viewer->addLine<pcl::PointXYZRGB> (cloud1->points[surface_match1[i]],cloud2->points[surface_match2[i]],1,0,0,buf);
//            continue;
//        }

        double di = (nx*diffX + ny*diffY + nz*diffZ)*rw;

        double prob = func->getProb(di);

        cloud_ptr->points.push_back(cloud1->points[surface_match1[i]]);
        cloud_ptr->points.back().r = 255.0*(1-prob);
        cloud_ptr->points.back().g = 255.0*prob;
        cloud_ptr->points.back().b = 0;
    }

    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_ptr, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cloud_ptr), "scloud");
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "scloud");
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cloud2), "dcloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "dcloud");
    viewer->spin();
}

EdgeData::EdgeData(DataNode * node1_, DataNode * node2_){
    score = 0;
    node1 = node1_;
    node2 = node2_;
    surface_match1.resize(0);
    surface_match2.resize(0);
    surface_rangew.resize(0);
    rematchPoseLast = Eigen::Matrix4d::Identity();
    optPoseLast = Eigen::Matrix4d::Identity();
    modelPoseLast = Eigen::Matrix4d::Identity();
    rematched = false;
    nr_matches_orig = 0;
}
EdgeData::~EdgeData(){}

bool EdgeData::needRefinement(Eigen::Matrix4d p, double convergence){
    if(rematched){return true;}
    double change = getChange(optPoseLast,p,node1->meandist + node2->meandist);
    //printf("change: %f convergence: %f\n",change,convergence);
    if(change < convergence){return false;}
    return true;
}

int EdgeData::surfaceRematch(Eigen::Matrix4d p,double convergence, bool force){
    if(!force && getChange(rematchPoseLast,p,node1->meandist + node2->meandist) < convergence){return -1;}
    rematchPoseLast = p;
    rematched = true;

    double * surface_p1             = node1->surface_p;
    double * surface_i1             = node1->surface_i;
    Tree3d * surface_trees3d1       = node1->surface_trees3d;
    unsigned int surface_nr_active1 = node1->surface_nr_active;
    unsigned int * surface_active1  = node1->surface_active;

    double * surface_p2             = node2->surface_p;
    double * surface_i2             = node2->surface_i;
    bool * surface_v2               = node2->surface_valid;
    Tree3d * surface_trees3d2       = node2->surface_trees3d;
    unsigned int surface_nrp2       = node2->surface_nrp;

    const double & m00 = p(0,0); const double & m01 = p(0,1); const double & m02 = p(0,2); const double & m03 = p(0,3);
    const double & m10 = p(1,0); const double & m11 = p(1,1); const double & m12 = p(1,2); const double & m13 = p(1,3);
    const double & m20 = p(2,0); const double & m21 = p(2,1); const double & m22 = p(2,2); const double & m23 = p(2,3);

    size_t ret_indexes      [1];
    double out_dists_sqr    [1];
    double tpoint           [3];
    nanoflann::KNNResultSet<double> resultSet(1);

    surface_match1.resize(surface_nr_active1);
    surface_match2.resize(surface_nr_active1);
    surface_rangew.resize(surface_nr_active1);

    int added = 0;
    for(unsigned int i = 0; i < surface_nr_active1; i++){
        unsigned int ind = surface_active1[i];
        double x  = surface_p1[3*ind+0];
        double y  = surface_p1[3*ind+1];
        double z  = surface_p1[3*ind+2];

        tpoint[0] = m00*x + m01*y + m02*z + m03;
        tpoint[1] = m10*x + m11*y + m12*z + m13;
        tpoint[2] = m20*x + m21*y + m22*z + m23;

        resultSet.init(ret_indexes, out_dists_sqr);
        surface_trees3d2->findNeighbors(resultSet,tpoint, nanoflann::SearchParams(10));
        int match_id = ret_indexes[0];
        if(match_id >= 0 && match_id < surface_nrp2){
            if(surface_v2[match_id]){
                surface_match1[added] = ind;
                surface_match2[added] = match_id;
                surface_rangew[added] = sqrt( 1.0/(1.0/surface_i1[ind]+1.0/surface_i2[match_id]));
                added++;
            }
        }
    }
    surface_match1.resize(added);
    surface_match2.resize(added);
    surface_rangew.resize(added);
    nr_matches_orig = added;
    return 0;
}

void EdgeData::computeSurfaceResiduals(Eigen::Matrix4d p, double * residuals, unsigned long & counter){
    double * surface_p1             = node1->surface_p;
    double * surface_n1             = node1->surface_n;

    double * surface_p2             = node2->surface_p;
    double * surface_n2             = node2->surface_n;

    const double & m00 = p(0,0); const double & m01 = p(0,1); const double & m02 = p(0,2); const double & m03 = p(0,3);
    const double & m10 = p(1,0); const double & m11 = p(1,1); const double & m12 = p(1,2); const double & m13 = p(1,3);
    const double & m20 = p(2,0); const double & m21 = p(2,1); const double & m22 = p(2,2); const double & m23 = p(2,3);

    unsigned int nr_matches = surface_match1.size();
    for(unsigned int i = 0; i < nr_matches; i++){
        unsigned int i1 = surface_match1[i];
        unsigned int i2 = surface_match2[i];
        double rw       = surface_rangew[i];

		double x  = surface_p1[3*i1+0];
		double y  = surface_p1[3*i1+1];
		double z  = surface_p1[3*i1+2];
		double tx = m00*x + m01*y + m02*z + m03;
		double ty = m10*x + m11*y + m12*z + m13;
		double tz = m20*x + m21*y + m22*z + m23;

//		double nx  = surface_n1[3*i1+0];
//		double ny  = surface_n1[3*i1+1];
//		double nz  = surface_n1[3*i1+2];
//		double tnx = m00*nx + m01*ny + m02*nz;
//		double tny = m10*nx + m11*ny + m12*nz;
//		double tnz = m20*nx + m21*ny + m22*nz;

        float dx = tx-surface_p2[3*i2+0];
        float dy = ty-surface_p2[3*i2+1];
        float dz = tz-surface_p2[3*i2+2];

        float qx = surface_n2[3*i2+0];
        float qy = surface_n2[3*i2+1];
        float qz = surface_n2[3*i2+2];

        float di = qx*dx + qy*dy + qz*dz;

//		const double & angle = qx*tnx+qy*tny+qz*tnz;
//		if(angle < 0){ continue;}

//		double d2 = dx*dx + dy*dy + dz*dz;
//		if(d2 > stopD2){continue;}

        residuals[counter] = di*rw;
        counter++;
    }
}

void EdgeData::addSurfaceOptimization(DistanceWeightFunction2 * func, Eigen::Matrix4d p1, Eigen::Matrix4d p2, Matrix6d & ATA, Vector6d & ATb ){
    rematched = false;
    optPoseLast = p2.inverse()*p1;

    double info = pow(func->getNoise(),-2);
    double stopD2 = pow(10.0*func->getNoise(),2);

    const double & m001 = p1(0,0); const double & m011 = p1(0,1); const double & m021 = p1(0,2); const double & m031 = p1(0,3);
    const double & m101 = p1(1,0); const double & m111 = p1(1,1); const double & m121 = p1(1,2); const double & m131 = p1(1,3);
    const double & m201 = p1(2,0); const double & m211 = p1(2,1); const double & m221 = p1(2,2); const double & m231 = p1(2,3);

    const double & m002 = p2(0,0); const double & m012 = p2(0,1); const double & m022 = p2(0,2); const double & m032 = p2(0,3);
    const double & m102 = p2(1,0); const double & m112 = p2(1,1); const double & m122 = p2(1,2); const double & m132 = p2(1,3);
    const double & m202 = p2(2,0); const double & m212 = p2(2,1); const double & m222 = p2(2,2); const double & m232 = p2(2,3);

    double * surface_p1             = node1->surface_p;
    double * surface_n1             = node1->surface_n;

    double * surface_p2             = node2->surface_p;
    double * surface_n2             = node2->surface_n;

    score = 0;

    unsigned int nr_matches = surface_match1.size();
    for(unsigned int i = 0; i < nr_matches; i++){
        unsigned int i1 = surface_match1[i];
        unsigned int i2 = surface_match2[i];
        double rw       = surface_rangew[i];

        double src_x  = surface_p1[3*i1+0];
        double src_y  = surface_p1[3*i1+1];
        double src_z  = surface_p1[3*i1+2];
        double src_nx = surface_n1[3*i1+0];
        double src_ny = surface_n1[3*i1+1];
        double src_nz = surface_n1[3*i1+2];

        const double & sx = m001*src_x + m011*src_y + m021*src_z + m031;
        const double & sy = m101*src_x + m111*src_y + m121*src_z + m131;
        const double & sz = m201*src_x + m211*src_y + m221*src_z + m231;
        const double & nx = m001*src_nx + m011*src_ny + m021*src_nz;
        const double & ny = m101*src_nx + m111*src_ny + m121*src_nz;
        const double & nz = m201*src_nx + m211*src_ny + m221*src_nz;

        double dst_x  = surface_p2[3*i2+0];
        double dst_y  = surface_p2[3*i2+1];
        double dst_z  = surface_p2[3*i2+2];
        double dst_nx = surface_n2[3*i2+0];
        double dst_ny = surface_n2[3*i2+1];
        double dst_nz = surface_n2[3*i2+2];

        const double & dx = m002*dst_x + m012*dst_y + m022*dst_z + m032;
        const double & dy = m102*dst_x + m112*dst_y + m122*dst_z + m132;
        const double & dz = m202*dst_x + m212*dst_y + m222*dst_z + m232;
        const double & dnx = m002*dst_nx + m012*dst_ny + m022*dst_nz;
        const double & dny = m102*dst_nx + m112*dst_ny + m122*dst_nz;
        const double & dnz = m202*dst_nx + m212*dst_ny + m222*dst_nz;

        double diffX = sx-dx;
        double diffY = sy-dy;
        double diffZ = sz-dz;

        double di = (nx*diffX + ny*diffY + nz*diffZ)*rw;

        double prob = func->getProb(di);
        if(i < nr_matches_orig){
            score += prob;
        }

		const double & angle = nx*dnx+ny*dny+nz*dnz;
		if(angle < 0){continue;}
		double d2 = diffX*diffX + diffY*diffY + diffZ*diffZ;
		if(d2 > stopD2){continue;}

        double weight = info * prob*rw*rw;

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
    score /= double(node1->surface_nr_active);
}



MassRegistrationPPR3::MassRegistrationPPR3(){

    GaussianDistribution * gd = new GaussianDistribution();
    reglib::DistanceWeightFunction2PPR3 * sfunc  = new reglib::DistanceWeightFunction2PPR3(gd);
    sfunc->debugg_print                          = false;
    sfunc->noise_min                             = 0.0001;
    sfunc->startreg                              = 0.01;
    sfunc->reg_shrinkage                         = 0.5;
    surface_func = sfunc;
    useSurfacePoints = true;
    convergence_mul = 0.1;

    func_setup = 0;
}

MassRegistrationPPR3::~MassRegistrationPPR3(){
    for(unsigned int i = 0; i < nodes.size(); i++){delete nodes[i];}
    for(unsigned int i = 0; i < edges.size(); i++){
        for(unsigned int j = 0; j < edges[i].size(); j++){
            if(edges[i][j] != 0){ delete edges[i][j]; }
        }
    }

    for(unsigned int i = 0; i < edge_surface_func.size(); i++){
        for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
            if(edge_surface_func[i][j] != 0){ delete edge_surface_func[i][j]; }
        }
    }

    delete surface_func;
}



int MassRegistrationPPR3::total_nonconverged(std::vector<Eigen::Matrix4d> before, std::vector<Eigen::Matrix4d> after){
    double convergence = getConvergence();
    int not_converged = 0;
    for(unsigned int i = 0; i < edges.size(); i++){
        for(unsigned int j = 0; j < edges[i].size(); j++){
            if(edges[i][j] != 0){

                Eigen::Matrix4d dbefore = before[j].inverse()*before[i];
                Eigen::Matrix4d dafter = after[j].inverse()*after[i];

                double change = getChange(dbefore,dafter,edges[i][j]->node1->meandist + edges[i][j]->node2->meandist);

                if(func_setup == 0){
                    //printf("model %i to %i -> change %f convergence %f\n",i,j,change,convergence);
                    if(change >= convergence){
                        not_converged++;
                    }
                }else if(func_setup == 1 ){
                    //printf("model %i to %i -> change %f convergence %f\n",i,j,change,convergence_mul*edge_surface_func[i][j]->getNoise());
                    if(change >= convergence_mul*edge_surface_func[i][j]->getNoise()){
                        not_converged++;
                    }
                }
            }
        }
    }
    //printf("not_converged: %i\n",not_converged);
    return not_converged;
}

double MassRegistrationPPR3::getConvergence(){return convergence_mul*surface_func->getNoise();}

void MassRegistrationPPR3::show(std::vector<Eigen::Matrix4d> poses, bool stop){
    viewer->removeAllShapes();
    viewer->removeAllPointClouds();
    char buf [1024];
    for(unsigned int i = 0; i < nodes.size(); i++){
        sprintf(buf,"cloud%i",i);
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld = nodes[i]->getPCLcloud(poses[i], nodes[i]->randr, nodes[i]->randg, nodes[i]->randb);
        viewer->addPointCloud<pcl::PointXYZRGB> (cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld), buf);
    }
    if(stop){    viewer->spin();}
    else{        viewer->spinOnce();}
}

void MassRegistrationPPR3::addModel(Model * model, int active){
    nodes.push_back( new DataNode(model,useSurfacePoints,active) );
    edges.push_back(std::vector< EdgeData * >());
    edges.back().resize(nodes.size()-1);
    for(unsigned int i = 0; i < edges.back().size(); i++){edges.back()[i] = 0;}
    for(unsigned int i = 0; i < edges.size(); i++){edges[i].push_back(0);}

    if(func_setup == 1){
        edge_surface_func.push_back(std::vector< DistanceWeightFunction2 * >());
        edge_surface_func.back().resize(nodes.size()-1);
        for(unsigned int i = 0; i < edge_surface_func.back().size(); i++){  edge_surface_func.back()[i] = 0;}
        for(unsigned int i = 0; i < edge_surface_func.size(); i++){         edge_surface_func[i].push_back(0);}
    }
}

void MassRegistrationPPR3::removeLastNode(){
    if(nodes.size() > 0){
        delete nodes.back();
        nodes.pop_back();
    }

    if(edges.size() > 0){
        for(unsigned int i = 0; i < edges.size(); i++){
           if(edges[i].back() != 0){delete edges[i].back();}
           edges[i].pop_back();
           if(edges.back()[i] != 0){delete edges.back()[i];}
        }
        edges.pop_back();
    }

    if(edge_surface_func.size() > 0){
        for(unsigned int i = 0; i < edge_surface_func.size(); i++){
           if(edge_surface_func[i].back() != 0){delete edge_surface_func[i].back();}
           edge_surface_func[i].pop_back();
           if(edge_surface_func.back()[i] != 0){delete edge_surface_func.back()[i];}
        }
        edge_surface_func.pop_back();
    }
}

unsigned long MassRegistrationPPR3::rematch(std::vector<Eigen::Matrix4d> poses, bool force, bool debugg_print){
    double convergence = getConvergence();
    unsigned long sum_surface_rematches = 0;
    for(unsigned int i = 0; i < edges.size(); i++){
        for(unsigned int j = 0; j < edges[i].size(); j++){
            if(edges[i][j] != 0){
                Eigen::Matrix4d p = poses[j].inverse()*poses[i];
                sum_surface_rematches += edges[i][j]->surfaceRematch(p,convergence,force) != -1;
            }
        }
    }

    for(unsigned int i = 0; i < edges.size(); i++){
        for(unsigned int j = i+1; j < edges[i].size(); j++){
            if(edges[i][j] != 0){
                Eigen::Matrix4d p = poses[j].inverse()*poses[i];

                int nr_matches1 = edges[i][j]->surface_match1.size();
                int nr_matches2 = edges[j][i]->surface_match1.size();

                //edges[i][j]->showMatches(viewer,p);
                //edges[j][i]->showMatches(viewer,p.inverse());


                for(unsigned int k = 0; k < nr_matches1; k++){
                    edges[j][i]->surface_match2.push_back(edges[i][j]->surface_match1[k]);
                    edges[j][i]->surface_match1.push_back(edges[i][j]->surface_match2[k]);
                    edges[j][i]->surface_rangew.push_back(edges[i][j]->surface_rangew[k]);
                }

                for(unsigned int k = 0; k < nr_matches2; k++){
                    edges[i][j]->surface_match2.push_back(edges[j][i]->surface_match1[k]);
                    edges[i][j]->surface_match1.push_back(edges[j][i]->surface_match2[k]);
                    edges[i][j]->surface_rangew.push_back(edges[j][i]->surface_rangew[k]);
                }

                //edges[i][j]->showMatches(viewer,p);
                //edges[j][i]->showMatches(viewer,p.inverse());
            }
        }
    }


    return sum_surface_rematches > 0;
}

unsigned long MassRegistrationPPR3::model(std::vector<Eigen::Matrix4d> poses, bool force, bool debugg_print){

    if(useSurfacePoints){
        if(func_setup == 0){
            unsigned long sum_surface_rematches = 0;
            for(unsigned int i = 0; i < edges.size(); i++){
                for(unsigned int j = 0; j < edges[i].size(); j++){
                    if(edges[i][j] != 0){
                        sum_surface_rematches += edges[i][j]->surface_match1.size();
                    }
                }
            }

            double * surfaceResiduals = new double[sum_surface_rematches];
            unsigned long surfaceCounter = 0;
            for(unsigned int i = 0; i < edges.size(); i++){
                for(unsigned int j = 0; j < edges[i].size(); j++){
                    if(edges[i][j] != 0){
                        Eigen::Matrix4d p = poses[j].inverse()*poses[i];
                        edges[i][j]->computeSurfaceResiduals(p,surfaceResiduals,surfaceCounter);
                    }
                }
            }
            surface_func->computeModel(surfaceResiduals,surfaceCounter,1);
            delete[] surfaceResiduals;
        }else if(func_setup == 1){
            unsigned long max_surface_rematches = 0;
            for(unsigned int i = 0; i < edges.size(); i++){
                for(unsigned int j = 0; j < edges[i].size(); j++){
                    if(edges[i][j] != 0){
                        max_surface_rematches = std::max(max_surface_rematches,edges[i][j]->surface_match1.size());
                    }
                }
            }

            double * surfaceResiduals = new double[max_surface_rematches];
            for(unsigned int i = 0; i < edges.size(); i++){
                for(unsigned int j = 0; j < edges[i].size(); j++){
                    if(edges[i][j] != 0 && edges[i][j]->surface_match1.size() != 0){
                        unsigned long surfaceCounter = 0;
                        Eigen::Matrix4d p = poses[j].inverse()*poses[i];
                        edges[i][j]->computeSurfaceResiduals(p,surfaceResiduals,surfaceCounter);
                        edge_surface_func[i][j]->computeModel(surfaceResiduals,surfaceCounter,1);
                    }
                }
            }
            delete[] surfaceResiduals;
        }
    }

    return 1;
}


unsigned long MassRegistrationPPR3::refine(std::vector<Eigen::Matrix4d> & poses, bool force, bool debugg_print){
    double convergence = getConvergence();
    for(long outer=0; outer < 150; ++outer) {

        int not_converged = 0;

        for(unsigned int i = 0; i < edges.size(); i++){
            for(unsigned int j = 0; j < edges[i].size(); j++){
                if(edges[i][j] != 0  && edges[i][j]->surface_match1.size() > 0){
                    if(func_setup == 0){
                        if(edges[i][j]->needRefinement(poses[j].inverse()*poses[i],convergence)){not_converged++;}
                    }else if(func_setup == 1){
                        if(edges[i][j]->needRefinement(poses[j].inverse()*poses[i],convergence_mul*edge_surface_func[i][j]->getNoise())){not_converged++;}
                    }
                }
            }
        }

        //printf("outer: %i not_converged: %i \n",outer,not_converged);
        if(!(force && outer == 0) && not_converged == 0){break;}


        for(unsigned int i = 0; i < edges.size(); i++){
            int not_conv = 0;
            for(unsigned int j = 0; j < edges[i].size(); j++){
                if(edges[i][j] != 0 && edges[i][j]->surface_match1.size() > 0){
                    if(func_setup == 0){
                        if(edges[i][j]->needRefinement(poses[j].inverse()*poses[i],convergence)){not_conv++;}
                    }else if(func_setup == 1){
                        if(edges[i][j]->needRefinement(poses[j].inverse()*poses[i],convergence_mul*edge_surface_func[i][j]->getNoise())){not_conv++;}
                    }
                }
                if(edges[j][i] != 0 && edges[j][i]->surface_match1.size() > 0){
                    if(func_setup == 0){
                        if(edges[j][i]->needRefinement(poses[i].inverse()*poses[j],convergence)){not_conv++;}
                    }else if(func_setup == 1){
                        if(edges[j][i]->needRefinement(poses[i].inverse()*poses[j],convergence_mul*edge_surface_func[j][i]->getNoise())){not_conv++;}
                    }
                }
            }

            if(!(force && outer == 0) && not_conv == 0){continue;}

            Matrix6d pATA;
            Vector6d pATb;
            pATA.setZero ();
            pATb.setZero ();
            for(unsigned int j = 0; j < edges[i].size(); j++){
                EdgeData * e = edges[i][j];
                DistanceWeightFunction2 * func = 0;
                if(func_setup == 0){        func = surface_func;}
                else if(func_setup == 1){   func = edge_surface_func[i][j];}
                if(e != 0 && func != 0){
                    e->addSurfaceOptimization(func , poses[i], poses[j],pATA,pATb );
                }
            }

            pATA.coeffRef (6)  = pATA.coeff (1);
            pATA.coeffRef (12) = pATA.coeff (2);
            pATA.coeffRef (13) = pATA.coeff (8);
            pATA.coeffRef (18) = pATA.coeff (3);
            pATA.coeffRef (19) = pATA.coeff (9);
            pATA.coeffRef (20) = pATA.coeff (15);
            pATA.coeffRef (24) = pATA.coeff (4);
            pATA.coeffRef (25) = pATA.coeff (10);
            pATA.coeffRef (26) = pATA.coeff (16);
            pATA.coeffRef (27) = pATA.coeff (22);
            pATA.coeffRef (30) = pATA.coeff (5);
            pATA.coeffRef (31) = pATA.coeff (11);
            pATA.coeffRef (32) = pATA.coeff (17);
            pATA.coeffRef (33) = pATA.coeff (23);
            pATA.coeffRef (34) = pATA.coeff (29);

            for(long d = 0; d < 6; d++){pATA(d,d) += 0.000000001;}



            Matrix6d nATA;
            Vector6d nATb;
            nATA.setZero ();
            nATb.setZero ();
//            for(unsigned int j = 0; j < edges[i].size(); j++){
//                if(edges[j][i] != 0){
//                    if(func_setup == 0){
//                        edges[j][j]->addSurfaceOptimization(surface_func            , poses[j], poses[i],nATA,nATb );
//                    }else if(func_setup == 1){
//                        edges[j][j]->addSurfaceOptimization(edge_surface_func[j][i] , poses[j], poses[i],nATA,nATb );
//                    }
//                }
//            }

            for(unsigned int j = 0; j < edges[i].size(); j++){
                EdgeData * e = edges[j][i];
                DistanceWeightFunction2 * func = 0;
                if(func_setup == 0){        func = surface_func;}
                else if(func_setup == 1){   func = edge_surface_func[j][i];}
                if(e != 0 && func != 0){
                    e->addSurfaceOptimization(func , poses[j], poses[i],nATA,nATb );
                }
            }

            nATA.coeffRef (6)  = nATA.coeff (1);
            nATA.coeffRef (12) = nATA.coeff (2);
            nATA.coeffRef (13) = nATA.coeff (8);
            nATA.coeffRef (18) = nATA.coeff (3);
            nATA.coeffRef (19) = nATA.coeff (9);
            nATA.coeffRef (20) = nATA.coeff (15);
            nATA.coeffRef (24) = nATA.coeff (4);
            nATA.coeffRef (25) = nATA.coeff (10);
            nATA.coeffRef (26) = nATA.coeff (16);
            nATA.coeffRef (27) = nATA.coeff (22);
            nATA.coeffRef (30) = nATA.coeff (5);
            nATA.coeffRef (31) = nATA.coeff (11);
            nATA.coeffRef (32) = nATA.coeff (17);
            nATA.coeffRef (33) = nATA.coeff (23);
            nATA.coeffRef (34) = nATA.coeff (29);

            for(long d = 0; d < 6; d++){nATA(d,d) += 0.000000001;}

            Matrix6d ATA = pATA+nATA;
            Vector6d ATb = pATb-nATb;

//            Vector6d nx = static_cast<Vector6d> (nATA.inverse () * nATb);
//            Eigen::Affine3d ntransformation = Eigen::Affine3d(constructTransformationMatrix(nx(0,0),nx(1,0),nx(2,0),nx(3,0),nx(4,0),nx(5,0)));
//            std::cout << "nATA\n" << nATA << std::endl << std::endl;
//            std::cout << "nATb\n" << nATb << std::endl << std::endl;
//            std::cout << "nx\n" << nx << std::endl << std::endl;
//            std::cout << "ntransformation\n" << ntransformation.matrix() << std::endl << std::endl;
//            std::cout << "ntransformation.inverse()\n" << ntransformation.matrix().inverse() << std::endl << std::endl;


//            Vector6d px = static_cast<Vector6d> (pATA.inverse () * pATb);
//            Eigen::Affine3d ptransformation = Eigen::Affine3d(constructTransformationMatrix(px(0,0),px(1,0),px(2,0),px(3,0),px(4,0),px(5,0)));
//            std::cout << "pATA\n" << pATA << std::endl << std::endl;
//            std::cout << "pATb\n" << pATb << std::endl << std::endl;
//            std::cout << "px\n" << px << std::endl << std::endl;
//            std::cout << "ptransformation\n" << ptransformation.matrix() << std::endl << std::endl;
//            std::cout << "ptransformation.inverse()\n" << ptransformation.matrix().inverse() << std::endl << std::endl;

            Vector6d x = static_cast<Vector6d> (ATA.inverse () * ATb);

            if(i == 0){
                Eigen::Affine3d transformation = Eigen::Affine3d(constructTransformationMatrix(-x(0,0),-x(1,0),-x(2,0),-x(3,0),-x(4,0),-x(5,0)));
                for(unsigned int j = 1; j < poses.size(); j++){
                    poses[j] = transformation*poses[j];
                }

            }else{
                Eigen::Affine3d transformation = Eigen::Affine3d(constructTransformationMatrix(x(0,0),x(1,0),x(2,0),x(3,0),x(4,0),x(5,0)));
                poses[i] = transformation*poses[i];
            }

//            Eigen::Affine3d transformation = Eigen::Affine3d(constructTransformationMatrix(x(0,0),x(1,0),x(2,0),x(3,0),x(4,0),x(5,0)));
//            std::cout << "ATA\n" << ATA << std::endl << std::endl;
//            std::cout << "ATb\n" << ATb << std::endl << std::endl;
//            std::cout << "x\n" << x << std::endl << std::endl;
//            std::cout << "transformation\n" << transformation.matrix() << std::endl << std::endl;
//            std::cout << "transformation.inverse()\n" << transformation.matrix().inverse() << std::endl << std::endl;



//poses[i] = transformation*poses[i];


            //show(poses,true);
            if(i == 0){normalizePoses(poses);}
        }
    }
	if(visualizationLvl > 4){show(poses,false);}
    return 1;
}

void MassRegistrationPPR3::normalizePoses(std::vector<Eigen::Matrix4d> & poses){
    Eigen::Matrix4d firstinv = poses.front().inverse();
    for(long i = 0; i < poses.size(); i++){poses[i] = firstinv*poses[i];}
}

MassFusionResults MassRegistrationPPR3::getTransforms(std::vector<Eigen::Matrix4d> poses){
    //srand(0);
    //printf("\n\n\n\n");
    if(visualizationLvl > 0){printf("start MassRegistrationPPR3::getTransforms(std::vector<Eigen::Matrix4d> poses)\n");}
    if(poses.size() != nodes.size()){
        printf("MassRegistrationPPR3::getTransforms ERROR:poses.size() != nodes.size()\n");
        return MassFusionResults(poses,-1);
    }

    for(unsigned int i = 0; i < edges.size(); i++){
        for(unsigned int j = 0; j < edges[i].size(); j++){
            if(i == j){continue;}
            if(edges[i][j] == 0){
                edges[i][j] = new EdgeData(nodes[i], nodes[j]);
            }
        }
    }

    surface_func->reset();

    if(func_setup == 1){
        for(unsigned int i = 0; i < edge_surface_func.size(); i++){
            for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
                if(i == j){continue;}
                if(edge_surface_func[i][j] == 0){
                    edge_surface_func[i][j] = surface_func->clone();
                }
                edge_surface_func[i][j]->reset();
            }
        }
    }

    //show(poses,true);
    rematch(poses,  true);
    model(poses,    true);
    refine(poses,    true);

    std::vector<Eigen::Matrix4d> prev = poses;

    for(int func = 0; func < 50; func++){
		if(visualizationLvl == 2){
			printf("func: %i \n",func);
			show(poses,false);
		}
        //printf("func: %i\n",func);
        //if(func == 3){return MassFusionResults(poses,1);}
        std::vector< std::vector<Eigen::Matrix4d> > prevs;
        //printf("=========================================================================\n");
        for(int i = 0; i < 150; i++){
            //printf("func: %i i:%i \n",func,i);
            rematch(poses,  false,func == 2);
            model(poses,    false,func == 2);
            refine(poses,    true,func == 2);

            int not_converged = total_nonconverged(prev,poses);
            for(unsigned int j = 0; j < prevs.size(); j++){
                not_converged = std::min(not_converged,int(total_nonconverged(prevs[j],poses)));
            }
            //printf("%i -> not_convereturn MassFusionResults(poses,maxscore);rged: %i\n",i, not_converged);
			if(visualizationLvl == 3 && i == 0){

                printf("Noise: %f\n",surface_func->getNoise());
                for(unsigned int i = 0; i < edges.size(); i++){
                    for(unsigned int j = 0; j < edges[i].size(); j++){
                        if(edge_surface_func[i][j] != 0){
                            printf("%5.5f ",1.0*edge_surface_func[i][j]->getNoise());
                        }else{printf("%5.5f ",0.0);}
                    }
                    printf("\n");
                }

                printf("conv: %i\n",not_converged);
                for(unsigned int i = 0; i < edges.size(); i++){
                    for(unsigned int j = 0; j < edges[i].size(); j++){
                        if(edge_surface_func[i][j] != 0){
                            Eigen::Matrix4d dbefore = prev[j].inverse()*prev[i];
                            Eigen::Matrix4d dafter = poses[j].inverse()*poses[i];
                            double change = getChange(dbefore,dafter,edges[i][j]->node1->meandist + edges[i][j]->node2->meandist);
                            printf("%3.3f ",change / (convergence_mul*edge_surface_func[i][j]->getNoise()));
                        }else{printf("%3.3f ",0.0);}
                    }
                    printf("\n");
                }

				double maxscore = 0;
				for(unsigned int i = 0; i < edges.back().size(); i++){
					if(edges.back()[i] != 0){
						maxscore = std::max(maxscore,edges.back()[i]->score);
					}
				}


				double maxnoise = 0;

				show(poses,true);
				for(unsigned int i = 0; i < edge_surface_func.size(); i++){
					for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
						if(edge_surface_func[i][j] != 0){
							maxnoise = std::max(edge_surface_func[i][j]->getNoise(),maxnoise);
						}
					}
				}


				for(unsigned int i = 0; i < edge_surface_func.size(); i++){
					for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
						if(edge_surface_func[i][j] != 0 && edge_surface_func[i][j]->getNoise() == maxnoise){
							edges[i][j]->showMatchesW(viewer,edge_surface_func[i][j],poses[i],poses[j]);
						}
					}
				}

				edges.back().front()->showMatchesW(viewer,edge_surface_func.back().front(),poses.back(),poses.front());
			}

            prevs.push_back(poses);
            prev = poses;
            if(not_converged == 0){break;}
        }



        if(func_setup == 0){
            double surface_noise_before = surface_func->getNoise();
            surface_func->update();
            double surface_noise_after = surface_func->getNoise();
            if(fabs(1.0 - surface_noise_after/surface_noise_before) < 0.01){break;}
        }else if(func_setup == 1){
            int not_conv = 0;
            for(unsigned int i = 0; i < edge_surface_func.size(); i++){
                for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
                    if(i == j){continue;}
                    if(edge_surface_func[i][j] != 0){
                        double surface_noise_before = edge_surface_func[i][j]->getNoise();
                        edge_surface_func[i][j]->update();
                        double surface_noise_after = edge_surface_func[i][j]->getNoise();
                        if(fabs(1.0 - surface_noise_after/surface_noise_before) >= 0.01){not_conv++;}
                    }
                }
            }
            if(not_conv == 0){break;}
        }
    }

    for(int sz = 0; sz < 0; sz++){
        std::vector<Eigen::Matrix4d> prev2 = poses;
        printf("sz: %i\n",sz);
        int not_full = 0;
        for(unsigned int i = 0; i < nodes.size(); i++){
            nodes[i]->surface_nr_active = std::min(2*nodes[i]->surface_nr_active,nodes[i]->surface_nrp);
            not_full += nodes[i]->surface_nr_active != nodes[i]->surface_nrp;
        }
        printf("surface_nr_active: %i\n",nodes[0]->surface_nr_active);
        show(poses,true);

        for(int i = 0; i < 15; i++){
            rematch(poses,  i == 0);
            model(poses,    false);
            refine(poses,    true);
            double convergence = getConvergence();
            int not_converged = 0;
            for(unsigned int i = 0; i < nodes.size(); i++){
                double change = getChange(poses[i],prev[i],2.0*nodes[i]->meandist);
                if(visualizationLvl > 0){printf("%i -> %f\n",i,change);}
                if(change > convergence){not_converged++;}
            }
            if(visualizationLvl > 0){printf("not_converged: %i\n",not_converged);}
            prev = poses;
        }
        if(not_full == 0){break;}

        int not_converged = 0;
        double convergence = getConvergence();
        for(unsigned int i = 0; i < nodes.size(); i++){
            double change = getChange(poses[i],prev2[i],2.0*nodes[i]->meandist);
            printf("%i -> %f\n",i,change);
            if(change > convergence){not_converged++;}
        }
        //if(visualizationLvl > 0){printf("not_converged: %i\n",not_converged);}
        //printf("not_converged: %i\n",not_converged);
        prev2 = poses;
        if(not_converged == 0){break;}
    }


    double maxscore = 0;
    for(unsigned int i = 0; i < edges.back().size(); i++){
        if(edges.back()[i] != 0){
            maxscore = std::max(maxscore,edges.back()[i]->score);
        }
    }


	if(visualizationLvl > 0){show(poses,true);};
	if(visualizationLvl == 1){

		for(unsigned int i = 0; i < edge_surface_func.size(); i++){
			for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
				if(edge_surface_func[i][j] != 0){
					edge_surface_func[i][j]->debugg_print = true;
				}
			}
		}

		model(poses,true);

		printf("Noise: %f\n",surface_func->getNoise());
		for(unsigned int i = 0; i < edges.size(); i++){
			for(unsigned int j = 0; j < edges[i].size(); j++){
				if(edge_surface_func[i][j] != 0){
					printf("%5.5f ",1.0*edge_surface_func[i][j]->getNoise());
				}else{printf("%5.5f ",0.0);}
			}
			printf("\n");
		}
		printf("maxscore: %f\n",maxscore);

		edges.back().front()->showMatchesW(viewer,edge_surface_func.back().front(),poses.back(),poses.front());

		for(unsigned int i = 0; i < edge_surface_func.size(); i++){
			for(unsigned int j = 0; j < edge_surface_func[i].size(); j++){
				if(edge_surface_func[i][j] != 0){
					edge_surface_func[i][j]->debugg_print = false;
				}
			}
		}
	}


    //std::cout<<poses.back()<<std::endl;
    //show(poses,true);
    return MassFusionResults(poses,maxscore);
}

}