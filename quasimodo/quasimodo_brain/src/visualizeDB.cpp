#include "ModelDatabase/ModelDatabase.h"
#include "ModelStorage/ModelStorage.h"
#include "Util/Util.h"

using namespace quasimodo_brain;
using namespace reglib;
using namespace std;
using namespace Eigen;

OcclusionScore addReprojections( DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Matrix4d p,
        vector< vector<double> > & distances, vector< vector<superpoint> > & points,
        RGBDFrame* frame, vector<superpoint> & spvec,
        boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, int debugg){

    double sum_olp = 0;
    double sum_ocl = 0;

    unsigned char * dst_detdata = (unsigned char*)(frame->det_dilate.data);

    Matrix4d pinv = p.inverse();

    std::vector<ReprojectionResult> rr_vec  = frame->getReprojections(spvec,p.inverse(),0,false);
    std::vector<superpoint> framesp			= frame->getSuperPoints();

    unsigned long nr_rr = rr_vec.size();
    for(unsigned long ind = 0; ind < nr_rr;ind++){
        ReprojectionResult & rr = rr_vec[ind];
        unsigned int src_ind = rr.src_ind;
        unsigned int dst_ind = rr.dst_ind;
        superpoint src_p = spvec[src_ind];

        if(src_p.z <= 0.0){continue;}

        superpoint dst_p = framesp[dst_ind];
        superpoint dp = dst_p;

        if(dst_detdata[dst_ind] != 0){continue;}

        dst_p.transform(p);
        points[src_ind].push_back(dst_p);
        distances[src_ind].push_back(rr.residualZ);

        double src_variance = 1.0/src_p.point_information;
        double dst_variance = 1.0/dst_p.point_information;
        double total_variance = src_variance+dst_variance;
        double total_stdiv = sqrt(total_variance);
        double rz = rr.residualZ;
        double d = rz/total_stdiv;
        double angle = rr.angle;

        double ocl = dfunc->getProbInfront(d);
        double olp = nfunc->getProb(1-angle) * dfunc->getProb(d);

        //double len = sqrt(dp.x*dp.x+dp.y*dp.y+dp.z*dp.z);
//        double reliable = fabs(dp.nx*dp.x+dp.ny*dp.y+dp.nz*dp.z)/len;
//        double ocl = reliable*(dfunc->getProbInfront(d) + (1-nfunc->getProb(1-angle))*(dfunc->getProb(d)));
//        double olp = reliable*nfunc->getProb(1-angle)*dfunc->getProb(d);

        sum_olp += olp;
        sum_ocl += ocl;
    }

    return OcclusionScore(sum_olp,sum_ocl);
}


OcclusionScore addResiduals(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc,
                  std::vector<superpoint> & spvec, Model * model, Eigen::Matrix4d pose,
                  vector< vector<double> > & distances, vector< vector<superpoint> > & points,
                  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, int debugg){
    if(debugg > 0){printf("addResiduals\n"); std::cout << pose << std::endl << std::endl;}

    OcclusionScore os;

    for(unsigned int i = 0; i < model->frames.size(); i++){
        if(debugg > 0){printf("relativeposes[%i]\n",i); std::cout << model->relativeposes[i] << std::endl << std::endl;}
        os.add( addReprojections(dfunc,nfunc,pose*model->relativeposes[i],distances,points,model->frames[i],spvec,viewer,debugg));
    }

    for(unsigned int i = 0; i < model->submodels.size(); i++){
        if(debugg > 0){printf("submodels_relativeposes[%i]\n",i); std::cout << model->submodels_relativeposes[i] << std::endl << std::endl;}
        os.add( addResiduals(dfunc,nfunc,spvec, model, pose*model->submodels_relativeposes[i],distances,points,viewer,debugg) );
    }
    return os;
}

OcclusionScore compare(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * model1, Model * model2, Matrix4d rp, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, int debugg){
    OcclusionScore os;

    std::vector<superpoint> spvec = model1->points;
    unsigned int nr_pixels = spvec.size();


    std::vector< std::vector<double> > distances;
    std::vector< std::vector<superpoint> > points;
    distances.resize(nr_pixels);
    points.resize(nr_pixels);

    os.add(addResiduals(dfunc,nfunc,spvec,model2,rp,distances,points,viewer,debugg));

    return os;
}

vector<vector < OcclusionScore > > computeOcclusionScore(vector<Model *> models, vector<Matrix4d> rps, int step, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, int debugg){
    unsigned int nr_models = models.size();
    vector<vector < OcclusionScore > > occlusionScores;
    occlusionScores.resize(nr_models);
    for(unsigned int i = 0; i < nr_models; i++){
        occlusionScores[i].resize(nr_models);
        for(unsigned int j = 0; j < nr_models; j++){
            occlusionScores[i][j] = OcclusionScore(0,0);
        }
    }

    DistanceWeightFunction2 * dfunc;
    DistanceWeightFunction2 * nfunc;

//    dfunc = new DistanceWeightFunction2();
//    dfunc->f = THRESHOLD;
//    dfunc->p = 0.005;

//    nfunc = new DistanceWeightFunction2();
//    nfunc->f = THRESHOLD;
//    nfunc->p = 0.50;

    dfunc = new DistanceWeightFunction2FILE("./dfunc.bin");
    nfunc = new DistanceWeightFunction2FILE("./nfunc.bin");

    for(unsigned int i = 0; i < models.size(); i++){
        for(unsigned int j = 0; j < models.size(); j++){
            if(i == j){continue;}
            if(debugg > 0){printf("compare: %i %i\n",i,j);}
            OcclusionScore os = compare(dfunc,nfunc,models[i],models[j], rps[i].inverse() * rps[j], viewer, debugg);
            occlusionScores[i][j].add(os);
            occlusionScores[j][i].add(os);
        }
    }

    delete dfunc;
    delete nfunc;

    return occlusionScores;
}

void getBox(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld, float & maxx, float & minx, float & maxy, float & miny, float & maxz, float & minz){
	maxx = cld->points.front().x;
	minx = cld->points.front().x;
	maxy = cld->points.front().y;
	miny = cld->points.front().y;
	maxz = cld->points.front().z;
	minz = cld->points.front().z;
	for(unsigned int i = 1; i < cld->points.size(); i++){
		maxx = std::max(cld->points[i].x,maxx);
		minx = std::min(cld->points[i].x,minx);

		maxy = std::max(cld->points[i].y,maxy);
		miny = std::min(cld->points[i].y,miny);

		maxz = std::max(cld->points[i].z,maxz);
		minz = std::min(cld->points[i].z,minz);
	}
}

void add(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld, float x, float y, float z){
	for(unsigned int i = 0; i < cld->points.size(); i++){
		cld->points[i].x += x;
		cld->points[i].y += y;
		cld->points[i].z += z;
	}
}

void addToViewer(pcl::PointXYZ prev, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, reglib::Model * model, std::string acc = "model", float startx = 0, float starty = 0, float startz = 0){
	model->recomputeModelPoints();
	printf("acc: %s\n",acc.c_str());
	printf("scores: %f %f %f %f\n",model->getScore(0),model->getScore(1),model->getScore(2),model->getScore(3));
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld	= model->getPCLcloud(1,false);
	float maxx, minx, maxy, miny, maxz, minz;
	getBox(cld,maxx,minx,maxy,miny,maxz, minz);
	add(cld,startx-minx,starty-miny,startz-minz);
	viewer->addPointCloud<pcl::PointXYZRGB> (cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld), acc);

	double x = 0;
	double y = 0;
	double z = 0;
	for(unsigned int i = 1; i < cld->points.size(); i++){
		x += cld->points[i].x;
		y += cld->points[i].y;
		z += cld->points[i].z;
	}

	pcl::PointXYZ curr;
	curr.x = x/double(cld->points.size());//+0.5*(maxx+minx);
	curr.y = y/double(cld->points.size());//+0.5*(maxy+miny);
	curr.z = z/double(cld->points.size());//+0.5*(maxz+minz);
	viewer->addLine<pcl::PointXYZ> (prev, curr,1,0,0,"line_"+acc);

	double prev_x = startx;
	for(unsigned int i = 0; i < model->submodels.size(); i++){
        addToViewer(curr,viewer,model->submodels[i],acc+"_"+std::to_string(i),startx+i,starty+1.3,0);
        //prev_x += 0.3+maxx-minx;
	}
}

int main(int argc, char **argv){

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
	viewer->addCoordinateSystem(0.1);
	viewer->setBackgroundColor(1.0,1.0,1.0);
//	viewer->addPointCloud<pcl::PointXYZRGB> (cld1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld1), "model1");
//	viewer->addPointCloud<pcl::PointXYZRGB> (cld2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld2), "model2");


	ModelStorageFile * storage = new ModelStorageFile();
	storage->print();
	storage->loadAllModels();
	storage->print();

	for (std::map<std::string,std::string>::iterator it=storage->keyPathMap.begin(); it!=storage->keyPathMap.end(); ++it){
		viewer->removeAllPointClouds();
		viewer->removeAllShapes();
		reglib::Model * model = storage->fetch(it->first);

        if(model->submodels.size() > 1 && model->getScore(3) > 1000 && model->getScore(2) > 250 && model->getScore(1) > 1000){
			pcl::PointXYZ curr;
			curr.x = 0;
			curr.y = 0;
			curr.z = 0;


            reglib::RegistrationRandom *	reg	= new reglib::RegistrationRandom(5);
            reglib::ModelUpdaterBasicFuse * mu	= new reglib::ModelUpdaterBasicFuse( model, reg);
            mu->occlusion_penalty               = 15;
            mu->viewer							= viewer;
            mu->show_scoring					= false;//fuse scoring show
            reg->visualizationLvl				= 0;

            std::vector<reglib::Model *> models;
            std::vector<Eigen::Matrix4d> rps;
            mu->addModelsToVector(models,rps,model,Eigen::Matrix4d::Identity());


            //Show alignment
            std::vector<std::vector < reglib::OcclusionScore > > ocs = mu->computeOcclusionScore(models,rps,1,mu->show_scoring);
            std::vector<std::vector < float > > scores = mu->getScores(ocs);
            std::vector<int> partition = mu->getPartition(scores,2,5,2);
            //std::vector<int> partition2 = getGroupings(scores, scores.size(),0,1,true);


            for(unsigned int i = 0; i < scores.size(); i++){
                for(unsigned int j = 0; j < scores.size(); j++){
                    if(scores[i][j] >= 0){printf(" ");}
                    printf("%5.5f ",0.00001*scores[i][j]);
                }
                printf("\n");
            }
            printf("partition  "); for(unsigned int i = 0; i < partition.size(); i++){printf("%i ", partition[i]);} printf("\n");
            //printf("partition2 "); for(unsigned int i = 0; i < partition2.size(); i++){printf("%i ", partition2[i]);} printf("\n");

            delete mu;
            delete reg;
            addToViewer(curr,viewer,model);
            viewer->spin();

//            viewer->removeAllPointClouds();
//            viewer->removeAllShapes();
//            //model->recomputeModelPoints(Matrix4d::Identity(),viewer);

//            vector<superpoint> spvec;
////            addAllSuperPoints( spvec,model,Matrix4d::Identity(),viewer);
////exit(0);

//            std::vector<std::vector < reglib::OcclusionScore > > ocs2 = computeOcclusionScore(models,rps,1, viewer, 0);
//            std::vector<std::vector < float > > scores2 = mu->getScores(ocs2);
//            printf("OCS2: \n");
//            for(unsigned int i = 0; i < scores2.size(); i++){
//                for(unsigned int j = 0; j < scores2.size(); j++){
//                    if(scores2[i][j] >= 0){printf(" ");}
//                    printf("%5.5f ",0.00001*scores2[i][j]);
//                }
//                printf("\n");
//            }

//            std::vector<ModelGrouping> gg = getGroupings(scores2, scores2.size(),0,1,true);

//            viewer->spin();
		}
		storage->fullHandback();
	}

	return 0;
}
