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
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld (new pcl::PointCloud<pcl::PointXYZRGB> ());
    if(model->frames.size() == 0){
        model->recomputeModelPoints();
        cld	= model->getPCLcloud(1,false);
    }else{
        for(unsigned int j = 0; j < model->frames.size(); j++){
            bool * mask = model->modelmasks[j]->maskvec;
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr X (new pcl::PointCloud<pcl::PointXYZRGB> ());
            pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr frame_cld =  model->frames[j]->getPCLcloud();
            int r = rand()%256;
            int g = rand()%256;
            int b = rand()%256;
            for(unsigned int k = 0; k < frame_cld->points.size(); k++){
                pcl::PointXYZRGBNormal p = frame_cld->points[k];
                if(mask[k] && !std::isnan(p.x)){
                    pcl::PointXYZRGB p2;
                    p2.r = r;
                    p2.g = g;
                    p2.b = b;
                    p2.x = p.x;
                    p2.y = p.y;
                    p2.z = p.z;
                    X->points.push_back(p2);
                }
            }

            pcl::PointCloud<pcl::PointXYZRGB>::Ptr Xt (new pcl::PointCloud<pcl::PointXYZRGB> ());
            pcl::transformPointCloud (*X, *Xt, model->relativeposes[j]);
            *cld += *Xt;
        }
    }

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
    curr.x = x/double(cld->points.size());
    curr.y = y/double(cld->points.size());
    curr.z = z/double(cld->points.size());
	viewer->addLine<pcl::PointXYZ> (prev, curr,1,0,0,"line_"+acc);

	for(unsigned int i = 0; i < model->submodels.size(); i++){
        addToViewer(curr,viewer,model->submodels[i],acc+"_"+std::to_string(i),startx+i,starty+1.3,0);
	}
}


pcl::PointXYZ addToViewer2(double startx, double & stopx, int depth, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, reglib::Model * model, std::string acc = "model"){


    std::vector<pcl::PointXYZ > mids;
    for(unsigned int i = 0; i < model->submodels.size(); i++){
        mids.push_back(addToViewer2(stopx,stopx,depth+1,viewer,model->submodels[i],acc+"_"+std::to_string(i)));
    }


    std::vector<pcl::PointXYZ > framemids;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cld (new pcl::PointCloud<pcl::PointXYZRGB> ());
    if(model->submodels.size() == 0){//Leaf

        model->recomputeModelPoints();
        cld	= model->getPCLcloud(1,false);

        for(unsigned int j = 0; j < model->frames.size(); j++){
            stopx++;
            bool * mask = model->modelmasks[j]->maskvec;
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr X (new pcl::PointCloud<pcl::PointXYZRGB> ());
            pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr frame_cld =  model->frames[j]->getPCLcloud();
            int r = rand()%256;
            int g = rand()%256;
            int b = rand()%256;
            double xsum = 0;
            double ysum = 0;
            double zsum = 0;
            for(unsigned int k = 0; k < frame_cld->points.size(); k++){
                pcl::PointXYZRGBNormal p = frame_cld->points[k];
                if(mask[k] && !std::isnan(p.x)){
                    pcl::PointXYZRGB p2;
                    p2.r = r;
                    p2.g = g;
                    p2.b = b;
                    p2.x = p.x+startx+j;
                    p2.y = p.y+1.3*double(depth+1);
                    p2.z = p.z;
                    X->points.push_back(p2);
                    xsum+= p.x;
                    ysum+= p.y;
                    zsum+= p.z;
                }
            }

            for(unsigned int k = 0; k < X->points.size(); k++){
                pcl::PointXYZRGB p = X->points[k];
                    p.x -= xsum/double(X->points.size());
                    p.y -= ysum/double(X->points.size());
                    p.z -= zsum/double(X->points.size());
            }
            viewer->addPointCloud<pcl::PointXYZRGB> (cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld), acc+"_frame_"+std::to_string(j));
            pcl::PointXYZ p;
            p.x = startx+j-xsum/double(X->points.size());
            p.y = 1.3*double(depth+1)-ysum/double(X->points.size());
            p.z = zsum/double(X->points.size());
            framemids.push_back(p);
        }

    }else{
        model->recomputeModelPoints();
        cld	= model->getPCLcloud(1,false);
    }

    double midX = 0.5*(startx+stopx);
    for(unsigned int j = 0; j < cld->points.size(); j++){
        cld->points[j].x += midX;
        cld->points[j].y += 1.3*double(depth);
    }

    //printf("%i %f %f\n",depth,startx,stopx);

    viewer->addPointCloud<pcl::PointXYZRGB> (cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld), acc);
//    viewer->spin();

    pcl::PointXYZ p;
    p.x = midX;
    p.y = 1.3*double(depth);
    p.z = 0;


    for(unsigned int i = 0; i < mids.size(); i++){
        viewer->addLine<pcl::PointXYZ> (mids[i],p,1,0,0,"line_"+acc+"_"+std::to_string(i));
    }

    for(unsigned int i = 0; i < framemids.size(); i++){
        viewer->addLine<pcl::PointXYZ> (framemids[i],p,1,0,0,"frameline_"+acc+"_"+std::to_string(i));
    }


    return p;
/*

    if(model->frames.size() == 0){
        model->recomputeModelPoints();
        cld	= model->getPCLcloud(1,false);
    }else{
        for(unsigned int j = 0; j < model->frames.size(); j++){
            bool * mask = model->modelmasks[j]->maskvec;
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr X (new pcl::PointCloud<pcl::PointXYZRGB> ());
            pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr frame_cld =  model->frames[j]->getPCLcloud();
            int r = rand()%256;
            int g = rand()%256;
            int b = rand()%256;
            for(unsigned int k = 0; k < frame_cld->points.size(); k++){
                pcl::PointXYZRGBNormal p = frame_cld->points[k];
                if(mask[k] && !std::isnan(p.x)){
                    pcl::PointXYZRGB p2;
                    p2.r = r;
                    p2.g = g;
                    p2.b = b;
                    p2.x = p.x;
                    p2.y = p.y;
                    p2.z = p.z;
                    X->points.push_back(p2);
                }
            }

            pcl::PointCloud<pcl::PointXYZRGB>::Ptr Xt (new pcl::PointCloud<pcl::PointXYZRGB> ());
            pcl::transformPointCloud (*X, *Xt, model->relativeposes[j]);
            *cld += *Xt;
        }
    }

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
    curr.x = x/double(cld->points.size());
    curr.y = y/double(cld->points.size());
    curr.z = z/double(cld->points.size());
    viewer->addLine<pcl::PointXYZ> (prev, curr,1,0,0,"line_"+acc);

    for(unsigned int i = 0; i < model->submodels.size(); i++){
        addToViewer(curr,viewer,model->submodels[i],acc+"_"+std::to_string(i),startx+i,starty+1.3,0);
    }
    */
}
class Compare
{
public:
	bool operator() (std::pair<double,std::string> a, std::pair<double,std::string> b){return a.first < b.first;}
};

int main(int argc, char **argv){

	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
    //viewer->addCoordinateSystem(0.1);
    viewer->setBackgroundColor(1.0,1.0,1.0);

	std::string filepath = "./";

	if(argc == 2){
		filepath = std::string(argv[1]);
	}



	std::priority_queue<std::pair<double,std::string>, std::vector<std::pair<double,std::string>>, Compare> pq;

	ModelStorageFile * storage = new ModelStorageFile(filepath);
	storage->print();
	storage->loadAllModels();
	storage->print();

	for (std::map<std::string,std::string>::iterator it=storage->keyPathMap.begin(); it!=storage->keyPathMap.end(); ++it){
		reglib::Model * model = storage->fetch(it->first);
        pq.push(std::make_pair(model->getScore(5),it->first));
		storage->fullHandback();
	}

	std::vector<int> dist;
	dist.resize(100);
	for(unsigned int i = 0; i < dist.size(); i++){dist[i] = 0;}

    int segments = 0;
    int mergescore = 0;

	while(pq.size() > 0){
		std::pair<double,std::string> pair = pq.top();
		printf("score: %f\n",pair.first);
		pq.pop();
		reglib::Model * model = storage->fetch(pair.second);

		printf("keyval:%s\n",model->submodels.front()->keyval.c_str());
		int subs = model->submodels.size();
        segments += subs;
        mergescore += (subs-1)*(subs-1);
		if(model->frames.size() > 0){subs++;}
		dist[subs]+=1;

		pcl::PointXYZ curr;
		curr.x = 0;
		curr.y = 0;
		curr.z = 0;



		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

        double startx = 0;
        double stopx = 0;
        addToViewer2(startx,stopx,0,viewer,model);
        //addToViewer(curr,viewer,model);
		viewer->spin();
		storage->fullHandback();
	}

	printf("dist = [");
	for(unsigned int i = 1; i < dist.size(); i++){printf("%i ",dist[i]);}
	printf("];\n");
    printf("nr_segments = %i mergescore: %i\n",segments,mergescore);

	return 0;
}
