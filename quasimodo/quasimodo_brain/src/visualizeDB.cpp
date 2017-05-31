#include "ModelDatabase/ModelDatabase.h"
#include "ModelStorage/ModelStorage.h"
#include "Util/Util.h"

using namespace quasimodo_brain;
using namespace reglib;
using namespace std;
using namespace Eigen;

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
					p2.r = p.r;
					p2.g = p.g;
					p2.b = p.b;
					p2.x = p.x;
					p2.y = p.y;
                    p2.z = p.z;
                    X->points.push_back(p2);
                    xsum+= p.x;
                    ysum+= p.y;
                    zsum+= p.z;
                }
            }

			for(unsigned int k = 0; k < X->points.size(); k++){
				pcl::PointXYZRGB & p = X->points[k];
				p.x -= xsum/double(X->points.size());
				p.y -= ysum/double(X->points.size());
				p.z -= zsum/double(X->points.size());
				p.x += startx+j;
				p.y += 1.3*double(depth+1);
				p.z += 0;
			}

			viewer->addPointCloud<pcl::PointXYZRGB> (X, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(X), acc+"_frame_"+std::to_string(j));
            pcl::PointXYZ p;
			p.x = startx+j;
			p.y = 1.3*double(depth+1);
			p.z = 0;
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

    viewer->addPointCloud<pcl::PointXYZRGB> (cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(cld), acc);

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
}
class Compare
{
public:
	bool operator() (std::pair<double,std::string> a, std::pair<double,std::string> b){return a.first < b.first;}
};

int main(int argc, char **argv){
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
    viewer->setBackgroundColor(1.0,1.0,1.0);

	std::string filepath = "./";

	if(argc == 2){filepath = std::string(argv[1]);}

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

		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

        double startx = 0;
        double stopx = 0;
        addToViewer2(startx,stopx,0,viewer,model);
		viewer->spin();
		storage->fullHandback();
	}

	printf("dist = [");
	for(unsigned int i = 1; i < dist.size(); i++){printf("%i ",dist[i]);}
	printf("];\n");
    printf("nr_segments = %i mergescore: %i\n",segments,mergescore);

	return 0;
}
