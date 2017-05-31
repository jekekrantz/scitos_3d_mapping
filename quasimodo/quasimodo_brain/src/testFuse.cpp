#include "ModelDatabase/ModelDatabase.h"
#include "ModelStorage/ModelStorage.h"
#include "Util/Util.h"
#include "reg2/RegistrationRandom2.cpp"
#include "mu2/ModelUpdater2.cpp"


using namespace quasimodo_brain;

int visualizationLvl = 0;

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

void print(std::vector<std::vector < float > > scores){
	printf("--------------------------\n");
	for(int i = 0; i < scores.size(); i++){
		for(int j = 0; j < scores.size(); j++){
			printf("%8.8i ",int(scores[i][j]));
		}
		printf("\n");
	}
}

Model * mergeModels(Model * model1, Model * model2, Matrix4d rp){
	if(model1->submodels.size() == 0 && model2->submodels.size() == 0){//both models are leaves, make new brach
		reglib::Model * newmodelHolder = new reglib::Model();
		model1->parrent = newmodelHolder;
		model2->parrent = newmodelHolder;

		newmodelHolder->submodels.push_back(model1);
		newmodelHolder->submodels_relativeposes.push_back(Eigen::Matrix4d::Identity());

		newmodelHolder->submodels.push_back(model2);
		newmodelHolder->submodels_relativeposes.push_back(rp);

		newmodelHolder->recomputeModelPoints();
		return newmodelHolder;
	}else if(model1->submodels.size() > model2->submodels.size()){

	}else{

	}
	return 0;
}

vector<int> getRepresentativeViews(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, std::vector< superpoint > spvec, vector<Matrix4d> cp, vector<RGBDFrame*> cf, double sum_olp = 0.9, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = 0){
	vector<int> ret;

	double nr_start = spvec.size();
	double removed = 0;
	vector<bool> taken;
	taken.resize(cp.size());
	for(unsigned int i = 0; i < cp.size(); i++){
		taken[i] = false;
	}

	while(true){
		int best = 0;
		std::vector<unsigned int> best_inds;
		for(unsigned int i = 0; i < cp.size(); i++){
			if(taken[i]){continue;}
			RGBDFrame * frame = cf[i];
			Matrix4d p = cp[i].inverse();

			unsigned int sum_olp = 0;
			std::vector<unsigned int> current_inds;

			unsigned short * dst_depthdata		= (unsigned short	*)(frame->depth.data);
			float		   * dst_normalsdata	= (float			*)(frame->normals.data);
			unsigned char  * dst_detdata        = (unsigned char    *)(frame->det_dilate.data);

			double m00 = p(0,0); double m01 = p(0,1); double m02 = p(0,2); double m03 = p(0,3);
			double m10 = p(1,0); double m11 = p(1,1); double m12 = p(1,2); double m13 = p(1,3);
			double m20 = p(2,0); double m21 = p(2,1); double m22 = p(2,2); double m23 = p(2,3);

			Camera * dst_camera				= frame->camera;
			const unsigned int dst_width	= dst_camera->width;
			const unsigned int dst_height	= dst_camera->height;
			const double dst_idepth			= dst_camera->idepth_scale;
			const double dst_cx				= dst_camera->cx;
			const double dst_cy				= dst_camera->cy;
			const double dst_fx				= dst_camera->fx;
			const double dst_fy				= dst_camera->fy;
			const double dst_ifx			= 1.0/dst_camera->fx;
			const double dst_ify			= 1.0/dst_camera->fy;
			const unsigned int dst_width2	= dst_camera->width  - 2;
			const unsigned int dst_height2	= dst_camera->height - 2;

			unsigned long nr_data = spvec.size();
			for(unsigned int src_ind = 0; src_ind < nr_data;++src_ind){
				const superpoint & sp = spvec[src_ind];

				double src_x = sp.x;
				double src_y = sp.y;
				double src_z = sp.z;
				double tz	= m20*src_x + m21*src_y + m22*src_z + m23;

				if(tz < 0){continue;}

				double tx	= m00*src_x + m01*src_y + m02*src_z + m03;
				double ty	= m10*src_x + m11*src_y + m12*src_z + m13;

				double itz	= 1.0/tz;
				double dst_w	= dst_fx*tx*itz + dst_cx;
				double dst_h	= dst_fy*ty*itz + dst_cy;

				if((dst_w > 0) && (dst_h > 0) && (dst_w < dst_width2) && (dst_h < dst_height2)){
					unsigned int dst_ind = unsigned(dst_h+0.5) * dst_width + unsigned(dst_w+0.5);

					double dst_z = dst_idepth*double(dst_depthdata[dst_ind]);
					double dst_nx = dst_normalsdata[3*dst_ind+0];
					if(dst_z > 0 && dst_nx != 2){
						if(dst_detdata[dst_ind] != 0){continue;}
						double dst_ny = dst_normalsdata[3*dst_ind+1];
						double dst_nz = dst_normalsdata[3*dst_ind+2];

						double dst_x = (double(int(0.5+dst_w)) - dst_cx) * dst_z * dst_ifx;
						double dst_y = (double(int(0.5+dst_h)) - dst_cy) * dst_z * dst_ify;

						double src_nx = sp.nx;
						double src_ny = sp.ny;
						double src_nz = sp.nz;

						double tnx	= m00*src_nx + m01*src_ny + m02*src_nz;
						double tny	= m10*src_nx + m11*src_ny + m12*src_nz;
						double tnz	= m20*src_nx + m21*src_ny + m22*src_nz;

						double residualZ = mysign(dst_z-tz)*fabs(tnx*(dst_x-tx) + tny*(dst_y-ty) + tnz*(dst_z-tz));//dst_z-tz;//mysign(dst_z-tz)*fabs(tnx*(dst_x-tx) + tny*(dst_y-ty) + tnz*(dst_z-tz));
						double angle = tnx*dst_nx + tny*dst_ny + tnz*dst_nz;


						double src_variance = 1.0/sp.point_information;
						double dst_variance = 1.0/getInformation(dst_z);
						double total_variance = src_variance+dst_variance;
						double total_stdiv = sqrt(total_variance);

						double d = residualZ/total_stdiv;

						double p_overlap_angle = nfunc->getProb(1-angle);
						double p_overlap = dfunc->getProb(d);

						p_overlap *= p_overlap_angle;

						if(p_overlap > 0.5){
							current_inds.push_back(src_ind);
						}
					}
				}
			}

			if(current_inds.size() > best_inds.size()){
				best = i;
				best_inds = current_inds;
			}
		}

		ret.push_back(best);
		taken[best] = true;

		double total_ratio = (double(best_inds.size())+removed)/nr_start;
		printf("%i total ratio: %f\n",best,total_ratio);
		if(total_ratio > sum_olp){break;}
		removed += best_inds.size();
		for(int i = best_inds.size()-1; i >= 0 ; i--){
			spvec[best_inds[i]] = spvec.back();
			spvec.pop_back();
		}

		if(cp.size() == ret.size()){break;}
	}




	return ret;
}


Model * mergeSubmodels(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * model, unsigned int max_size = 10, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = 0){
	printf("mergeSubmodels\n");
	if(model->submodels.size() <= max_size){return model;}

	reglib::RegistrationRandom2 *	reg			= new reglib::RegistrationRandom2(5);
	reglib::ModelUpdater2 * mu					= new reglib::ModelUpdater2();
	mu->occlusion_penalty						= 15;
	mu->viewer									= viewer;
	mu->show_scoring							= false;//fuse scoring show

	vector<Model *> models						= model->submodels;
	vector<Matrix4d> rps						= model->submodels_relativeposes;

	vector<vector < OcclusionScore > > ocs		= mu->computeOcclusionScore2(models,rps,1,false);
	std::vector<std::vector < float > > scores	= mu->getScores(ocs);

	print(scores);

//	viewer->removeAllPointClouds();
//	viewer->removeAllShapes();
//	double startx = 0;
//	double stopx = 0;
//	addToViewer2(startx,stopx,0,viewer,model);
//	viewer->spin();


//    unsigned int nr_models = models.size();
//    vector<vector < OcclusionScore > > occlusionScores;
//    occlusionScores.resize(nr_models);
//    for(unsigned int i = 0; i < nr_models; i++){
//        occlusionScores[i].resize(nr_models);
//        for(unsigned int j = 0; j < nr_models; j++){
//            occlusionScores[i][j] = OcclusionScore(0,0);
//        }
//	}

//    for(unsigned int i = 0; i < models.size(); i++){
//        Model * model1 = models[i];
//		for(unsigned int j = 0; j < models.size(); j++){
//            Model * model2 = models[j];
//			if(i != j){
//                Eigen::Matrix4d rp = rps[j].inverse() * rps[i];
//				occlusionScores[i][j] = computeOcclusionScore2(dfunc,nfunc,model1,model2, rp.inverse(), 1,debugg);
//				//printf("%i %i -> ",i,j); occlusionScores[i][j].print();
//            }
//        }
//	}
//    delete nfunc;
//    delete dfunc;

	while(scores.size() > max_size){
		float maxval = scores[0][0];
		unsigned int maxi = 0;
		unsigned int maxj = 0;
		for(unsigned int i = 0; i < scores.size(); i++){
			for(unsigned int j = i+1; j < scores.size(); j++){
				if(maxval < scores[i][j]){
					maxval = scores[i][j];
					maxi = i;
					maxj = j;
				}
			}
		}

		if(maxval < 0){break;}

		Matrix4d rpi = rps[maxi];
		Model * modeli = models[maxi];
		Matrix4d rpj = rps[maxj];
		Model * modelj = models[maxj];

		if(modeli->submodels.size() == 0 && modelj->submodels.size() == 0){//both models are leaves, make new brach
			printf("both are leaves\n");
			reglib::Model * newmodelHolder = new reglib::Model();
			modeli->parrent = newmodelHolder;
			modelj->parrent = newmodelHolder;

			newmodelHolder->submodels.push_back(modeli);
			newmodelHolder->submodels_relativeposes.push_back(Eigen::Matrix4d::Identity());
			newmodelHolder->rep_frames = modeli->rep_frames;
			newmodelHolder->rep_modelmasks = modeli->rep_modelmasks;
			newmodelHolder->rep_relativeposes = modeli->rep_relativeposes;


			Matrix4d m = rpi.inverse()*rpj;
			newmodelHolder->submodels.push_back(modelj);
			newmodelHolder->submodels_relativeposes.push_back(m);
//			for(unsigned int i = 0; i < modelj->rep_frames.size(); i++){
//				newmodelHolder->rep_frames.push_back(modelj->rep_frames[i]);
//				newmodelHolder->rep_modelmasks.push_back(modelj->rep_modelmasks[i]);
//				newmodelHolder->rep_relativeposes.push_back(m*modelj->rep_relativeposes[i]);
//			}
//			newmodelHolder->rep_frames.push_back(modelj->rep_frames);
//			newmodelHolder->rep_modelmasks.push_back(modelj->rep_modelmasks);
//			std::vector<Matrix4d> tmp = modelj->rep_relativeposes;
//			newmodelHolder->rep_relativeposes.push_back(tmp);

			newmodelHolder->recomputeModelPoints();

			//remove from old list
			int offseti = 0;
			for(unsigned i = 0; i < models.size(); i++){
				if(i != maxi && i != maxj){
					models[i-offseti] = models[i];
					rps[i-offseti] = rps[i];
				}else{
					offseti++;
				}
			}
			models.pop_back();
			models.pop_back();
			models.push_back(newmodelHolder);
			rps.pop_back();
			rps.pop_back();
			rps.push_back(rpi);
			//add att the end of old list
		}
		else if(modeli->submodels.size() > modeli->submodels.size()){
			modelj->parrent = modeli;
			modeli->submodels.push_back(modelj);
			modeli->submodels_relativeposes.push_back(rpi.inverse()*rpj);
			modeli->recomputeModelPoints();

			int offseti = 0;
			for(unsigned i = 0; i < models.size(); i++){
				if(i != maxi && i != maxj){
					models[i-offseti] = models[i];
					rps[i-offseti] = rps[i];
				}else{
					offseti++;
				}
			}
			models.pop_back();
			models.pop_back();
			models.push_back(modeli);
			rps.pop_back();
			rps.pop_back();
			rps.push_back(rpi);
		}else{
			modeli->parrent = modelj;
			modelj->submodels.push_back(modeli);
			modelj->submodels_relativeposes.push_back(rpj.inverse()*rpi);
			modelj->recomputeModelPoints();

			int offseti = 0;
			for(unsigned i = 0; i < models.size(); i++){
				if(i != maxi && i != maxj){
					models[i-offseti] = models[i];
					rps[i-offseti] = rps[i];
				}else{
					offseti++;
				}
			}
			models.pop_back();
			models.pop_back();
			models.push_back(modelj);
			rps.pop_back();
			rps.pop_back();
			rps.push_back(rpj);
		}


		model->submodels =  models;
		model->submodels_relativeposes = rps;

		ocs		= mu->computeOcclusionScore2(models,rps,1,false);
		scores	= mu->getScores(ocs);
		print(scores);
	}


	//mu->computeOcclusionScore2(models,rps,1,true);

	delete mu;
	delete reg;

	for(unsigned int i = 0; i < model->submodels.size(); i++){
		mergeSubmodels(dfunc,nfunc,model->submodels[i],max_size,viewer);
	}


	return model;
/*


			printf("models: %i\n",models.size());

	*/
/*

	print(scores);

	std::vector<std::vector < float > > scores_start = scores;
	std::vector < int > inds_start;
	std::vector < int > inds_new;
	while(scores.size() > max_size){
		float maxval = scores[0][0];
		unsigned int max_i = 0;
		unsigned int max_j = 0;
		for(unsigned int i = 0; i < scores.size(); i++){
			for(unsigned int j = i+1; j < scores.size(); j++){
				if(maxval < scores[i][j]){
					maxval = scores[i][j];
					max_i = i;
					max_j = j;
				}
			}
		}
		printf("%i %i -> %f\n",max_i,max_j,maxval);


		int nr_new = scores.size()-1;
		std::vector<std::vector < float > > scores_new;
		scores_new.resize(nr_new);
		for(unsigned int i = 0; i < nr_new; i++){scores_new[i].resize(nr_new);}

		int offseti = 0;
		for(int i = 0; i < scores.size(); i++){
			if(i != max_i && i != max_j){
				int offsetj = offseti;
				for(int j = i+1; j < scores.size(); j++){
					if(j != max_i && j != max_j){//Not part of the new object
						scores_new[i-offseti][j-offsetj] = scores[i][j];
					}else{
						offsetj++;
					}
				}
			}else{
				offseti++;
			}
		}
		scores = scores_new;
		print(scores_new);
	}
	exit(0);
	*/
}


int main(int argc, char **argv){


	DistanceWeightFunction2 * dfunc = new DistanceWeightFunction2();
	dfunc->f = THRESHOLD;
	dfunc->p = 0.005;

	DistanceWeightFunction2 * nfunc = new DistanceWeightFunction2();
	nfunc->f = THRESHOLD;
	nfunc->p = 0.50;


	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = boost::shared_ptr<pcl::visualization::PCLVisualizer>(new pcl::visualization::PCLVisualizer ("Modelserver Viewer"));
	viewer->setBackgroundColor(1.0,1.0,1.0);
//	int v1(0);
//	viewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
//	viewer->setBackgroundColor (0.5, 0.5, 0.5, v1);

//	int v2(0);
//	viewer->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
//	viewer->setBackgroundColor (0.7, 0.7, 0.7, v2);


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

	double startx = 0;
	double stopx = 0;

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


		for(int i = 0; i < model->submodels.size(); i++){
			printf("---------------------\n");
			vector<superpoint > spvec = model->submodels[i]->points;//, vector<Matrix4d> cp, vector<RGBDFrame*> cf, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = 0){
			vector<Matrix4d> cp		  = model->submodels[i]->relativeposes;
			vector<RGBDFrame*> cf	  = model->submodels[i]->frames;
			vector<int> rv = getRepresentativeViews(dfunc,nfunc,spvec,cp,cf, 0.95,viewer);
			for(int j = 0; j < rv.size(); j++){
				model->submodels[i]->rep_frames.push_back(model->submodels[i]->frames[rv[j]]);
				model->submodels[i]->rep_modelmasks.push_back(model->submodels[i]->modelmasks[rv[j]]);
				model->submodels[i]->rep_relativeposes.push_back(model->submodels[i]->relativeposes[rv[j]]);
			}
		}

//		viewer->removeAllPointClouds();
//		viewer->removeAllShapes();
//		double startx = 0;
//		double stopx = 0;
//		addToViewer2(startx,stopx,0,viewer,model);
//		viewer->spin();

/*
		reglib::RegistrationRandom2 *	reg	= new reglib::RegistrationRandom2(5);
		reglib::ModelUpdater2 * mu			= new reglib::ModelUpdater2();
		mu->occlusion_penalty               = 15;
		mu->viewer							= viewer;
		mu->show_scoring					= false;//fuse scoring show

		vector<Model *> models;
		vector<Matrix4d> rps;

		mu->addModelsToVector(models,rps,model,Eigen::Matrix4d::Identity());

		printf("models: %i\n",models.size());
		vector<vector < OcclusionScore > > ocs		= mu->computeOcclusionScore2(models,rps,1,false);
		mu->occlusion_penalty						= 15;
		std::vector<std::vector < float > > scores	= mu->getScores(ocs);

		float maxval = 0;
		for(int i = 0; i < scores.size(); i++){
			for(int j = 0; j < scores.size(); j++){
				maxval = std::max(maxval,fabs(scores[i][j]));
			}
		}

		for(int i = 0; i < scores.size(); i++){
			for(int j = 0; j < scores.size(); j++){
				printf("%8.8f ",scores[i][j]/maxval);
			}
			printf("\n");
		}
		std::vector<int> partition					= mu->getPartition(scores,2,5,2);
		printf("partition = [");
		for(int i = 0; i < partition.size(); i++){
			printf("%i ",partition[i]);
		}
		printf("];\n");

		delete mu;
		delete reg;
*/
		mergeSubmodels(dfunc,nfunc,model,5,viewer);

//		viewer->removeAllPointClouds();
//		viewer->removeAllShapes();

		addToViewer2(stopx,stopx,0,viewer,model,model->keyval);

		viewer->spin();
		storage->fullHandback();
	}

	printf("dist = [");
	for(unsigned int i = 1; i < dist.size(); i++){printf("%i ",dist[i]);}
	printf("];\n");
	printf("nr_segments = %i mergescore: %i\n",segments,mergescore);

	return 0;
}
