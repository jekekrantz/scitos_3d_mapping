#include "ModelUpdater2.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sys/time.h>
#include <emmintrin.h>

#include "registration/MassRegistration.h"

namespace reglib
{

//------------MODEL UPDATER-----------
ModelUpdater2::ModelUpdater2(Registration * registration_){
	registration = registration_;
	model = new reglib::Model();

    show_init_lvl = 0;//init show
    show_refine_lvl = 0;//refine show
    show_scoring = false;//fuse scoring show
}

ModelUpdater2::ModelUpdater2(Model * model_, Registration * registration_){
	registration = registration_;
	model = model_;

    show_init_lvl = 0;//init show
    show_refine_lvl = 0;//refine show
    show_scoring = false;//fuse scoring show
}

ModelUpdater2::~ModelUpdater2(){
	//printf("deleting ModelUpdater2\n");
}

FusionResults ModelUpdater2::registerModel(Model * model2, Eigen::Matrix4d guess, double uncertanity){

	if(model->points.size() > 0 && model2->points.size() > 0){
		registration->viewer	= viewer;

		std::vector<reglib::superpoint> model1_points = model->points;
		std::vector<reglib::superpoint> model2_points = model2->points;

//		Eigen::Affine3d startrot = Eigen::Affine3d::Identity();
//		startrot = Eigen::AngleAxisd(-30*2*M_PI/360.0, Eigen::Vector3d::UnitX());

//		for(unsigned long i = 0; i < model1_points.size(); i++){
//			model1_points[i].transform(startrot.matrix());
//		}

//		for(unsigned long i = 0; i < model2_points.size(); i++){
//			model2_points[i].transform(startrot.matrix());
//		}

		FusionResults fr;
		if(model->points.size() > model2->points.size() ){
			registration->setDst(model1_points);
			registration->setSrc(model2_points);
			fr = registration->getTransform(guess);
		}else{
			registration->setDst(model2_points);
			registration->setSrc(model1_points);
			fr = registration->getTransform(guess.inverse());
			for(unsigned int ca = 0; ca < fr.candidates.size(); ca++){
				fr.candidates[ca] = fr.candidates[ca].inverse();
			}
		}


        double best = -99999999999999;
        int best_id = -1;

		vector<Model *> testmodels;
		vector<Matrix4d> testrps;
		addModelsToVector(testmodels,testrps,model,Eigen::Matrix4d::Identity());
		addModelsToVector(testmodels,testrps,model2,Eigen::Matrix4d::Identity());

		int todo = fr.candidates.size();
		double expectedCost = double(todo)*computeOcclusionScoreCosts(testmodels);

		int step = 0.5 + expectedCost/11509168.5;// ~1 sec predicted
		step = std::max(1,step);

		for(unsigned int ca = 0; ca < todo && ca < 5; ca++){
			printf("ca: %i / %i -> %f \n",ca+1,todo,fr.scores[ca]);
			Eigen::Matrix4d pose = fr.candidates[ca];

			vector<Model *> models;
			vector<Matrix4d> rps;

			addModelsToVector(models,rps,model,Eigen::Matrix4d::Identity());
			unsigned int nr_models = models.size();
			addModelsToVector(models,rps,model2,pose);

			//Show alignment
			vector<vector < OcclusionScore > > ocs = computeOcclusionScore(models,rps,step,show_scoring);
			std::vector<std::vector < float > > scores = getScores(ocs);
			std::vector<int> partition = getPartition(scores,2,5,2);

			double sumscore1 = 0;
			for(unsigned int i = 0; i < models.size(); i++){
				for(unsigned int j = 0; j < models.size(); j++){
					if(i < nr_models && j < nr_models){sumscore1 += scores[i][j];}
					if(i >= nr_models && j >= nr_models){sumscore1 += scores[i][j];}
				}
			}

			double sumscore2 = 0;
			for(unsigned int i = 0; i < scores.size(); i++){
				for(unsigned int j = 0; j < scores.size(); j++){
					if(partition[i] == partition[j]){sumscore2 += scores[i][j];}
				}
			}

			double improvement = sumscore2-sumscore1;


//			printf("improvement: %10.10f ratio: %10.10f ",improvement*0.001,ocs[1][0].score/ocs[1][0].occlusions);
//			ocs[1][0].print();
////			computeOcclusionScore(models,rps,step,true);


//			if(improvement > 0){
//				//printf("improvement: %10.10f ",improvement*0.001);
//				//ocs[1][0].print();
//				//computeOcclusionScore(models,rps,step,true);
//			}
			if(show_scoring){
				printf("improvement: %10.10f ",improvement*0.001);
				ocs[1][0].print();
				for(unsigned int i = 0; i < scores.size(); i++){
					for(unsigned int j = 0; j < scores.size(); j++){
						if(scores[i][j] >= 0){printf(" ");}
						printf("%5.5f ",0.00001*scores[i][j]);
					}
					printf("\n");
				}
			}
			printf("partition "); for(unsigned int i = 0; i < partition.size(); i++){printf("%i ", partition[i]);} printf("\n");

			if(improvement > best){
				best = improvement;
				best_id = ca;
			}
		}

        if(best_id != -1){
            fr.score = 9999999;
            fr.guess = fr.candidates[best_id];
        }
		return fr;
	}
	return FusionResults();
}

void ModelUpdater2::fuse(Model * model2, Eigen::Matrix4d guess, double uncertanity){}


UpdatedModels ModelUpdater2::fuseData(FusionResults * f, Model * model1, Model * model2){
	UpdatedModels retval = UpdatedModels();
	Eigen::Matrix4d pose = f->guess;

//	printf("MODEL1 ");
//	model1->print();
//	printf("MODEL2 ");
//	model2->print();
//	printf("input pose\n");
//	std::cout << pose << std::endl << std::endl;

//	std::vector<Eigen::Matrix4d>	current_poses;
//	std::vector<RGBDFrame*>			current_frames;
//	std::vector<ModelMask*>			current_modelmasks;

//	for(unsigned int i = 0; i < model1->frames.size(); i++){
//		current_poses.push_back(				model1->relativeposes[i]);
//		current_frames.push_back(				model1->frames[i]);
//		current_modelmasks.push_back(			model1->modelmasks[i]);
//	}

//	for(unsigned int i = 0; i < model2->frames.size(); i++){
//		current_poses.push_back(	pose	*	model2->relativeposes[i]);
//		current_frames.push_back(				model2->frames[i]);
//		current_modelmasks.push_back(			model2->modelmasks[i]);
//	}

	vector<Model *> models;
	vector<Matrix4d> rps;
	addModelsToVector(models,rps,model,Eigen::Matrix4d::Identity());
	unsigned int nr_models1 = models.size();
	addModelsToVector(models,rps,model2,pose);

	double expectedCost = computeOcclusionScoreCosts(models);
//	printf("expectedCost: %f\n",expectedCost);
	int step = 0.5 + expectedCost/(10.0*11509168.5);// ~10 sec predicted max time
	step = std::max(1,step);

	vector<vector < OcclusionScore > > ocs = computeOcclusionScore(models,rps,step,show_scoring);
	std::vector<std::vector < float > > scores = getScores(ocs);
	std::vector<int> partition = getPartition(scores,2,5,2);

//	for(unsigned int i = 0; i < scores.size(); i++){
//		for(unsigned int j = 0; j < scores.size(); j++){
//			if(scores[i][j] >= 0){printf(" ");}
//			printf("%5.5f ",0.00001*scores[i][j]);
//		}
//		printf("\n");
//	}
//	printf("partition "); for(unsigned int i = 0; i < partition.size(); i++){printf("%i ", partition[i]);} printf("\n");


	double sumscore1 = 0;
	for(unsigned int i = 0; i < models.size(); i++){
		for(unsigned int j = 0; j < models.size(); j++){
			if(i < nr_models1 && j < nr_models1){sumscore1 += scores[i][j];}
			if(i >= nr_models1 && j >= nr_models1){sumscore1 += scores[i][j];}
		}
	}

	double sumscore2 = 0;
	for(unsigned int i = 0; i < scores.size(); i++){
		for(unsigned int j = 0; j < scores.size(); j++){
			if(partition[i] == partition[j]){sumscore2 += scores[i][j];}
		}
	}

	double improvement = sumscore2-sumscore1;

	for(unsigned int i = 0; i < scores.size(); i++){
		for(unsigned int j = 0; j < scores.size(); j++){
			if(scores[i][j] > 0){printf(" ");}
			printf("%5.5f ",0.0001*scores[i][j]);
		}
		printf("\n");
	}
	printf("partition "); for(unsigned int i = 0; i < partition.size(); i++){printf("%i ", partition[i]);} printf("\n");
	printf("sumscore before part: %f\n",sumscore1);
	printf("sumscore after  part: %f\n",sumscore2);
	printf("improvement:          %f\n",improvement);
//exit(0);
	std::vector<int> count;
	for(unsigned int i = 0; i < partition.size(); i++){
		if(int(count.size()) <= partition[i]){count.resize(partition[i]+1);}
		count[partition[i]]++;
	}

	int minpart = count[0];
	for(unsigned int i = 1; i < count.size(); i++){minpart = std::min(minpart,count[i]);}

	if(count.size() == 1){
		model1->merge(model2,pose);
		retval.updated_models.push_back(model1);
		retval.deleted_models.push_back(model2);
	}else if(improvement > 1){//Cannot fully fuse... separating...

		int c = 0;

		int model1_ind = partition.front();
		bool model1_same = true;
		for(unsigned int i = 0; i < model1->submodels.size(); i++){
			if(partition[c] != model1_ind){model1_same = false;}
			c++;
		}

		int model2_ind = partition.back();
		bool model2_same = true;
		for(unsigned int i = 0; i < model2->submodels.size(); i++){
			if(partition[c] != model2_ind){model2_same = false;}
			c++;
		}

		if(!model1_same || !model2_same){//If something changed, update models
			for(unsigned int i = 0; i < count.size(); i++){retval.new_models.push_back(new Model());}

			for(unsigned int i = 0; i < partition.size(); i++){
				retval.new_models[partition[i]]->submodels.push_back(models[i]);
				retval.new_models[partition[i]]->submodels_relativeposes.push_back(rps[i]);
				retval.new_models[partition[i]]->submodels.back()->parrent = retval.new_models[partition[i]];
			}

			for(unsigned int part = 0; part < retval.new_models.size(); part++){
				retval.new_models[part]->recomputeModelPoints();
			}

			retval.deleted_models.push_back(model1);
			retval.deleted_models.push_back(model2);
		}else{
			retval.unchanged_models.push_back(model1);
			retval.unchanged_models.push_back(model2);
		}
		return retval;
	}

	retval.unchanged_models.push_back(model1);
	retval.unchanged_models.push_back(model2);
	return retval;
}

void ModelUpdater2::computeMassRegistration(std::vector<Eigen::Matrix4d> current_poses, std::vector<RGBDFrame*> current_frames,std::vector<cv::Mat> current_masks){
	printf("void ModelUpdater2::computeMassRegistration\n");
	printf("WARNING: THIS METHOD NOT IMPLEMENTED\n");
	exit(0);
}

void ModelUpdater2::setRegistration( Registration * registration_){
	if(registration != 0){delete registration;}
	registration = registration;
}

OcclusionScore ModelUpdater2::computeOcclusionScore2(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, vector<superpoint> & spvec, Matrix4d p, RGBDFrame* frame, ModelMask* modelmask, int step,  bool debugg){
	OcclusionScore oc;

	double sum_ocl = 0;
	double sum_olp = 0;

	unsigned char  * dst_maskdata		= (unsigned char	*)(modelmask->mask.data);
	unsigned char  * dst_rgbdata		= (unsigned char	*)(frame->rgb.data);
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
	for(unsigned long src_ind = 0; src_ind < nr_data;++src_ind){
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
				double angle = tnx*dst_nx + tny*dst_ny + tnz*dst_nz;//dst_z-tz;//


				double src_variance = 1.0/sp.point_information;
				double dst_variance = 1.0/getInformation(dst_z);
				double total_variance = src_variance+dst_variance;
				double total_stdiv = sqrt(total_variance);

				double d = residualZ/total_stdiv;
				double d2 = sqrt((tx-dst_x)*(tx-dst_x) + (ty-dst_y)*(ty-dst_y) + (tz-dst_x)*(tz-dst_z))/total_stdiv;

				double p_overlap_angle = nfunc->getProb(1-angle);
				double p_overlap = std::min(dfunc->getProb(d),dfunc->getProb(0.1*d2));
				double p_occlusion = dfunc->getProbInfront(d);

				p_overlap *= p_overlap_angle;

				sum_ocl += p_occlusion;
				sum_olp += p_overlap;
			}
		}
	}

	sum_olp = std::max(sum_olp-0.05*nr_data,0.0);

	return OcclusionScore(sum_olp,sum_ocl);
}

OcclusionScore ModelUpdater2::computeOcclusionScore2(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * mod, vector<Matrix4d> cp, vector<RGBDFrame*> cf, vector<ModelMask*> cm, Matrix4d rp, int step,  bool debugg){
    OcclusionScore ocs;

    for(unsigned int i = 0; i < cp.size(); i++){
		ocs.add(computeOcclusionScore2(dfunc,nfunc,mod->points,cp[i].inverse()*rp,cf[i],cm[i],step,debugg));
    }
    return ocs;
}


OcclusionScore ModelUpdater2::computeOcclusionScore2(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * model1, Model * model2, Matrix4d rp, int step, bool debugg){
    OcclusionScore ocs;
	if(model2->rep_frames.size() > 0){
		ocs.add(computeOcclusionScore2(dfunc,nfunc,model1, model2->rep_relativeposes,model2->rep_frames,model2->rep_modelmasks,rp.inverse(),step,debugg));
	}else if(model2->frames.size() > 0){
		ocs.add(computeOcclusionScore2(dfunc,nfunc,model1, model2->relativeposes,model2->frames,model2->modelmasks,rp.inverse(),step,debugg));
	}else{
		for(unsigned int i = 0; i < model2->submodels.size(); i++){
		//	ocs.add(computeOcclusionScore2(dfunc,nfunc,model1,model2->submodels[i],model2->submodels_relativeposes[i] * rp,step,debugg));
		}
	}
	//    ocs.add(computeOcclusionScore(dfunc,nfunc,model2, model1->relativeposes,model1->frames,model1->modelmasks,rp,step,debugg));
    return ocs;
}

vector<vector < OcclusionScore > > ModelUpdater2::computeOcclusionScore2(vector<Model *> models, vector<Matrix4d> rps, int step, bool debugg){

	if(debugg){
        viewer->removeAllPointClouds();
        char buf [1024];
        for(unsigned int i = 0; i < models.size(); i++){
            sprintf(buf,"model%i",i);
            pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X = getPointCloudFromVector(models[i]->points,3,rand()%256,rand()%256,rand()%256);
            pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Xt (new pcl::PointCloud<pcl::PointXYZRGBNormal> ());
            pcl::transformPointCloud (*X, *Xt, rps[i].cast<float>());

            viewer->addPointCloud<pcl::PointXYZRGBNormal> (Xt, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Xt), buf);
            viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, buf);
        }
        viewer->spin();
        viewer->removeAllPointClouds();
/*
		for(unsigned int i = 0; i < models.size(); i++){
			Model * model1 = models[i];
			printf("model: %i\n",i);
			viewer->removeAllPointClouds();
			char buf [1024];
			for(unsigned int j = 0; j < model1->frames.size(); j++){
				sprintf(buf,"frame%i",j);

				bool * mask = model1->modelmasks[j]->maskvec;
				pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X (new pcl::PointCloud<pcl::PointXYZRGBNormal> ());
				pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr frame_cld =  model1->frames[j]->getPCLcloud();
				int r = rand()%256;
				int g = rand()%256;
				int b = rand()%256;
				for(unsigned int k = 0; k < frame_cld->points.size(); k++){
					pcl::PointXYZRGBNormal p = frame_cld->points[k];
					if(mask[k] && !std::isnan(p.x)){
						p.r = r;
						p.g = g;
						p.b = b;
						X->points.push_back(p);
					}
				}

				pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Xt (new pcl::PointCloud<pcl::PointXYZRGBNormal> ());
				pcl::transformPointCloud (*X, *Xt, model1->relativeposes[j]);

				viewer->addPointCloud<pcl::PointXYZRGBNormal> (Xt, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Xt), buf);
				viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, buf);
			}
			viewer->spin();
			viewer->removeAllPointClouds();

			for(unsigned int j = 0; j < models.size(); j++){
				Model * model2 = models[j];
				if(i != j){
					Eigen::Matrix4d rp = rps[j].inverse() * rps[i];

					printf("models: %i %i\n",i,j);
					viewer->removeAllPointClouds();
					char buf [1024];

					sprintf(buf,"model%i",j);
					pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X = getPointCloudFromVector(model2->points,3,255,0,0);
					viewer->addPointCloud<pcl::PointXYZRGBNormal> (X, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X), buf);
					viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, buf);

					sprintf(buf,"model%i",i);
					pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X0 = getPointCloudFromVector(model1->points,3,0,255,0);
					pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Xt (new pcl::PointCloud<pcl::PointXYZRGBNormal> ());
					pcl::transformPointCloud (*X0, *Xt, rp);

					viewer->addPointCloud<pcl::PointXYZRGBNormal> (Xt, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Xt), buf);
					viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, buf);

					viewer->spin();
					viewer->removeAllPointClouds();
				}
			}
		}
		*/
	}

    DistanceWeightFunction2 * dfunc = new DistanceWeightFunction2();
    dfunc->f = THRESHOLD;
    dfunc->p = 0.005;

    DistanceWeightFunction2 * nfunc = new DistanceWeightFunction2();
    nfunc->f = THRESHOLD;
    nfunc->p = 0.50;

    unsigned int nr_models = models.size();
    vector<vector < OcclusionScore > > occlusionScores;
    occlusionScores.resize(nr_models);
    for(unsigned int i = 0; i < nr_models; i++){
        occlusionScores[i].resize(nr_models);
        for(unsigned int j = 0; j < nr_models; j++){
            occlusionScores[i][j] = OcclusionScore(0,0);
        }
	}

    for(unsigned int i = 0; i < models.size(); i++){
        Model * model1 = models[i];
		for(unsigned int j = 0; j < models.size(); j++){
            Model * model2 = models[j];
			if(i != j){
                Eigen::Matrix4d rp = rps[j].inverse() * rps[i];
				occlusionScores[i][j] = computeOcclusionScore2(dfunc,nfunc,model1,model2, rp.inverse(), 1,debugg);
				//printf("%i %i -> ",i,j); occlusionScores[i][j].print();
            }
        }
	}

    delete nfunc;
    delete dfunc;
    return occlusionScores;
}


float ModelUpdater2::recursive_split2(std::vector<Graph*> * graphs_out,std::vector<std::vector<int>> * graphinds_out, Graph * graph, std::vector<int> graph_inds){
	printf("%s::%i\n",__PRETTY_FUNCTION__,__LINE__);
	if(boost::num_vertices(*graph) == 1){
		graphs_out->push_back(graph);
		graphinds_out->push_back(graph_inds);
		return 0;
	}

	std::vector<Graph*> second_graphs;
	std::vector<std::vector<int>> second_graphinds;
	float w = graph_cut(second_graphs,second_graphinds,*graph,graph_inds);
	printf("line %i -> %f\n",__LINE__,w);
	if(w <= 0){
		delete graph;
		return 2*w + recursive_split2(graphs_out, graphinds_out,second_graphs.front(),second_graphinds.front()) + recursive_split2(graphs_out, graphinds_out, second_graphs.back(),second_graphinds.back());
	}else{
		graphs_out->push_back(graph);
		graphinds_out->push_back(graph_inds);
		delete second_graphs.front();
		delete second_graphs.back();
		return 0;
	}
}

std::vector<int> ModelUpdater2::partition_graph2(std::vector< std::vector< float > > & scores){
	printf("%s::%i\n",__PRETTY_FUNCTION__,__LINE__);
	printf("size: %i\n",scores.size());

	int nr_data = scores.size();
	Graph* graph = new Graph(nr_data);
	std::vector<int> graph_inds;
	graph_inds.resize(nr_data);

	typename boost::property_map<Graph, boost::vertex_name_t>::type vertex_name = boost::get(boost::vertex_name, *graph);

	float sum = 0;
	for(int i = 0; i < nr_data; i++){
		graph_inds[i] = i;
		for(int j = i+1; j < nr_data; j++){
			float weight = scores[i][j];
			if(weight != 0){
				sum += 2*weight;
				edge_weight_property e = weight;
				boost::add_edge(i, j, e, *graph);
			}
		}
	}

	std::vector<Graph*> * graphs_out = new std::vector<Graph*>();
	std::vector<std::vector<int>> * graphinds_out = new std::vector<std::vector<int>>();

	printf("%s::%i\n",__PRETTY_FUNCTION__,__LINE__);
	if(boost::num_vertices(*graph) == 1){
		graphs_out->push_back(graph);
		graphinds_out->push_back(graph_inds);
		printf("exit at line %i\n",__LINE__);
		exit(0);
		//return 0;
	}

	std::vector<Graph*> second_graphs;
	std::vector<std::vector<int>> second_graphinds;
	float w = graph_cut(second_graphs,second_graphinds,*graph,graph_inds);
	printf("line %i -> %f\n",__LINE__,w);


	std::vector<int> part;
	part.resize(nr_data);
	for(unsigned int i = 0; i < graphinds_out->size(); i++){
		for(unsigned int j = 0; j < graphinds_out->at(i).size(); j++){
			part[graphinds_out->at(i).at(j)] = i;
		}
	}

	printf("partition = [");
	for(unsigned int i = 0; i < nr_data; i++){
		printf("%i ");
	}
	printf("];\n");


//	if(w <= 0){
//		delete graph;
//		return 2*w + recursive_split2(graphs_out, graphinds_out,second_graphs.front(),second_graphinds.front()) + recursive_split2(graphs_out, graphinds_out, second_graphs.back(),second_graphinds.back());
//	}else{
//		graphs_out->push_back(graph);
//		graphinds_out->push_back(graph_inds);
//		delete second_graphs.front();
//		delete second_graphs.back();
//		return 0;
//	}

//	float best = sum-recursive_split2(graphs_out,graphinds_out, graph,graph_inds );

//	std::vector<int> part;
//	part.resize(nr_data);
//	for(unsigned int i = 0; i < graphinds_out->size(); i++){
//		for(unsigned int j = 0; j < graphinds_out->at(i).size(); j++){
//			part[graphinds_out->at(i).at(j)] = i;
//		}
//	}

	return std::vector<int>();

//	int nr_data = scores.size();
//	Graph* graph = new Graph(nr_data);
//	std::vector<int> graph_inds;
//	graph_inds.resize(nr_data);

//	typename boost::property_map<Graph, boost::vertex_name_t>::type vertex_name = boost::get(boost::vertex_name, *graph);

//	float sum = 0;
//	for(int i = 0; i < nr_data; i++){
//		graph_inds[i] = i;
//		for(int j = i+1; j < nr_data; j++){
//			float weight = scores[i][j];
//			if(weight != 0){
//				sum += 2*weight;
//				edge_weight_property e = weight;
//				boost::add_edge(i, j, e, *graph);
//			}
//		}
//	}

//	std::vector<Graph*> * graphs_out = new std::vector<Graph*>();
//	std::vector<std::vector<int>> * graphinds_out = new std::vector<std::vector<int>>();
//	float best = sum-recursive_split2(graphs_out,graphinds_out, graph,graph_inds );

//	std::vector<int> part;
//	part.resize(nr_data);
//	for(unsigned int i = 0; i < graphinds_out->size(); i++){
//		for(unsigned int j = 0; j < graphinds_out->at(i).size(); j++){
//			part[graphinds_out->at(i).at(j)] = i;
//		}
//	}
//	return part;
}

std::vector<int> ModelUpdater2::getPartition2(std::vector< std::vector< float > > & scores){
	printf("%s::%i\n",__PRETTY_FUNCTION__,__LINE__);
	if(scores.size() < 20){
		printf("getGroupings2 \n");
		std::vector<int> part = getGroupings2(scores,scores.size());
		printf("partition = [");
		for(unsigned int i = 0; i < part.size(); i++){
			printf("%i ",part[i]);
		}
		printf("];\n");

		return part;
	}else{
		return partition_graph2(scores);
	}
}

}
