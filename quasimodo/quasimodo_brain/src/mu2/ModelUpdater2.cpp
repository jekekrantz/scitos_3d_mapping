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

vector<vector < OcclusionScore > > ModelUpdater2::computeOcclusionScore(vector<Model *> models, vector<Matrix4d> rps, int step, bool debugg){

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
	}

	DistanceWeightFunction2 * dfunc;
	DistanceWeightFunction2 * nfunc;
	dfunc = new DistanceWeightFunction2();
	dfunc->f = THRESHOLD;
	dfunc->p = 0.005;

	nfunc = new DistanceWeightFunction2();
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

	double startTime = getTime();
/*
	std::vector<double> dvec;
	std::vector<double> nvec;

	int tot_nr_pixels = 0;
	std::vector<int> offsets;

	for(unsigned int i = 0; i < models.size(); i++){
		offsets.push_back(tot_nr_pixels);
		unsigned int nr_pixels = models[i]->points.size();
		tot_nr_pixels += nr_pixels;
	}


	for(unsigned int i = 0; i < models.size(); i++){
		Model * model1 = models[i];
		for(unsigned int j = 0; j < models.size(); j++){
			if(i == j){continue;}
			Model * model2 = models[j];
			for(unsigned int k = 0; k < model2->relativeposes.size(); k++){
				Eigen::Matrix4d p;
				p = model2->relativeposes[k].inverse()*(rps[j].inverse() * rps[i]);
				testgetDynamicWeights(true,dvec,nvec,dfunc,nfunc,p, model1->points,0,0,0,model2->frames[k]);
			}
		}
	}

	double dstdval = 0;
	for(unsigned int i = 0; i < dvec.size(); i++){dstdval += dvec[i]*dvec[i];}
	dstdval = sqrt(dstdval/double(dvec.size()-1));



	double noiseWeight = dfunc->getNoise();
*/
	for(unsigned int i = 0; i < models.size(); i++){
		Model * model1 = models[i];
/*
		unsigned int nr_pixels = model1->points.size();
		double * current_overlaps		= new double[nr_pixels];
		double * current_occlusions		= new double[nr_pixels];
		double * current_notocclusions	= new double[nr_pixels];

		double * bg_overlaps			= new double[nr_pixels];
		double * bg_occlusions			= new double[nr_pixels];
		double * bg_notocclusions		= new double[nr_pixels];

		for(unsigned int m = 0; m < nr_pixels; m++){
			bg_overlaps[m] = 0;
			bg_occlusions[m] = 0.0;
			bg_notocclusions[m] = 0.0;
		}

		if(debugg){
			printf("SELF OCCLUSION\n");
		}
		for(unsigned int k = 0; k < model1->relativeposes.size(); k++){
			Eigen::Matrix4d p = model1->relativeposes[k].inverse();
			//testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,bg_overlaps, bg_occlusions, bg_notocclusions,model1->frames[k],debugg);
			testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,bg_overlaps, bg_occlusions, bg_notocclusions,model1->frames[k],false);
		}

		for(unsigned int m = 0; m < nr_pixels; m++){
			bg_occlusions[m]		= std::min(0.99999,bg_occlusions[m]/std::max(1.0,bg_notocclusions[m]));
		}

		for(unsigned int j = 0; j < models.size(); j++){

			if(i == j){continue;}

			if(debugg){
				printf("COMPARE %i to %i\n",i,j);
			}

			for(unsigned int m = 0; m < nr_pixels; m++){
				current_overlaps[m] = 0;
				current_occlusions[m] = 0.0;
				current_notocclusions[m] = 0.0;
			}

			Model * model2 = models[j];
			//For all frames
			for(unsigned int k = 0; k < model2->relativeposes.size(); k++){
				Eigen::Matrix4d p;
				//p = model2->relativeposes[k].inverse()*(rps[i].inverse() * rps[j]);
				//testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,current_overlaps, current_occlusions, current_notocclusions,model2->frames[k],debugg);

				p = model2->relativeposes[k].inverse()*(rps[j].inverse() * rps[i]);
				testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,current_overlaps, current_occlusions, current_notocclusions,model2->frames[k],debugg);

				//p = model2->relativeposes[k].inverse()*(rps[i] * rps[j].inverse());
				//testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,current_overlaps, current_occlusions, current_notocclusions,model2->frames[k],debugg);

				//p = model2->relativeposes[k].inverse()*(rps[j].inverse() * rps[i].inverse());
				//testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,current_overlaps, current_occlusions, current_notocclusions,model2->frames[k],debugg);
			}

			for(unsigned int m = 0; m < nr_pixels; m++){
				current_occlusions[m]	= std::min(0.99999,current_occlusions[m]/std::max(1.0,current_notocclusions[m]));
			}

			double occlusion_count = 0;
			double overlap_count = 0;

			for(unsigned int ind = 0; ind < nr_pixels; ind++){
				double ocl = std::max(0.0,current_occlusions[ind]-bg_occlusions[ind]);
				double olp = std::min(1-ocl,current_overlaps[ind]);

				occlusion_count += ocl;
				overlap_count += olp;
			}

			OcclusionScore oc (overlap_count/noiseWeight,occlusion_count/noiseWeight);
			occlusionScores[i][j].add(oc);
			occlusionScores[j][i].add(oc);

			if(debugg){
				printf("dfunc->getNoise() = %f\n",dfunc->getNoise());
				pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cld = getPointCloudFromVector(model1->points);
				for(unsigned int ind = 0; ind < nr_pixels; ind++){
					double ocl = std::max(0.0,current_occlusions[ind]-bg_occlusions[ind]);
					double olp = std::min(1-ocl,current_overlaps[ind]);
					cld->points[ind].r = 255*ocl;
					cld->points[ind].g = 255*olp;
					cld->points[ind].b = 0;
				}
				viewer->setBackgroundColor(1.0,1.0,1.0);
				viewer->removeAllPointClouds();
				viewer->addPointCloud<pcl::PointXYZRGBNormal> (cld, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(cld), "cld");
				viewer->spin();
			}
		}

		delete[] current_overlaps;
		delete[] current_occlusions;
		delete[] current_notocclusions;

		delete[] bg_overlaps;
		delete[] bg_occlusions;
		delete[] bg_notocclusions;
		*/
	}

	delete nfunc;
	delete dfunc;
	return occlusionScores;
}

}


