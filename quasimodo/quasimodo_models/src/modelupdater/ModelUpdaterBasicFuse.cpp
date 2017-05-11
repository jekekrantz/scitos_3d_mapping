#include "modelupdater/ModelUpdaterBasicFuse.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sys/time.h>
#include <emmintrin.h>

#include "registration/MassRegistration.h"

namespace reglib
{

//------------MODEL UPDATER-----------
ModelUpdaterBasicFuse::ModelUpdaterBasicFuse(Registration * registration_){
	registration = registration_;
	model = new reglib::Model();

    show_init_lvl = 0;//init show
    show_refine_lvl = 0;//refine show
    show_scoring = false;//fuse scoring show
}

ModelUpdaterBasicFuse::ModelUpdaterBasicFuse(Model * model_, Registration * registration_){
	registration = registration_;
	model = model_;

    show_init_lvl = 0;//init show
    show_refine_lvl = 0;//refine show
    show_scoring = false;//fuse scoring show
}

ModelUpdaterBasicFuse::~ModelUpdaterBasicFuse(){
	//printf("deleting ModelUpdaterBasicFuse\n");
}

FusionResults ModelUpdaterBasicFuse::registerModel(Model * model2, Eigen::Matrix4d guess, double uncertanity){

	if(model->points.size() > 0 && model2->points.size() > 0){
		registration->viewer	= viewer;
		registration->setDst(model->points);
		registration->setSrc(model2->points);

//        for(unsigned int i = 0; i < 2000; i++){
//			FusionResults fr1 = registration->getTransform(guess);
//			printf("%i :: score :: %f\n",i,fr1.scores.front());
//		}
//		exit(0);
		FusionResults fr = registration->getTransform(guess);
//		delete cd1;
//		delete cd2;

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

		for(unsigned int ca = 0; ca < todo && ca < 10; ca++){
			printf("ca: %i / %i \n",ca+1,todo);
			Eigen::Matrix4d pose = fr.candidates[ca];

			vector<Model *> models;
			vector<Matrix4d> rps;

			addModelsToVector(models,rps,model,Eigen::Matrix4d::Identity());
			unsigned int nr_models = models.size();
			addModelsToVector(models,rps,model2,pose);

			//Show alignment
			vector<vector < OcclusionScore > > ocs = computeOcclusionScore(models,rps,step,false);
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


			if(improvement > best){
//				for(unsigned int i = 0; i < scores.size(); i++){
//					for(unsigned int j = 0; j < scores.size(); j++){
//						if(scores[i][j] >= 0){printf(" ");}
//						printf("%5.5f ",0.00001*scores[i][j]);
//					}
//					printf("\n");
//				}
//				printf("partition "); for(unsigned int i = 0; i < partition.size(); i++){printf("%i ", partition[i]);} printf("\n");
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

void ModelUpdaterBasicFuse::fuse(Model * model2, Eigen::Matrix4d guess, double uncertanity){}


UpdatedModels ModelUpdaterBasicFuse::fuseData(FusionResults * f, Model * model1, Model * model2){
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

void ModelUpdaterBasicFuse::computeMassRegistration(std::vector<Eigen::Matrix4d> current_poses, std::vector<RGBDFrame*> current_frames,std::vector<cv::Mat> current_masks){
	printf("void ModelUpdaterBasicFuse::computeMassRegistration\n");
	printf("WARNING: THIS METHOD NOT IMPLEMENTED\n");
	exit(0);
}

void ModelUpdaterBasicFuse::setRegistration( Registration * registration_){
	if(registration != 0){delete registration;}
	registration = registration;
}

}


