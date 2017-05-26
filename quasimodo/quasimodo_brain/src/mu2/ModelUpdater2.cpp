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

OcclusionScore ModelUpdater2::computeOcclusionScore(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, vector<superpoint> & spvec, Matrix4d p, RGBDFrame* frame, ModelMask* modelmask, int step,  bool debugg){
    OcclusionScore oc;

    Matrix4d invp = p.inverse();



    unsigned char  * dst_maskdata		= (unsigned char	*)(modelmask->mask.data);
    unsigned char  * dst_rgbdata		= (unsigned char	*)(frame->rgb.data);
    unsigned short * dst_depthdata		= (unsigned short	*)(frame->depth.data);
    float		   * dst_normalsdata	= (float			*)(frame->normals.data);
    unsigned char  * dst_detdata        = (unsigned char    *)(frame->det_dilate.data);

    //std::vector<superpoint> framesp = cf->getSuperPoints();

    std::vector<ReprojectionResult> rr_vec	= frame->getReprojections(spvec,p,modelmask->maskvec,false);
    std::vector<superpoint> framesp			= frame->getSuperPoints(p);
    std::vector<superpoint> framesp1		= frame->getSuperPoints(p);
    std::vector<superpoint> framesp2		= frame->getSuperPoints(invp);
    unsigned long nr_rr = rr_vec.size();

    if(debugg){
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();



        vector<superpoint> spvec2 = spvec;
        for(unsigned int i = 0; i < spvec2.size(); i++){
//            spvec2[i].y -= 0.2*3.0;
//            spvec2[i].x += 0.15*3.0;
            spvec2[i].r = 255;
            spvec2[i].g =   0;
            spvec2[i].b = 255;
        }

        for(unsigned long ind = 0; ind < nr_rr; ind+= 1){//std::max(1,int(nr_rr/300))){

            char buf [1024];
            sprintf(buf,"line%i",ind);
            ReprojectionResult & rr = rr_vec[ind];
            superpoint & src_p  =   spvec[rr.src_ind];
            superpoint & dst_p  = framesp2[rr.dst_ind];
            if(src_p.z < 0.0){continue;}

            double src_variance = 1.0/src_p.point_information;
            double dst_variance = 1.0/dst_p.point_information;
            double total_variance = src_variance+dst_variance;
            double total_stdiv = sqrt(total_variance);
            double rz = rr.residualZ;
            double d = rz/total_stdiv;
            double surface_angle = rr.angle;

            double p_overlap_angle = nfunc->getProb(1-surface_angle);
            double p_overlap = dfunc->getProb(d);
            double p_occlusion = dfunc->getProbInfront(d);

            p_overlap *= p_overlap_angle;

            double p_other = 1.0 - p_overlap - p_occlusion;

            superpoint & src_p2 =   spvec2[rr.src_ind];
            src_p2.r = 255*p_occlusion;
            src_p2.g = 255*p_overlap;
            src_p2.b = 255*p_other;

            //viewer->addLine<pcl::PointXYZRGBNormal>(X2->points[rr.src_ind],Y2->points[rr.dst_ind],255.0*p_occlusion,255.0*p_overlap,255.0*p_other, buf,2);
        }

        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X1 = getPointCloudFromVector(spvec,    3,255,0,0);
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr X2 = getPointCloudFromVector(spvec2,   0,255,0,0);
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y1 = getPointCloudFromVector(framesp2, 0,255,0,0);
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr Y2 = getPointCloudFromVector(framesp2, 0,255,0,0);
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (X1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X1), "model1",1);
        viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "model1",1);
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (X2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(X2), "model2",2);
        viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "model2",2);
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y1, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y1), "frame1",1);
        viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "frame1",1);
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (Y2, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(Y2), "frame2",2);
        viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "frame2",2);
        viewer->spin();
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();
    }

    return oc;
}

OcclusionScore ModelUpdater2::computeOcclusionScore(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * mod, vector<Matrix4d> cp, vector<RGBDFrame*> cf, vector<ModelMask*> cm, Matrix4d rp, int step,  bool debugg){
    OcclusionScore ocs;
    for(unsigned int i = 0; i < cp.size(); i++){
        printf("frame %i\n",i);
        ocs.add(computeOcclusionScore(dfunc,nfunc,mod->points,cp[i].inverse()*rp,cf[i],cm[i],step,debugg));
    }
    return ocs;
}


OcclusionScore ModelUpdater2::computeOcclusionScore(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * model1, Model * model2, Matrix4d rp, int step, bool debugg){
    OcclusionScore ocs;
    ocs.add(computeOcclusionScore(dfunc,nfunc,model1, model2->relativeposes,model2->frames,model2->modelmasks,rp.inverse(),step,debugg));
//    ocs.add(computeOcclusionScore(dfunc,nfunc,model2, model1->relativeposes,model1->frames,model1->modelmasks,rp,step,debugg));
    return ocs;
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

        if(debugg){
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
        }

        for(unsigned int j = 0; false && j < models.size(); j++){
            Model * model2 = models[j];
            if(i != j){
                Eigen::Matrix4d rp = rps[j].inverse() * rps[i];

                if(debugg){
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
                OcclusionScore os = computeOcclusionScore(dfunc,nfunc,model1,model2, rp.inverse(), 1,debugg);
            }
        }
    }

    delete nfunc;
    delete dfunc;
    return occlusionScores;
/*    std::vector<double> dvec;
    std::vector<double> nvec;
    DistanceWeightFunction2 * dfunc;
    DistanceWeightFunction2 * nfunc;

    unsigned int nr_models = models.size();
    vector<vector < OcclusionScore > > occlusionScores;
    occlusionScores.resize(nr_models);
    for(unsigned int i = 0; i < nr_models; i++){
        occlusionScores[i].resize(nr_models);
        for(unsigned int j = 0; j < nr_models; j++){
            occlusionScores[i][j] = OcclusionScore(0,0);
        }
    }

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

    dfunc = new DistanceWeightFunction2();
    dfunc->f = THRESHOLD;
    dfunc->p = 0.005;

    nfunc = new DistanceWeightFunction2();
    nfunc->f = THRESHOLD;
    nfunc->p = 0.50;

    double noiseWeight = dfunc->getNoise();

    for(unsigned int i = 0; i < models.size(); i++){
        Model * model1 = models[i];

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
                Eigen::Matrix4d p = model2->relativeposes[k].inverse()*(rps[j].inverse() * rps[i]);
                testgetDynamicWeights(false,dvec,nvec,dfunc,nfunc,p, model1->points,current_overlaps, current_occlusions, current_notocclusions,model2->frames[k],debugg);
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
    }

    return occlusionScores;
    */
}

/*
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
            Model * model2 = models[i];
            if(i != j){
                Eigen::Matrix4d rp = rps[j].inverse() * rps[i];
                OcclusionScore os = computeOcclusionScore(model1,model2, rp, 1,debugg);
            }
        }
	}

	delete nfunc;
	delete dfunc;
	return occlusionScores;
}
*/
/*
OcclusionScore ModelUpdater2::computeOcclusionScore(vector<superpoint> & spvec, Matrix4d cp, RGBDFrame* cf, ModelMask* cm, int step,  bool debugg){
    OcclusionScore oc;
printf("%s :: %i\n",__PRETTY_FUNCTION__,__LINE__);
    unsigned char  * dst_maskdata		= (unsigned char	*)(cm->mask.data);
    unsigned char  * dst_rgbdata		= (unsigned char	*)(cf->rgb.data);
    unsigned short * dst_depthdata		= (unsigned short	*)(cf->depth.data);
    float		   * dst_normalsdata	= (float			*)(cf->normals.data);

    float m00 = cp(0,0); float m01 = cp(0,1); float m02 = cp(0,2); float m03 = cp(0,3);
    float m10 = cp(1,0); float m11 = cp(1,1); float m12 = cp(1,2); float m13 = cp(1,3);
    float m20 = cp(2,0); float m21 = cp(2,1); float m22 = cp(2,2); float m23 = cp(2,3);

    Camera * dst_camera				= cf->camera;
    const unsigned int dst_width	= dst_camera->width;
    const unsigned int dst_height	= dst_camera->height;
    const float dst_idepth			= dst_camera->idepth_scale;
    const float dst_cx				= dst_camera->cx;
    const float dst_cy				= dst_camera->cy;
    const float dst_fx				= dst_camera->fx;
    const float dst_fy				= dst_camera->fy;
    const float dst_ifx				= 1.0/dst_camera->fx;
    const float dst_ify				= 1.0/dst_camera->fy;
    const unsigned int dst_width2	= dst_camera->width  - 2;
    const unsigned int dst_height2	= dst_camera->height - 2;

    int nr_data = spvec.size();

    std::vector<float> residuals;
    std::vector<int> debugg_src_inds;
    std::vector<int> debugg_dst_inds;
    std::vector<float> weights;
    residuals.reserve(nr_data);
    if(debugg){
        debugg_src_inds.reserve(nr_data);
        debugg_dst_inds.reserve(nr_data);
    }
    weights.reserve(nr_data);

    for(unsigned int ind = 0; ind < spvec.size();ind+=step){
        superpoint & sp = spvec[ind];

        float src_x = sp.x;
        float src_y = sp.y;
        float src_z = sp.z;
        float tz	= m20*src_x + m21*src_y + m22*src_z + m23;

        if(tz < 0){continue;}

        float src_nx = sp.nx;
        float src_ny = sp.ny;
        float src_nz = sp.nz;

        float point_information = sp.point_information;


        float tx	= m00*src_x + m01*src_y + m02*src_z + m03;
        float ty	= m10*src_x + m11*src_y + m12*src_z + m13;

        float itz	= 1.0/tz;
        float dst_w	= dst_fx*tx*itz + dst_cx;
        float dst_h	= dst_fy*ty*itz + dst_cy;

        if((dst_w > 0) && (dst_h > 0) && (dst_w < dst_width2) && (dst_h < dst_height2)){
            unsigned int dst_ind = unsigned(dst_h+0.5) * dst_width + unsigned(dst_w+0.5);

            float dst_z = dst_idepth*float(dst_depthdata[dst_ind]);
            float dst_nx = dst_normalsdata[3*dst_ind+0];
            if(dst_z > 0 && dst_nx != 2){
                //if(dst_detdata[dst_ind] != 0){continue;}
                float dst_ny = dst_normalsdata[3*dst_ind+1];
                float dst_nz = dst_normalsdata[3*dst_ind+2];

                float dst_x = (float(dst_w) - dst_cx) * dst_z * dst_ifx;
                float dst_y = (float(dst_h) - dst_cy) * dst_z * dst_ify;

                float tnx	= m00*src_nx + m01*src_ny + m02*src_nz;
                float tny	= m10*src_nx + m11*src_ny + m12*src_nz;
                float tnz	= m20*src_nx + m21*src_ny + m22*src_nz;

                double d = mysign(dst_z-tz)*fabs(tnx*(dst_x-tx) + tny*(dst_y-ty) + tnz*(dst_z-tz));
                double dst_noise = dst_z * dst_z;
                double point_noise = 1.0/sqrt(point_information);

                double compare_mul = sqrt(dst_noise*dst_noise + point_noise*point_noise);
                d *= compare_mul;

                double dist_dst = sqrt(dst_x*dst_x+dst_y*dst_y+dst_z*dst_z);
                double angle_dst = fabs((dst_x*dst_nx+dst_y*dst_ny+dst_z*dst_nz)/dist_dst);

                residuals.push_back(d);
                weights.push_back(angle_dst*angle_dst*angle_dst);
                if(debugg){
                    debugg_src_inds.push_back(ind);
                    debugg_dst_inds.push_back(dst_ind);
                }
            }
        }
    }

    //	DistanceWeightFunction2PPR2 * func = new DistanceWeightFunction2PPR2();
    //	func->maxp			= 1.0;
    //	func->update_size	= true;
    //	func->zeromean      = true;
    //	func->startreg		= 0.0001;
    //	func->debugg_print	= debugg;
    //	func->bidir			= true;
    //	func->maxnoise      = pred;
    //	func->reset();

    DistanceWeightFunction2 * func = new DistanceWeightFunction2();
    func->f = THRESHOLD;
    func->p = 0.02;

    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(1,residuals.size());
    for(unsigned int i = 0; i < residuals.size(); i++){X(0,i) = residuals[i];}
    func->computeModel(X);

    Eigen::VectorXd  W = func->getProbs(X);

    delete func;

    for(unsigned int i = 0; i < residuals.size(); i++){
        float r = residuals[i];
        float weight = W(i);
        float ocl = 0;
        if(r > 0){ocl += 1-weight;}
        oc.score		+= weight*weights.at(i);
        oc.occlusions	+= ocl*weights.at(i);
    }


    if(debugg){
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr dcloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        dcloud->points.resize(dst_width*dst_height);
        for(unsigned int dst_w = 0; dst_w < dst_width; dst_w++){
            for(unsigned int dst_h = 0; dst_h < dst_height;dst_h++){
                unsigned int dst_ind = dst_h*dst_width+dst_w;
                float z = dst_idepth*float(dst_depthdata[dst_ind]);
                if(z > 0){
                    float x = (float(dst_w) - dst_cx) * z * dst_ifx;
                    float y = (float(dst_h) - dst_cy) * z * dst_ify;
                    dcloud->points[dst_ind].x = x;
                    dcloud->points[dst_ind].y = y;
                    dcloud->points[dst_ind].z = z;
                    dcloud->points[dst_ind].r = dst_rgbdata[3*dst_ind+2];
                    dcloud->points[dst_ind].g = dst_rgbdata[3*dst_ind+1];
                    dcloud->points[dst_ind].b = dst_rgbdata[3*dst_ind+0];
                    if(dst_maskdata[dst_ind] == 255){
                        dcloud->points[dst_ind].r = 0;
                        dcloud->points[dst_ind].g = 255;
                        dcloud->points[dst_ind].b = 0;
                    }
                }
            }
        }

        oc.print();
        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr scloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
        for(unsigned int i = 0; i < spvec.size(); i++){
            superpoint & sp = spvec[i];
            pcl::PointXYZRGBNormal p;
            p.x = sp.x;
            p.y = sp.y;
            p.z = sp.z;

            float tx	= m00*p.x + m01*p.y + m02*p.z + m03;
            float ty	= m10*p.x + m11*p.y + m12*p.z + m13;
            float tz	= m20*p.x + m21*p.y + m22*p.z + m23;

            p.x = tx;
            p.y = ty;
            p.z = tz;

            p.normal_x = sp.nx;
            p.normal_y = sp.ny;
            p.normal_z = sp.nz;
            p.r = 0;//sp.r;
            p.g = 0;//sp.g;
            p.b = 255;//sp.b;

            scloud->points.push_back(p);
        }

        viewer->removeAllPointClouds();
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (scloud, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(scloud), "scloud");
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (dcloud, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(dcloud), "dcloud");
        viewer->spin();

        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

        for(unsigned int i = 0; i < residuals.size(); i++){
            float r = residuals[i];
            float weight = W(i);
            float ocl = 0;
            if(r > 0){ocl += 1-weight;}
            if(debugg){
                unsigned int src_ind = debugg_src_inds[i];
                unsigned int dst_ind = debugg_dst_inds[i];
                if(ocl > 0.01 || weight > 0.01){
                    scloud->points[src_ind].r = 255.0*ocl*weights.at(i);
                    scloud->points[src_ind].g = 255.0*weight*weights.at(i);
                    scloud->points[src_ind].b = 0;

                    if(i % 300 == 0){
                        char buf [1024];
                        sprintf(buf,"line%i",i);
                        viewer->addLine<pcl::PointXYZRGBNormal> (scloud->points[src_ind], dcloud->points[dst_ind],buf);
                    }
                }else{
                    scloud->points[src_ind].x = 0;
                    scloud->points[src_ind].y = 0;
                    scloud->points[src_ind].z = 0;
                }
            }
        }

        viewer->removeAllPointClouds();
        viewer->addPointCloud<pcl::PointXYZRGBNormal> (scloud, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(scloud), "scloud");
        viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 7, "scloud");

        viewer->addPointCloud<pcl::PointXYZRGBNormal> (dcloud, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(dcloud), "dcloud");
        viewer->spin();
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

    }
    //printf("stop :: %s::%i\n",__FILE__,__LINE__);
    return oc;
}
*/
}
////void ModelUpdater::computeMovingDynamicStatic(std::vector<cv::Mat> & movemask, std::vector<cv::Mat> & dynmask, vector<Matrix4d> bgcp, vector<RGBDFrame*> bgcf, vector<Matrix4d> poses, vector<RGBDFrame*> frames, bool debugg, std::string savePath){
////	static int segment_run_counter = -1;
////	segment_run_counter++;

////	double computeMovingDynamicStatic_startTime = getTime();

////	SegmentationResults sr;

////	int tot_nr_pixels = 0;
////	std::vector<int> offsets;

////	for(unsigned int i = 0; i < frames.size(); i++){
////		offsets.push_back(tot_nr_pixels);
////		unsigned int nr_pixels = frames[i]->camera->width * frames[i]->camera->height;
////		tot_nr_pixels += nr_pixels;
////	}

////	std::vector<unsigned char> labels;
////	labels.resize(tot_nr_pixels);

////	//Graph...
////	std::vector< std::vector<int> > interframe_connectionId;
////	std::vector< std::vector<float> > interframe_connectionStrength;
////	interframe_connectionId.resize(tot_nr_pixels);
////	interframe_connectionStrength.resize(tot_nr_pixels);

////	int current_point         = 0;
////	float * priors            = new float[3*tot_nr_pixels];
////	float * prior_weights     = new float[2*tot_nr_pixels];
////	bool * valids             = new bool[tot_nr_pixels];

////	std::vector<double> dvec;
////	std::vector<double> nvec;
////	DistanceWeightFunction2 * dfunc = 0;
////	DistanceWeightFunction2 * nfunc = 0;

////	double startTime = getTime();

////	std::vector< std::vector<superpoint> > framesp_test;
////	std::vector< std::vector<superpoint> > framesp;
////	for(unsigned int i = 0; i < frames.size(); i++){
////		framesp_test.push_back(frames[i]->getSuperPoints(Eigen::Matrix4d::Identity(),10,false));
////		framesp.push_back(frames[i]->getSuperPoints());
////	}
////	printf("frames init time: %5.5fs\n",getTime()-startTime);

////	startTime = getTime();
////    std::vector< std::vector<superpoint> > bgsp_test;
////	std::vector< std::vector<superpoint> > bgsp;
////	for(unsigned int i = 0; i < bgcf.size(); i++){
////        bgsp_test.push_back(bgcf[i]->getSuperPoints(Eigen::Matrix4d::Identity(),10,false));
////		bgsp.push_back(bgcf[i]->getSuperPoints());
////	}
////	printf("bg init time:     %5.5fs\n",getTime()-startTime);

////    startTime = getTime();


////    for(unsigned int i = 0; i < frames.size(); i++){
////        std::vector<superpoint> & framesp1_test = framesp_test[i];
////        std::vector<superpoint> & framesp1		= framesp[i];
////        for(unsigned int j = 0; j < frames.size(); j++){
////            if(i == j){continue;}
////            Eigen::Matrix4d p = poses[i].inverse() * poses[j];
////            std::vector<superpoint> & framesp2 = framesp[j];
////            getDynamicWeights(true,dvec,nvec,dfunc,nfunc,p.inverse(),frames[i],framesp1_test,framesp1, 0, 0, 0, frames[j],framesp2,0,0,interframe_connectionId,interframe_connectionStrength,false);
////        }
////    }

////    for(unsigned int i = 0; i < bgcf.size(); i++){
////        std::vector<superpoint> & framesp1_test = bgsp_test[i];
////        std::vector<superpoint> & framesp1		= bgsp[i];
////        for(unsigned int j = 0; j < bgcf.size(); j++){
////            if(i == j){continue;}
////            Eigen::Matrix4d p = bgcp[i].inverse() * bgcp[j];
////            std::vector<superpoint> & framesp2 = bgsp[j];
////            getDynamicWeights(true,dvec,nvec,dfunc,nfunc,p.inverse(),bgcf[i],framesp1_test,framesp1, 0, 0, 0, bgcf[j],framesp2,0,0,interframe_connectionId,interframe_connectionStrength,false);
////        }
////    }

////    for(unsigned int i = 0; i < frames.size(); i++){
////        std::vector<superpoint> & framesp1_test = framesp_test[i];
////        std::vector<superpoint> & framesp1		= framesp[i];
////        for(unsigned int j = 0; j < bgcf.size(); j++){
////            Eigen::Matrix4d p = poses[i].inverse() * bgcp[j];
////            std::vector<superpoint> & framesp2 = bgsp[j];
////            //getDynamicWeights(true,dvec,nvec,dfunc,nfunc,p.inverse(),frames[i],framesp1_test,framesp1, 0, 0, 0, bgcf[j],framesp2,0,0,interframe_connectionId,interframe_connectionStrength,false);
////        }
////    }


////	double dstdval = 0;
////	for(unsigned int i = 0; i < dvec.size(); i++){dstdval += dvec[i]*dvec[i];}
////	dstdval = sqrt(dstdval/double(dvec.size()-1));

////    //GeneralizedGaussianDistribution * dggdnfunc	= new GeneralizedGaussianDistribution(true,true,false,true,true);
////    //dggdnfunc->nr_refineiters					= 4;

////	GeneralizedGaussianDistribution * dggdnfunc	= new GeneralizedGaussianDistribution(true,true,false);
////    dggdnfunc->nr_refineiters					= 4;
////    DistanceWeightFunction2PPR3 * dfuncTMP		= new DistanceWeightFunction2PPR3(dggdnfunc,0.1,1000);
////	dfunc = dfuncTMP;
////	dfuncTMP->startreg				= 0.000;
////	dfuncTMP->max_under_mean		= false;
////    dfuncTMP->debugg_print			= false;
////	dfuncTMP->bidir					= true;
////	dfuncTMP->zeromean				= false;
////	dfuncTMP->maxp					= 0.9999;
////    dfuncTMP->maxd					= 0.5;
////    dfuncTMP->histogram_size		= 1000;
////	dfuncTMP->fixed_histogram_size	= false;
////	dfuncTMP->startmaxd				= dfuncTMP->maxd;
////	dfuncTMP->starthistogram_size	= dfuncTMP->histogram_size;
////    dfuncTMP->blurval				= 0.5;
////    dfuncTMP->blur                  = 0.01;
////	dfuncTMP->maxnoise				= dstdval;
////	dfuncTMP->compute_infront		= true;
////	dfuncTMP->ggd					= true;
////	dfuncTMP->reset();

////	if(savePath.size() != 0){
////		dfuncTMP->savePath = std::string(savePath)+"/segment_"+std::to_string(segment_run_counter)+"_dfunc.m";
////	}

////	dfunc->computeModel(dvec);

////	GeneralizedGaussianDistribution * ggdnfunc	= new GeneralizedGaussianDistribution(true,true);
////	ggdnfunc->nr_refineiters					= 4;
////	DistanceWeightFunction2PPR3 * nfuncTMP		= new DistanceWeightFunction2PPR3(ggdnfunc);
////	nfunc = nfuncTMP;
////    nfuncTMP->startreg				= 0.0;
////    nfuncTMP->debugg_print			= false;
////	nfuncTMP->bidir					= false;
////	nfuncTMP->zeromean				= true;
////	nfuncTMP->maxp					= 0.9999;
////    nfuncTMP->maxd					= 2.0;
////    nfuncTMP->histogram_size		= 1000;
////    nfuncTMP->fixed_histogram_size	= false;
////	nfuncTMP->startmaxd				= nfuncTMP->maxd;
////	nfuncTMP->starthistogram_size	= nfuncTMP->histogram_size;
////    nfuncTMP->blurval				= 0.5;
////    nfuncTMP->blur                  = 0.01;
////	nfuncTMP->stdval2				= 1;
////	nfuncTMP->maxnoise				= 1;
////	nfuncTMP->ggd					= true;
////	nfuncTMP->reset();
////	nfunc->computeModel(nvec);

////	if(savePath.size() != 0){
////		nfuncTMP->savePath = std::string(savePath)+"/segment_"+std::to_string(segment_run_counter)+"_nfunc.m";
////	}

////    dfunc->saveSimple(1.0,-1.0, 10000 , "./dfunc.bin");
////    nfunc->saveSimple(2.0, 0.0, 10000 , "./nfunc.bin");



////	printf("training time:     %5.5fs\n",getTime()-startTime);

////	long frameConnections = 0;
////	std::vector< std::vector< std::vector<float> > > pixel_weights;

////	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud  (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

////	double total_priortime = 0;
////	double total_connectiontime = 0;
////	double total_alloctime = 0;
////	double total_dealloctime = 0;
////	double total_Dynw = 0;

////	double maxprob_same = 0.999999999999999999;

////	for(unsigned int i = 0; i < frames.size(); i++){
////		if(debugg != 0){printf("currently workin on frame %i\n",i);}
////		int offset = offsets[i];
////		RGBDFrame * frame = frames[i];
////		float		   * normalsdata	= (float			*)(frame->normals.data);
////		startTime = getTime();
////		std::vector< std::vector<float> > probs = frame->getImageProbs();

////		total_connectiontime += getTime()-startTime;
////		pixel_weights.push_back(probs);
////		unsigned short * depthdata	= (unsigned short	*)(frame->depth.data);

////		Camera * camera				= frame->camera;
////		const unsigned int width	= camera->width;
////		const unsigned int height	= camera->height;
////		const float idepth			= camera->idepth_scale;
////		const float cx				= camera->cx;
////		const float cy				= camera->cy;
////		const float ifx				= 1.0/camera->fx;
////		const float ify				= 1.0/camera->fy;

////		startTime = getTime();
////		unsigned int nr_pixels = width * height;
////		double * current_overlaps	= new double[nr_pixels];
////		double * current_occlusions		= new double[nr_pixels];
////		double * current_notocclusions		= new double[nr_pixels];
////		for(unsigned int j = 0; j < nr_pixels; j++){
////			current_overlaps[j] = 1.0;
////			current_occlusions[j] = 1.0;
////			current_notocclusions[j] = 1.0;
////		}

////		double * bg_overlaps			= new double[nr_pixels];
////		double * bg_occlusions			= new double[nr_pixels];
////		double * bg_notocclusions		= new double[nr_pixels];
////		for(unsigned int j = 0; j < nr_pixels; j++){
////			bg_overlaps[j]	= 1.0;
////			bg_occlusions[j] = 1.0;
////			bg_notocclusions[j] = 1.0;
////		}
////		total_alloctime += getTime()-startTime;


////		startTime = getTime();
//////		for(unsigned int j = 0; j < frames.size(); j++){
//////			if(i == j){continue;}
//////			Eigen::Matrix4d p = poses[i].inverse() * poses[j];
//////			getDynamicWeights(false,dvec,nvec,dfunc,nfunc,p.inverse(),frames[i],framesp_test[i],framesp[i], current_overlaps, current_occlusions, current_notocclusions, frames[j],framesp[j],offsets[i],offsets[j],interframe_connectionId,interframe_connectionStrength,false);
//////		}


//////		for(unsigned int j = 0; j < bgcf.size(); j++){
//////			Eigen::Matrix4d p = poses[i].inverse() * bgcp[j];
//////			getDynamicWeights(false,dvec,nvec,dfunc,nfunc,p.inverse(),frames[i],framesp_test[i],framesp[i], bg_overlaps, bg_occlusions, bg_notocclusions, bgcf[j],bgsp[j],-1,-1,interframe_connectionId,interframe_connectionStrength,false);
//////		}



////		std::vector< std::vector<double> > bg_distances;
////		std::vector< std::vector<superpoint> > bg_points;
////		bg_distances.resize(640*480);
////		bg_points.resize(640*480);
////		for(unsigned int j = 0; j < bgcf.size(); j++){
////			Eigen::Matrix4d p = poses[i].inverse() * bgcp[j];
////			addReprojections(p.inverse(),dfunc,nfunc,bg_distances,bg_points,frames[i],framesp_test[i],framesp[i], bgcf[j],bgsp[j],-1,-1,interframe_connectionId,interframe_connectionStrength,0);
////		}

////		std::vector< std::vector<double> > current_distances;
////		std::vector< std::vector<superpoint> > current_points;
////		current_distances.resize(640*480);
////		current_points.resize(640*480);
////		for(unsigned int j = 0; j < frames.size(); j++){
////			if(i == j){continue;}
////			Eigen::Matrix4d p = poses[i].inverse() * poses[j];
////			addReprojections(p.inverse(),dfunc,nfunc,current_distances,current_points,frames[i],framesp_test[i],framesp[i],frames[j],framesp[j],offsets[i],offsets[j],interframe_connectionId,interframe_connectionStrength,0);
////		}

////		//printf("start agg\n");
////		double agg_start = getTime();
////		aggregateProbs(dfunc,nfunc,framesp[i],bg_distances, bg_points,bg_overlaps, bg_occlusions, bg_notocclusions);
////		printf("stop agg: %5.5fs\n",getTime()-agg_start);
////		//printf("start agg\n");
////		agg_start = getTime();
////		aggregateProbs(dfunc,nfunc,framesp[i],current_distances, current_points,current_overlaps, current_occlusions, current_notocclusions);
////		printf("stop agg: %5.5fs\n",getTime()-agg_start);

////		for(unsigned int j = 0; j < nr_pixels; j++){
////			bg_occlusions[j]		= std::min(0.99999,1-bg_occlusions[j]);
////			bg_overlaps[j]			= std::min(0.99999,1-bg_overlaps[j]);

////			current_occlusions[j]	= std::min(0.99999,1-current_occlusions[j]);
////			current_overlaps[j]		= std::min(0.99999,1-current_overlaps[j]);
////		}

////		total_Dynw += getTime()-startTime;

////		startTime = getTime();
////		unsigned char * detdata = (unsigned char*)(frame->det_dilate.data);
////		for(unsigned int h = 0; h < height;h++){
////			for(unsigned int w = 0; w < width;w++){
////				int ind = h*width+w;

////				valids[offset+ind] = detdata[ind] == 0 && normalsdata[3*ind] != 2;
////				setupPriors(3,
////                        current_occlusions[ind],current_overlaps[ind],bg_occlusions[ind],bg_overlaps[ind],valids[offset+ind],
////						priors[3*(offset+ind)+0], priors[3*(offset+ind)+1],priors[3*(offset+ind)+2],
////						prior_weights[2*(offset+ind)+0], prior_weights[2*(offset+ind)+1]);

////				if(probs[0][ind] > 0.00000001){frameConnections++;}
////				if(probs[1][ind] > 0.00000001){frameConnections++;}

////				current_point++;
////			}
////        }

////		double start_inf = getTime();
////		gc::Graph<double,double,double> * g = new gc::Graph<double,double,double>(nr_pixels,2*nr_pixels);
////		for(unsigned long ind = 0; ind < nr_pixels;ind++){
////			g -> add_node();
////			double weightFG = prior_weights[2*(offset+ind)+0];
////			double weightBG = prior_weights[2*(offset+ind)+1];
////			g -> add_tweights( ind, weightFG, weightBG );
////		}


////		for(unsigned int w = 0; w < width;w++){
////			for(unsigned int h = 0; h < height;h++){
////				int ind = h*width+w;
////				if(w > 0 && probs[0][ind] > 0.00000001){
////					int other = ind-1;
////					double p_same = std::min(double(probs[0][ind]),maxprob_same);
////					double not_p_same = 1-p_same;
////					double weight = -log(not_p_same);
////					g -> add_edge( ind, other, weight, weight );
////				}

////				if(h > 0 && probs[1][ind] > 0.00000001){
////					int other = ind-width;
////					double p_same = std::min(double(probs[1][ind]),maxprob_same);
////					double not_p_same = 1-p_same;
////					double weight = -log(not_p_same);
////					g -> add_edge( ind, other, weight, weight );
////				}
////			}
////		}

////		g -> maxflow();
////		for(unsigned long ind = 0; ind < nr_pixels;ind++){labels[offset+ind] = g->what_segment(ind);}

////		if(debugg != 0){printf("local inference time: %10.10fs\n\n",getTime()-start_inf);}




////		//		cv::Mat detrgb;
////		//		detrgb.create(height,width,CV_8UC3);

////		//		for(unsigned long ind = 0; ind < nr_pixels;ind++){
////		//			int dd = detdata[ind] == 0;
////		//			detrgb.data[3*ind+0] = dd * frame->rgb.data[3*ind+0];
////		//			detrgb.data[3*ind+1] = dd * frame->rgb.data[3*ind+1];
////		//			detrgb.data[3*ind+2] = dd * frame->rgb.data[3*ind+2];
////		//		}

////		//		cv::namedWindow( "edgeimg"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "edgeimg",		edgeimg );
////		//		cv::namedWindow( "det_dilate"	, cv::WINDOW_AUTOSIZE );		cv::imshow( "det_dilate",	frame->det_dilate);
////		//		cv::namedWindow( "det_dilate2"	, cv::WINDOW_AUTOSIZE );		cv::imshow( "det_dilate2",	detrgb);
////		//		cv::namedWindow( "rgb"			, cv::WINDOW_AUTOSIZE );		cv::imshow( "rgb",			frame->rgb );
////		//		cv::namedWindow( "depth"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "depth",		frame->depth );
////		//		cv::namedWindow( "priors"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "priors",		priorsimg );
////		//		cv::namedWindow( "labelimg"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "labelimg",		labelimg );
////		//		cv::waitKey(0);



////		if(savePath.size() != 0){
////			cv::Mat edgeimg;
////			edgeimg.create(height,width,CV_8UC3);
////			unsigned char * edgedata = (unsigned char*)edgeimg.data;

////			for(int j = 0; j < width*height; j++){
////				edgedata[3*j+0] = 0;
////				edgedata[3*j+1] = 255.0*(1-probs[0][j]);
////				edgedata[3*j+2] = 255.0*(1-probs[1][j]);
////			}

////			cv::Mat labelimg;
////			labelimg.create(height,width,CV_8UC3);
////			unsigned char * labelimgdata = (unsigned char*)labelimg.data;
////			for(unsigned long ind = 0; ind < nr_pixels;ind++){
////				double label = g->what_segment(ind);
////				labelimgdata[3*ind+0] = 255*label;
////				labelimgdata[3*ind+1] = 255*label;
////				labelimgdata[3*ind+2] = 255*label;
////			}

////			cv::Mat priorsimg;
////			priorsimg.create(height,width,CV_8UC3);
////			unsigned char * priorsdata = (unsigned char*)priorsimg.data;
////			for(unsigned long ind = 0; ind < nr_pixels;ind++){
////				priorsdata[3*ind+0]			= 255.0*priors[3*(offset+ind)+2];
////				priorsdata[3*ind+1]			= 255.0*priors[3*(offset+ind)+1];
////				priorsdata[3*ind+2]			= 255.0*priors[3*(offset+ind)+0];
////			}

////			cv::Mat current_overlapsimg;
////			current_overlapsimg.create(height,width,CV_8UC3);
////			unsigned char * current_overlapsdata = (unsigned char*)current_overlapsimg.data;
////			for(unsigned long ind = 0; ind < nr_pixels;ind++){
////				current_overlapsdata[3*ind+0]			= 255.0*current_overlaps[ind];
////				current_overlapsdata[3*ind+1]			= 255.0*current_overlaps[ind];
////				current_overlapsdata[3*ind+2]			= 255.0*current_overlaps[ind];
////			}

////			cv::Mat bg_overlapsimg;
////			bg_overlapsimg.create(height,width,CV_8UC3);
////			unsigned char * bg_overlapsdata = (unsigned char*)bg_overlapsimg.data;
////			for(unsigned long ind = 0; ind < nr_pixels;ind++){
////				bg_overlapsdata[3*ind+0]			= 255.0*bg_overlaps[ind];
////				bg_overlapsdata[3*ind+1]			= 255.0*bg_overlaps[ind];
////				bg_overlapsdata[3*ind+2]			= 255.0*bg_overlaps[ind];
////			}

////			cv::Mat current_occlusionsimg;
////			current_occlusionsimg.create(height,width,CV_8UC3);
////			unsigned char * current_occlusionsdata = (unsigned char*)current_occlusionsimg.data;
////			for(unsigned long ind = 0; ind < nr_pixels;ind++){
////				current_occlusionsdata[3*ind+0]			= 255.0*current_occlusions[ind];
////				current_occlusionsdata[3*ind+1]			= 255.0*current_occlusions[ind];
////				current_occlusionsdata[3*ind+2]			= 255.0*current_occlusions[ind];
////			}

////			cv::Mat bg_occlusionsimg;
////			bg_occlusionsimg.create(height,width,CV_8UC3);
////			unsigned char * bg_occlusionsdata = (unsigned char*)bg_occlusionsimg.data;
////			for(unsigned long ind = 0; ind < nr_pixels;ind++){
////				bg_occlusionsdata[3*ind+0]			= 255.0*bg_occlusions[ind];
////				bg_occlusionsdata[3*ind+1]			= 255.0*bg_occlusions[ind];
////				bg_occlusionsdata[3*ind+2]			= 255.0*bg_occlusions[ind];
////			}
////			if(debugg){
////				cv::imshow( "current_overlapsimg.png", current_overlapsimg );
////				cv::imshow( "bg_overlapsimg.png", bg_overlapsimg );
////				cv::imshow( "current_occlusionsimg.png", current_occlusionsimg );
////				cv::imshow( "bg_occlusionsimg.png", bg_occlusionsimg );


////				cv::namedWindow( "edgeimg"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "edgeimg",		edgeimg );
////				cv::namedWindow( "rgb"			, cv::WINDOW_AUTOSIZE );		cv::imshow( "rgb",			frame->rgb );
////				cv::namedWindow( "depth"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "depth",		frame->depth );
////				cv::namedWindow( "priors"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "priors",		priorsimg );
////				cv::namedWindow( "labelimg"		, cv::WINDOW_AUTOSIZE );		cv::imshow( "labelimg",		labelimg );
////				cv::waitKey(0);
////			}
////			if(savePath.size() != 0){
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_current_overlapsimg.png", current_overlapsimg );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_bg_overlapsimg.png", bg_overlapsimg );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_current_occlusionsimg.png", current_occlusionsimg );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_bg_occlusionsimg.png", bg_occlusionsimg );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_edgeimg.png", edgeimg );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_rgb.png", frame->rgb  );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_priors.png", priorsimg );
////				cv::imwrite( savePath+"/segment_"+std::to_string(segment_run_counter)+"_frame_"+std::to_string(i)+"_labelimg.png", labelimg );
////			}
////		}
////		delete g;

////		Eigen::Matrix4d p = poses[i];
////		float m00 = p(0,0); float m01 = p(0,1); float m02 = p(0,2); float m03 = p(0,3);
////		float m10 = p(1,0); float m11 = p(1,1); float m12 = p(1,2); float m13 = p(1,3);
////		float m20 = p(2,0); float m21 = p(2,1); float m22 = p(2,2); float m23 = p(2,3);
////		for(unsigned int h = 0; h < height;h++){
////			for(unsigned int w = 0; w < width;w++){
////				int ind = h*width+w;
////				float z = idepth*float(depthdata[ind]);
////				float x = (float(w) - cx) * z * ifx;
////				float y = (float(h) - cy) * z * ify;
////				pcl::PointXYZRGBNormal point;
////				point.x = m00*x + m01*y + m02*z + m03;
////				point.y = m10*x + m11*y + m12*z + m13;
////				point.z = m20*x + m21*y + m22*z + m23;
////				point.r = frame->rgb.data[3*ind+2];//priors[3*(offset+ind)+0]*255.0;
////				point.g = frame->rgb.data[3*ind+1];//priors[3*(offset+ind)+1]*255.0;
////				point.b = frame->rgb.data[3*ind+0];//priors[3*(offset+ind)+2]*255.0;
////				cloud->points.push_back(point);
////			}
////		}
////		total_priortime += getTime()-startTime;

////		startTime = getTime();
////		delete[] current_occlusions;
////		delete[] current_notocclusions;
////		delete[] current_overlaps;
////		delete[] bg_occlusions;
////		delete[] bg_notocclusions;
////		delete[] bg_overlaps;
////		total_dealloctime += getTime()-startTime;
////	}

////	delete dfuncTMP;
////	delete nfuncTMP;

////	printf("total_priortime        = %5.5fs\n",		total_priortime);
////	printf("total_connectiontime   = %5.5fs\n",		total_connectiontime);
////	printf("total_alloctime        = %5.5fs\n",		total_alloctime);
////	printf("total_dealloctime      = %5.5fs\n",		total_dealloctime);
////	printf("total_Dynw             = %5.5fs\n",		total_Dynw);

////	long interframeConnections = 0;
////	for(unsigned int i = 0; i < interframe_connectionId.size();i++){interframeConnections += interframe_connectionId[i].size();}

////	double start_inf = getTime();
////	gc::Graph<double,double,double> * g = new gc::Graph<double,double,double>(current_point,frameConnections+interframeConnections);
////	for(unsigned long i = 0; i < current_point;i++){
////		g -> add_node();
////		double weightFG = prior_weights[2*i+0];
////		double weightBG = prior_weights[2*i+1];

////		g -> add_tweights( i, weightFG, weightBG );
////	}

////	//double maxprob_same = 0.99999999999;
////	for(unsigned int i = 0; i < frames.size(); i++){
////		int offset = offsets[i];
////		Camera * camera				= frames[i]->camera;
////		const unsigned int width	= camera->width;
////		const unsigned int height	= camera->height;
////		std::vector< std::vector<float> > & probs = pixel_weights[i];
////		for(unsigned int w = 0; w < width;w++){
////			for(unsigned int h = 0; h < height;h++){
////				int ind = h*width+w;
////				if(w > 0 && probs[0][ind] > 0.00000001 && w < width-1){
////					int other = ind-1;
////					double p_same = std::min(double(probs[0][ind]),maxprob_same);
////					double not_p_same = 1-p_same;
////					double weight = -log(not_p_same);
////					g -> add_edge( ind+offset, other+offset, weight, weight );
////				}

////				if(h > 0 && probs[1][ind] > 0.00000001 && h < height-1) {
////					int other = ind-width;
////					double p_same = std::min(double(probs[1][ind]),maxprob_same);
////					double not_p_same = 1-p_same;
////					double weight = -log(not_p_same);
////					g -> add_edge( ind+offset, other+offset, weight, weight );
////				}
////			}
////		}
////	}

////	std::vector<std::vector<unsigned long> > interframe_connectionId_added;
////	interframe_connectionId_added.resize(interframe_connectionId.size());

////	std::vector<std::vector<double> > interframe_connectionStrength_added;
////	interframe_connectionStrength_added.resize(interframe_connectionStrength.size());

////	double initAdded = 0;
////	for(unsigned int i = 0; i < interframe_connectionId.size();i++){
////		for(unsigned int j = 0; j < interframe_connectionId[i].size();j++){
////			double weight = interframe_connectionStrength[i][j];
////			unsigned long other = interframe_connectionId[i][j];
////			if(weight > 0.01 && labels[i] != labels[other]){
////				g -> add_edge( i, other, weight, weight );

////				interframe_connectionId_added[i].push_back(other);
////				interframe_connectionStrength_added[i].push_back(weight);

////				interframe_connectionStrength[i][j] = interframe_connectionStrength[i].back();
////				interframe_connectionStrength[i].pop_back();

////				interframe_connectionId[i][j] = interframe_connectionId[i].back();
////				interframe_connectionId[i].pop_back();
////				j--;

////				initAdded++;
////			}
////		}
////	}

////	g -> maxflow();
////	for(unsigned long ind = 0; ind < current_point;ind++){labels[ind] = g->what_segment(ind);}

////	double tot_inf = 0;
////	for(unsigned int it = 0; it < 40; it++){
////		double start_inf1 = getTime();

////		double diffs = 0;
////		for(unsigned int i = 0; i < interframe_connectionId.size();i++){
////			for(unsigned int j = 0; j < interframe_connectionId[i].size();j++){
////				double weight = interframe_connectionStrength[i][j];
////				unsigned long other = interframe_connectionId[i][j];
////				if(weight > 0.01 && labels[i] != labels[other]){diffs++;}
////			}
////		}


////		double adds = 100000;
////		double prob = std::min(adds / diffs,1.0);
////		printf("diffs: %f adds: %f prob: %f ",diffs,adds,prob);
////		printf("ratio of total diffs: %f\n",diffs/double(frameConnections+interframeConnections));

////		if(diffs == 0){break;}

////		double trueadds = 0;
////		for(unsigned int i = 0; i < interframe_connectionId.size();i++){
////			for(unsigned int j = 0; j < interframe_connectionId[i].size();j++){
////				double weight = interframe_connectionStrength[i][j];
////				unsigned long other = interframe_connectionId[i][j];
////				if(weight > 0.01 && labels[i] != labels[other]){
////					if(rand() <= prob*RAND_MAX){
////						trueadds++;
////						g -> add_edge( i, other, weight, weight );

////						interframe_connectionId_added[i].push_back(other);
////						interframe_connectionStrength_added[i].push_back(weight);

////						interframe_connectionStrength[i][j] = interframe_connectionStrength[i].back();
////						interframe_connectionStrength[i].pop_back();

////						interframe_connectionId[i][j] = interframe_connectionId[i].back();
////						interframe_connectionId[i].pop_back();
////						j--;
////					}
////				}
////			}
////		}

////		g -> maxflow();
////		for(unsigned long ind = 0; ind < current_point;ind++){labels[ind] = g->what_segment(ind);}

////		tot_inf += getTime()-start_inf1;
////		if(debugg != 0){printf("static inference1 time: %10.10fs total: %10.10f\n\n",getTime()-start_inf1,tot_inf);}

////		if(tot_inf > 90){break;}
////	}

////	double interfrace_constraints_added = 0;
////	for(unsigned int i = 0; i < interframe_connectionId.size();i++){
////		interfrace_constraints_added += interframe_connectionId_added[i].size();
////		for(unsigned int j = 0; j < interframe_connectionId_added[i].size();j++){
////			interframe_connectionStrength[i].push_back(interframe_connectionStrength_added[i][j]);
////			interframe_connectionId[i].push_back(interframe_connectionId_added[i][j]);
////		}
////	}

////	delete g;
////	double end_inf = getTime();
////	printf("static inference time: %10.10fs interfrace_constraints added ratio: %f\n",end_inf-start_inf,interfrace_constraints_added/double(interframeConnections));

////	const unsigned int nr_frames		= frames.size();
////	const unsigned int width			= frames[0]->camera->width;
////	const unsigned int height			= frames[0]->camera->height;
////	const unsigned int pixels_per_image	= width*height;
////	const unsigned int nr_pixels		= nr_frames*pixels_per_image;

////	double probthresh = 0.5;
////	double str_probthresh = -log(probthresh);
////	unsigned int number_of_dynamics = 0;
////	unsigned int nr_obj_dyn = 0;
////	unsigned int nr_obj_mov = 0;
////	std::vector<unsigned int> objectlabel;
////	std::vector<int> labelID;
////	labelID.push_back(0);
////	objectlabel.resize(nr_pixels);

////	for(unsigned long i = 0; i < nr_pixels; i++){objectlabel[i] = 0;}
////	for(unsigned long ind = 0; ind < nr_pixels; ind++){
////		if(valids[ind] && objectlabel[ind] == 0 && labels[ind] != 0){
////			unsigned int current_label = labels[ind];
////			number_of_dynamics++;
////			objectlabel[ind] = number_of_dynamics;
////			unsigned long todocounter = 0;
////			std::vector< unsigned long > todo;
////			todo.push_back(ind);

////            double score0 = 0;//dynamic
////            double score1 = 0;//Moving

////			double pscore0 = 0;
////			double pscore1 = 0;

////			double nscore0 = 0;
////			double nscore1 = 0;

////			double totsum = 0;
////			while(todocounter < todo.size()){
////				unsigned long cind = todo[todocounter++];
////				unsigned long frameind = cind / pixels_per_image;


////				unsigned long iind = cind % pixels_per_image;
////				unsigned long w = iind % width;
////				unsigned long h = iind / width;

////				double p0 = priors[3*cind+0];
////				double p1 = priors[3*cind+1];
////				double p2 = priors[3*cind+2];

////				if(valids[cind]){
////                    if(p0 > p1){score0 += p0 - p1;}//Probably moving
////                    else{       score1 += p1 - p0;}//Probably dynamic

//                    totsum++;
//                }

//                float * dedata = (float*)(frames[frameind]->de.data);
//                unsigned short * depthdata = (unsigned short *)(frames[frameind]->depth.data);
//                if(depthdata[iind] == 0){printf("big giant WTF... file %i line %i\n",__FILE__,__LINE__);continue;}

//                int dir;
//                dir = -1;
//                if( w > 0 && objectlabel[cind+dir] == 0 && labels[cind+dir] == current_label && depthdata[iind+dir] != 0 && (dedata[3*(iind+dir)+1]+dedata[3*(iind+dir)+2]) < probthresh){
//                    objectlabel[cind+dir] = number_of_dynamics;
//                    todo.push_back(cind+dir);
//                }

//                dir = 1;
//                if( w < (width-1) && objectlabel[cind+dir] == 0 && labels[cind+dir] == current_label && depthdata[iind+dir] != 0 && (dedata[3*(iind+dir)+1]+dedata[3*(iind+dir)+2]) < probthresh){
//                    objectlabel[cind+dir] = number_of_dynamics;
//                    todo.push_back(cind+dir);
//                }

//                dir = -width;
//                if( h > 0 && objectlabel[cind+dir] == 0 && labels[cind+dir] == current_label && depthdata[iind+dir] != 0 && (dedata[3*(iind+dir)+1]+dedata[3*(iind+dir)+2]) < probthresh){
//                    objectlabel[cind+dir] = number_of_dynamics;
//                    todo.push_back(cind+dir);
//                }

//                dir = width;
//                if( h < (height-1) && objectlabel[cind+dir] == 0 && labels[cind+dir] == current_label && depthdata[iind+dir] != 0 && (dedata[3*(iind+dir)+1]+dedata[3*(iind+dir)+2]) < probthresh){
//                    objectlabel[cind+dir] = number_of_dynamics;
//                    todo.push_back(cind+dir);
//                }

//                for(unsigned long j = 0; j < interframe_connectionId[cind].size();j++){
//                    unsigned long other = interframe_connectionId[cind][j];
//                    if(interframe_connectionStrength[cind][j] > str_probthresh && objectlabel[other] == 0 && labels[other] == current_label){
//                        objectlabel[other] = number_of_dynamics;
//                        todo.push_back(other);
//                    }
//                }
//            }

////			score0 = pscore0+nscore0;
////			score1 = pscore1+nscore1;


//            labelID.push_back(0);
//            if(debugg != 0){
//                if(totsum > 100){
//                    printf("---------------------------\n");
//                    printf("score0: %10.10f score1: %10.10f ",score0,score1);
//                    printf("totsum: %10.10f\n",totsum);

//                    printf("pscore0: %10.10f nscore0: %10.10f ",pscore0,nscore0);
//                    printf("pscore1: %10.10f nscore1: %10.10f\n",pscore1,nscore1);
//                }
//            }

//            if(std::max(score0,score1) < 100){continue;}

//            score0 *= 2.0;

//            if(score1 > score0){
//                labelID.back() = ++nr_obj_dyn;
//                if(debugg != 0){printf("Dynamic: %f -> %f\n",score1,totsum);}
//                sr.component_dynamic.push_back(todo);
//                sr.scores_dynamic.push_back(score1);
//                sr.total_dynamic.push_back(totsum);
//            }else{
//                labelID.back() = --nr_obj_mov;
//                if(debugg != 0){printf("Moving: %f -> %f\n",score0,totsum);}
//                sr.component_moving.push_back(todo);
//                sr.scores_moving.push_back(score0);
//                sr.total_moving.push_back(totsum);
//            }
//        }
//    }

//    for(unsigned long ind = 0; ind < nr_pixels; ind++){
//        unsigned int ol = objectlabel[ind];
//        if(ol != 0 && labelID[ol] == 0){
//            objectlabel[ind] = 0;
//            labels[ind] = 0;
//        }
//    }

//    for(unsigned int i = 0; i < interframe_connectionId.size();i++){interframe_connectionId[i].clear();}
//    interframe_connectionId.clear();

//    for(unsigned int i = 0; i < interframe_connectionStrength.size();i++){interframe_connectionStrength[i].clear();}
//    interframe_connectionStrength.clear();

//    printf("connectedComponent: %5.5fs\n",getTime()-start_inf);

//    int current = 0;
//    for(unsigned long i = 0; i < frames.size(); i++){
//        Camera * camera				= frames[i]->camera;
//        const unsigned int width	= camera->width;
//        const unsigned int height	= camera->height;

//        cv::Mat m;
//        m.create(height,width,CV_8UC1);
//        unsigned char * mdata = (unsigned char*)m.data;

//        cv::Mat d;
//        d.create(height,width,CV_8UC1);
//        unsigned char * ddata = (unsigned char*)d.data;

//        cv::Mat d2;
//        d2.create(height,width,CV_8UC1);
//        unsigned char * ddata2 = (unsigned char*)d2.data;


//        for(int j = 0; j < width*height; j++){
//            mdata[j] = 0;
//            ddata[j] = 0;
//            ddata2[j] = labels[current];
//            unsigned int label = objectlabel[current];
//            int lid = labelID[label];
//            if(lid >  0){
//                ddata[j] = lid;
//            }else if(lid < 0){
//                mdata[j] = -lid;
//            }
//            current++;
//        }
//        movemask.push_back(m);
//        dynmask.push_back(d);
//    }

//    if(savePath.size() != 0){
//        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_sample (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
//        cloud_sample->points.resize(current_point);
//        cloud_sample->width = current_point;
//        cloud_sample->height = 1;
//        for(unsigned int i = 0; i < current_point; i++){
//            cloud_sample->points[i]	= cloud->points[i];
//            cloud_sample->points[i].r = priors[3*i+0]*255.0;
//            cloud_sample->points[i].g = priors[3*i+1]*255.0;
//            cloud_sample->points[i].b = priors[3*i+2]*255.0;
//        }
//        if(cloud_sample->points.size()){
//            pcl::io::savePCDFileBinaryCompressed (savePath+"/segment_"+std::to_string(segment_run_counter)+"_priors.pcd", *cloud_sample);
//        }

//        for(unsigned int i = 0; i < current_point; i++){
//            cloud_sample->points[i].r = 0;
//            cloud_sample->points[i].g = 0;
//            cloud_sample->points[i].b = 255;
//        }

//        for(unsigned int c = 0; c < sr.component_dynamic.size(); c++){
//            for(unsigned int i = 0; i < sr.component_dynamic[c].size(); i++){
//                cloud_sample->points[sr.component_dynamic[c][i]].r = 0;
//                cloud_sample->points[sr.component_dynamic[c][i]].g = 255;
//                cloud_sample->points[sr.component_dynamic[c][i]].b = 0;
//            }
//        }

//        for(unsigned int c = 0; c < sr.component_moving.size(); c++){
//            for(unsigned int i = 0; i < sr.component_moving[c].size(); i++){
//                cloud_sample->points[sr.component_moving[c][i]].r = 255;
//                cloud_sample->points[sr.component_moving[c][i]].g = 0;
//                cloud_sample->points[sr.component_moving[c][i]].b = 0;
//            }
//        }

//        if(cloud_sample->points.size()){
//            pcl::io::savePCDFileBinaryCompressed (savePath+"/segment_"+std::to_string(segment_run_counter)+"_classes.pcd", *cloud_sample);
//        }

//        for(unsigned int i = 0; i < current_point; i++){
//            cloud_sample->points[i].r = 0;
//            cloud_sample->points[i].g = 0;
//            cloud_sample->points[i].b = 255;
//        }

//        for(unsigned int c = 0; c < sr.component_dynamic.size(); c++){
//            int randr = rand()%256;
//            int randg = rand()%256;
//            int randb = rand()%256;

//            for(unsigned int i = 0; i < sr.component_dynamic[c].size(); i++){
//                cloud_sample->points[sr.component_dynamic[c][i]].r = randr;
//                cloud_sample->points[sr.component_dynamic[c][i]].g = randg;
//                cloud_sample->points[sr.component_dynamic[c][i]].b = randb;
//            }
//        }

//        for(unsigned int c = 0; c < sr.component_moving.size(); c++){
//            int randr = rand()%256;
//            int randg = rand()%256;
//            int randb = rand()%256;

//            for(unsigned int i = 0; i < sr.component_moving[c].size(); i++){
//                cloud_sample->points[sr.component_moving[c][i]].r = randr;
//                cloud_sample->points[sr.component_moving[c][i]].g = randg;
//                cloud_sample->points[sr.component_moving[c][i]].b = randb;
//            }
//        }

//        if(cloud_sample->points.size()){
//            pcl::io::savePCDFileBinaryCompressed (savePath+"/segment_"+std::to_string(segment_run_counter)+"_clusters.pcd", *cloud_sample);
//        }

//        cloud->width = cloud->points.size();
//        cloud->height = 1;

//        if(cloud->points.size()){
//            pcl::io::savePCDFileBinaryCompressed (savePath+"/segment_"+std::to_string(segment_run_counter)+"_full.pcd", *cloud);
//        }

//        cloud_sample->points.resize(0);
//        for(unsigned int c = 0; c < sr.component_dynamic.size(); c++){
//            for(unsigned int i = 0; i < sr.component_dynamic[c].size(); i++){
//                cloud_sample->points.push_back(cloud->points[sr.component_dynamic[c][i]]);
//            }
//        }
//        cloud_sample->width = cloud_sample->points.size();
//        cloud_sample->height = 1;

//        if(cloud_sample->points.size()){
//            pcl::io::savePCDFileBinaryCompressed (savePath+"/segment_"+std::to_string(segment_run_counter)+"_dynamicobjects.pcd", *cloud_sample);
//        }
//    }


//    printf("computeMovingDynamicStatic total time: %5.5fs\n",getTime()-computeMovingDynamicStatic_startTime);
//    if(debugg){
//        pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_sample (new pcl::PointCloud<pcl::PointXYZRGBNormal>);

//        for(unsigned int i = 0; i < current_point; i++){
//            cloud->points[i].r = priors[3*i+0]*255.0;
//            cloud->points[i].g = priors[3*i+1]*255.0;
//            cloud->points[i].b = priors[3*i+2]*255.0;
//        }

//        cloud_sample->points.clear();
//        for(unsigned int i = 0; i < current_point; i++){
//            if(rand() % 4 == 0){
//                cloud_sample->points.push_back(cloud->points[i]);
//            }
//        }
//        viewer->removeAllPointClouds();
//        viewer->addPointCloud<pcl::PointXYZRGBNormal> (cloud_sample, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(cloud_sample), "cloud");
//        viewer->spin();

//        cloud_sample->points.clear();
//        for(unsigned int i = 0; i < current_point; i++){
//            if(rand() % 4 == 0){

//                double p_fg = exp(-prior_weights[2*i+0]);

//                cloud_sample->points.push_back(cloud->points[i]);
//                cloud_sample->points.back().r = p_fg*255.0;//(priors[3*i+0]+priors[3*i+2])*255.0;
//                cloud_sample->points.back().g = p_fg*255.0;
//                cloud_sample->points.back().b = p_fg*255.0;//;
//            }
//        }
//        viewer->removeAllPointClouds();
//        viewer->addPointCloud<pcl::PointXYZRGBNormal> (cloud_sample, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(cloud_sample), "cloud");
//        viewer->spin();


//        cloud_sample->points.clear();
//        for(unsigned int i = 0; i < current_point; i++){
//            if(rand() % 4 == 0 && labels[i] == 0){
//                cloud_sample->points.push_back(cloud->points[i]);
//                cloud_sample->points.back().r = 0;
//                cloud_sample->points.back().g = 0;
//                cloud_sample->points.back().b = 255;
//            }
//        }

//        for(unsigned int c = 0; c < sr.component_dynamic.size(); c++){
//            int randr = rand()%156;
//            int randg = rand()%256;
//            int randb = rand()%256;
//            for(unsigned int i = 0; i < sr.component_dynamic[c].size(); i++){
//                cloud_sample->points.push_back(cloud->points[sr.component_dynamic[c][i]]);
//                cloud_sample->points.back().r = randr;
//                cloud_sample->points.back().g = randg;
//                cloud_sample->points.back().b = randb;
//            }
//        }

//        for(unsigned int c = 0; c < sr.component_moving.size(); c++){
//            int randr = rand()%256;
//            int randg = rand()%256;
//            int randb = rand()%256;
//            for(unsigned int i = 0; i < sr.component_moving[c].size(); i++){
//                cloud_sample->points.push_back(cloud->points[sr.component_moving[c][i]]);
//                cloud_sample->points.back().r = 255;//randr;
//                cloud_sample->points.back().g = 0;//randg;
//                cloud_sample->points.back().b = 0;//randb;
//            }
//        }
//        viewer->removeAllPointClouds();
//        viewer->addPointCloud<pcl::PointXYZRGBNormal> (cloud_sample, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(cloud_sample), "cloud");
//        viewer->spin();

//        cloud_sample->points.clear();
//        for(unsigned int i = 0; i < current_point; i++){
//            if(rand() % 1 == 0){
//                cloud_sample->points.push_back(cloud->points[i]);
//            }
//        }
//        viewer->removeAllPointClouds();
//        viewer->addPointCloud<pcl::PointXYZRGBNormal> (cloud_sample, pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal>(cloud_sample), "cloud");
//        viewer->spin();
//    }

//    delete[] valids;
//    delete[] priors;
//    delete[] prior_weights;
//    //printf("computeMovingDynamicStatic total time: %5.5fs\n",getTime()-computeMovingDynamicStatic_startTime);
//}



