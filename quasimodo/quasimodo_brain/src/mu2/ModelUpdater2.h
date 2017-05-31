#ifndef reglibModelUpdater2_H
#define reglibModelUpdater2_H

#include "../../../quasimodo_models/include/modelupdater/ModelUpdater.h"

#include "registration/nanoflann.hpp"

namespace reglib
{
	class ModelUpdater2: public ModelUpdater{
		public:

		std::vector<superpoint> fusedmodel;
		Registration * registration;

		//ModelGraph * graph;

		ModelUpdater2(Registration * registration_ = new Registration());
		ModelUpdater2(Model * model_,Registration * registration_ = new Registration());
		~ModelUpdater2();

		virtual FusionResults registerModel(Model * model2, Eigen::Matrix4d guess = Eigen::Matrix4d::Identity(), double uncertanity = -1);
		//virtual void addFrame(RGBDFrame * frame, Eigen::Matrix4d pose);
		virtual void setRegistration( Registration * registration_);
		virtual void fuse(Model * model2, Eigen::Matrix4d guess = Eigen::Matrix4d::Identity(), double uncertanity = -1);
		virtual UpdatedModels fuseData(FusionResults * f, Model * model1, Model * model2);
		//virtual void refine();
		virtual	void computeMassRegistration(std::vector<Eigen::Matrix4d> current_poses, std::vector<RGBDFrame*> current_frames,std::vector<cv::Mat> current_masks);

		virtual OcclusionScore computeOcclusionScore2(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, vector<superpoint> & spvec, Matrix4d cp, RGBDFrame* cf, ModelMask* cm, int step = 1,  bool debugg = false);
		virtual OcclusionScore computeOcclusionScore2(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * mod, vector<Matrix4d> cp, vector<RGBDFrame*> cf, vector<ModelMask*> cm, Matrix4d rp = Matrix4d::Identity(), int step = 1, bool debugg = false);
		virtual OcclusionScore computeOcclusionScore2(DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, Model * model1, Model * model2, Matrix4d rp = Matrix4d::Identity(),int step = 1, bool debugg = false);
		virtual vector<vector < OcclusionScore > > computeOcclusionScore2(vector<Model *> models, vector<Matrix4d> rps, int step = 1, bool debugg = false);
		virtual std::vector<int> getPartition2(std::vector< std::vector< float > > & scores);
		virtual std::vector<int> partition_graph2(std::vector< std::vector< float > > & scores);
		virtual float recursive_split2(std::vector<Graph*> * graphs_out, std::vector<std::vector<int>> * graphinds_out, Graph * graph, std::vector<int> graph_inds);

	};

}

#endif // reglibModelUpdater2_H
