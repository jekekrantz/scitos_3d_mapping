#ifndef RegistrationRandom2_H
#define RegistrationRandom2_H

#include "../../../quasimodo_models/include/registration/Registration.h"
#include <time.h>

#include "../../../quasimodo_models/include/registration/RegistrationRefinement.h"

namespace reglib
{
    class RegistrationOptimization
    {
        public:
		double overlap;
		double info;

        Eigen::Matrix<double, 6, 1> ATb;
        Eigen::Matrix<double, 6, 6> ATA;

        RegistrationOptimization(Eigen::Matrix<double, 6, 1> ATb_ = Eigen::Matrix<double, 6, 1>::Zero(), Eigen::Matrix<double, 6, 6> ATA_ = Eigen::Matrix<double, 6, 6>::Identity()){
            ATb = ATb_;
            ATA = ATA_;
        }

        Eigen::Affine3d solve(){
            Vector6d x = static_cast< Eigen::Matrix<double, 6, 1> > (ATA.inverse () * ATb);
            return Eigen::Affine3d(constructTransformationMatrix(-x(0,0),-x(1,0),-x(2,0),-x(3,0),-x(4,0),-x(5,0)));
        }

		double getScore(){
			VectorXcd eivals = ATA.eigenvalues();
			double score = eivals(0).real();
			for(unsigned int d = 1; d < 6; d++){score = std::min(eivals(d).real(),score);}
			return score;
		}
    };

	class InternalFusionResults
	{
		public:
		RegistrationOptimization ro;
		FusionResults fr;
		InternalFusionResults(FusionResults fr_){
			fr = fr_;
		}
	};

	class RegistrationRandom2 : public Registration
	{
		public:

		Eigen::Matrix4d initTransform;

		unsigned int nr_src_target;
		unsigned int nr_dst_target;

		DistanceWeightFunction2 * dfunc_;
		DistanceWeightFunction2 * nfunc_;

		std::vector<superpoint> src_random;
		std::vector<superpoint> dst_random;

		std::vector<superpoint> Tsrc;
		std::vector<superpoint> Tdst;

		unsigned int steprx;
		unsigned int stepry;
		unsigned int steprz;
		double start_rx;
		double start_ry;
		double start_rz;
		double stop_rx;
		double stop_ry;
		double stop_rz;

		unsigned int steptx;
		unsigned int stepty;
		unsigned int steptz;
		double start_tx;
		double start_ty;
		double start_tz;
		double stop_tx;
		double stop_ty;
		double stop_tz;

		int dst_nr_arraypoints;
		double * dst_arraypoints;
		Tree3d * dst_trees3d;
		ArrayData3D<double> * dst_a3d;

		int src_nr_arraypoints;
		double * src_arraypoints;
		Tree3d * src_trees3d;
		ArrayData3D<double> * src_a3d;

		double rotationweight;

		RegistrationRandom2(unsigned int steps = 4);
		~RegistrationRandom2();


		virtual double getMeanDist(std::vector<superpoint> & data);
		virtual void deMean(std::vector<superpoint> & Tsrc);


		virtual void setSrc(std::vector<superpoint> & src_);
		virtual void setDst(std::vector<superpoint> & dst_);


        RegistrationOptimization getOptimization(std::vector<double> & standard_div, DistanceWeightFunction2 * dfunc, DistanceWeightFunction2 * nfunc, std::vector<reglib::superpoint> & A, std::vector<reglib::superpoint> & Q);

		void transformSuperPoints(std::vector<superpoint> & spvec, Eigen::Matrix4d cp);

		std::vector<reglib::superpoint> getMatches(const std::vector<reglib::superpoint> & X, const std::vector<reglib::superpoint> & Y, Tree3d * Y_tree);

		void pruneDuplicates(std::vector<InternalFusionResults> & fr_X, double threshold, double rotationWeight);

		double getStdDivVec(const std::vector<double> & V);
		std::vector<double> getStdDiv(std::vector<reglib::superpoint> & A, std::vector<reglib::superpoint> & Q);
		void getResiduals(std::vector<double> & residuals_d, std::vector<double> & residuals_a, std::vector<double> & standard_div, std::vector<reglib::superpoint> & A, std::vector<reglib::superpoint> & Q);
		std::vector<superpoint> getRandom(std::vector<superpoint> points);

		InternalFusionResults refine(InternalFusionResults guess, unsigned int nrp = 1000, bool regularize = true);
		FusionResults getTransform(Eigen::MatrixXd guess);

		Eigen::Matrix4d getTransform(std::vector<superpoint> & X, std::vector<superpoint> & Y);


		void show(std::vector<reglib::superpoint> & X1, std::vector<reglib::superpoint> & Y1,std::vector<reglib::superpoint> & X2, std::vector<reglib::superpoint> & Y2 );

	};
}

#endif
