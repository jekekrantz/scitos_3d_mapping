#ifndef DistanceWeightFunction2PPR3_H
#define DistanceWeightFunction2PPR3_H

#include <cmath>
#include <sys/time.h>
#include "DistanceWeightFunction2.h"

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres/iteration_callback.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>

#include "weightfunctions/Distribution.h"
#include "core/Util.h"
#include "weightfunctions/SignalProcessing.h"


using namespace Eigen;
namespace reglib{

class DistanceWeightFunction2PPR3 : public DistanceWeightFunction2
{
public:

    SignalProcessing * sp;
    Distribution * dist;

	bool useIRLSreweight;

    bool fixed_histogram_size;

    double stdval;
	double stdval2;
	double mulval;
	double meanval;
	double meanval2;

    double maxnoise;

    bool zeromean;

	double costpen;

	double maxp;

	bool first;

	bool update_size;
	double target_length;
	double data_per_bin;
	double meanoffset;
	double blur;

	double startmaxd;
	double maxd;
	double mind;
    double noise_min;

	int histogram_size;
	int starthistogram_size;
	double blurval;
	double stdgrow;

	double noiseval;
	double startreg;

	bool ggd;
	bool compute_infront;

//    std::vector<float> prob;
//    std::vector<float> infront;
//    std::vector<float> histogram;
//    std::vector<float> blur_histogram;
//    std::vector<float> noise;
//    std::vector<float> noisecdf;

	float * prob;
	float * irls;
    float * infront;
    float * histogram;
    float * blur_histogram;
    float * noise;
    float * noisecdf;

	int nr_refineiters;
	bool refine_mean;
	bool refine_mul;
	bool refine_std;

	bool threshold;
	bool uniform_bias;
	bool scale_convergence;
	double nr_inliers;

	bool rescaling;

	bool interp;

	bool max_under_mean;

	bool bidir;
	int iter;

	double histogram_mul;
	double histogram_mul2;

    int min_histogram_size;
    int max_histogram_size;

    double reg_shrinkage;

	double tinl;
	double toth;

    bool tune;


    virtual DistanceWeightFunction2 * clone();

    DistanceWeightFunction2PPR3(Distribution * dist,	double maxd_	= 0.25, int histogram_size_ = 25000);

    DistanceWeightFunction2PPR3(	double maxd_	= 0.25, int histogram_size_ = 25000);
    ~DistanceWeightFunction2PPR3();

    void recomputeHistogram(float * hist, double * vec, unsigned int nr_data);
    virtual void recomputeHistogram(float * hist, MatrixXd & mat);
    virtual void recomputeProbs();

    void initComputeModel();
    void setInitialNoise(double stdval);
    virtual void computeModel(double * vec, unsigned int nr_data, unsigned int dim = 1);
	virtual void computeModel(MatrixXd mat);
	virtual VectorXd getProbs(MatrixXd mat);
	virtual double getProb(double d, bool debugg = false);
	virtual double getProbInp(double d, bool debugg = false);
	virtual double getIRLS(double d, bool debugg = false);
	virtual double getInfront(double d, bool debugg = false);
	virtual double getProbInfront(double d, bool debugg = false);
    virtual double getProbInfront(double start, double stop, bool debugg = false);
	virtual double getNoise();
	virtual double getMean();

	virtual double getConvergenceThreshold();
	virtual bool update();
	virtual void reset();
	virtual std::string getString();
	virtual double getInd(double d, bool debugg = false);
    virtual double getDfromInd(double ind, bool debugg = false);

	virtual double getWeight(double invstd, double d,double & infoweight, double & prob, bool debugg = false);
    virtual VectorXd getWeights(std::vector<double > invstd, MatrixXd mat, bool debugg = false);
    virtual void print();
    virtual void setTune();
    virtual double getPower();
};

}

#endif // DistanceWeightFunction2PPR3test_H
