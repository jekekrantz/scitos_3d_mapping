#ifndef reglibRGBDFrame_H
#define reglibRGBDFrame_H

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h> 
#include <chrono>

#include <Eigen/Dense>

// PCL specific includes
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/features/integral_image_normal.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>

#include "Util.h"
#include "Camera.h"
#include "../weightfunctions/DistanceWeightFunction2.h"

namespace reglib
{
	class RGBDFrame{
		public:
		std::string keyval;
		Camera * camera;
		unsigned long id;
		double capturetime;
		Eigen::Matrix4d pose;
		int sweepid;

		//float * rgbdata;
		cv::Mat det_dilate;
		cv::Mat rgb;
		cv::Mat depth;
		cv::Mat normals;
		cv::Mat depthedges;
		int * labels;
		int nr_labels;

		std::vector< std::vector<double> > connections;
		std::vector< std::vector<double> > intersections;

		RGBDFrame();
		RGBDFrame(Camera * camera_,cv::Mat rgb_, cv::Mat depth_, double capturetime_ = 0, Eigen::Matrix4d pose_ = Eigen::Matrix4d::Identity(), bool compute_normals = true);
		~RGBDFrame();

		void show(bool stop = false);
		pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr getPCLcloud();
        void savePCD(std::string path = "cloud.pcd", Eigen::Matrix4d pose = Eigen::Matrix4d::Identity());

        void save(std::string path = "");
		RGBDFrame * clone();
		static RGBDFrame * load(Camera * cam, std::string path);
		std::vector<ReprojectionResult> getReprojections(std::vector<superpoint> & spvec, Eigen::Matrix4d cp);
	};
}

#endif // reglibRGBDFrame_H
