#ifndef brainUtil_H
#define brainUtil_H

#include <thread>
#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>

#include "ros/ros.h"
#include "std_msgs/String.h"
#include <sensor_msgs/PointCloud2.h>
#include <string.h>

#include <cv_bridge/cv_bridge.h>
#include <sensor_msgs/image_encodings.h>

#include "eigen_conversions/eigen_msg.h"
#include "tf_conversions/tf_eigen.h"


#include "quasimodo_msgs/segment_model.h"

#include <image_geometry/pinhole_camera_model.h>
#include <sensor_msgs/CameraInfo.h>

#include "modelupdater/ModelUpdater.h"
#include "core/RGBDFrame.h"

#include "metaroom_xml_parser/simple_xml_parser.h"
#include <metaroom_xml_parser/load_utilities.h>

#include <ros/ros.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>

#include <soma_llsd/GetScene.h>
#include <soma_llsd/GetSegment.h>
#include <soma_llsd/InsertScene.h>
#include <soma_llsd_msgs/Segment.h>
#include <quasimodo_conversions/conversions.h>

#include <tf2_ros/transform_broadcaster.h>
#include <tf2_ros/transform_listener.h>

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "std_msgs/String.h"

#include "eigen_conversions/eigen_msg.h"
#include <tf_conversions/tf_eigen.h>

#include "quasimodo_msgs/model.h"
#include "quasimodo_msgs/rgbd_frame.h"
#include "quasimodo_msgs/model_from_frame.h"
#include "quasimodo_msgs/index_frame.h"
#include "quasimodo_msgs/fuse_models.h"
#include "quasimodo_msgs/get_model.h"
#include "quasimodo_msgs/retrieval_query_result.h"
#include "quasimodo_msgs/retrieval_query.h"
#include "quasimodo_msgs/retrieval_result.h"
#include "quasimodo_msgs/recognize.h"
#include "soma_msgs/SOMAObject.h"
#include "soma_manager/SOMAInsertObjs.h"

#include "modelupdater/ModelUpdater.h"
#include "core/RGBDFrame.h"
#include <sensor_msgs/PointCloud2.h>
#include <string.h>

#include "metaroom_xml_parser/simple_xml_parser.h"
#include "metaroom_xml_parser/simple_summary_parser.h"

#include <image_transport/image_transport.h>
#include <cv_bridge/cv_bridge.h>
#include <sensor_msgs/image_encodings.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <map>
#include <thread>
#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include <boost/filesystem.hpp>

namespace quasimodo_brain {

//reglib::Model * processAV(std::string path, bool compute_edges = true, std::string savePath = "");

void readViewXML(std::string roomLogName, std::string xmlFile, std::vector<reglib::RGBDFrame *> & frames, std::vector<Eigen::Matrix4d> & poses, bool compute_edges = true, std::string savePath = "");

void setLargeStack();

reglib::RGBDFrame * getFrame(std::string soma_id, ros::NodeHandle * np);

void recomputeSomaFrame(reglib::RGBDFrame * & frame, ros::NodeHandle * np, std::string savePath = "");

void guaranteeFolder(std::string filepath);

bool isNumber(std::string str);

bool fileExists(std::string path);
int getdirs (std::string dir, std::vector<std::string> & files);
int getdir (std::string dir, std::vector<std::string> & files);
std::string getPoseString(Eigen::Matrix4d pose);

Eigen::Matrix4d getPoseFromString(std::string str);

std::string initSegment(ros::NodeHandle& n, reglib::Model * model);
reglib::Model * getModelFromSegment(ros::NodeHandle& n, std::string segment_id);

soma_llsd_msgs::Scene getScene(ros::NodeHandle& n, reglib::RGBDFrame * frame, std::string current_waypointid = "", std::string roomRunNumber = "");

reglib::Camera * getCam(sensor_msgs::CameraInfo & info);
reglib::RGBDFrame * getFrame(soma_llsd_msgs::Scene & scene);


std::vector<reglib::superpoint> getSuperPoints(std::string path);

std::vector<reglib::superpoint> getRoomSuperPoints(std::string path, std::string savePath);

void transformSuperPoints(std::vector<reglib::superpoint> & spvec, Eigen::Matrix4d cp);

void saveSuperPoints(std::string path, std::vector<reglib::superpoint> & spvec, Eigen::Matrix4d pose, float ratio_keep = 0.1);

std::vector<Eigen::Matrix4d> readPoseXML(std::string xmlFile);

void savePoses(std::string xmlFile, std::vector<Eigen::Matrix4d> poses, int maxposes = -1);

Eigen::Matrix4d getPose(QXmlStreamReader * xmlReader);

int readNumberOfViews(std::string xmlFile);

void writeXml(std::string xmlFile, std::vector<reglib::RGBDFrame *> & frames, std::vector<Eigen::Matrix4d> & poses);

void writePose(QXmlStreamWriter* xmlWriter, Eigen::Matrix4d pose);

void remove_old_seg(std::string sweep_folder, bool backwards = false);

std::string replaceAll(std::string str, std::string from, std::string to);

void cleanPath(std::string & path);
void sortXMLs(std::vector<std::string> & sweeps);

std::vector<Eigen::Matrix4f> getRegisteredViewPosesFromFile(std::string poses_file, int no_transforms);
reglib::Model * loadFromRaresFormat(std::string path);

double getTime();
reglib::Model * getModelFromMSG(quasimodo_msgs::model & msg, bool compute_edges = true);

void addToModelMSG(quasimodo_msgs::model & msg, reglib::Model * model, Eigen::Affine3d rp = Eigen::Affine3d::Identity(), bool addClouds = false);
quasimodo_msgs::model getModelMSG(reglib::Model * model, bool addClouds = false);

std::vector<Eigen::Matrix4f> getRegisteredViewPoses(const std::string& poses_file, const int& no_transforms);

Eigen::Matrix4d getMat(tf::StampedTransform tf);
reglib::Model * load_metaroom_model(std::string sweep_xml, std::string savePath = "");

void segment(std::vector< reglib::Model * > bgs, std::vector< reglib::Model * > models, std::vector< std::vector< cv::Mat > > & internal,
             std::vector< std::vector< cv::Mat > > & external, std::vector< std::vector< cv::Mat > > & dynamic, int debugg = 0,
             std::string savePath = "", std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d> >* model_relative_poses = NULL);
std::vector<reglib::Model *> loadModelsXML(std::string path);
std::vector<reglib::Model *> loadModelsPCDs(std::string path);

}

#endif
