#include <ros/ros.h>
#include <mongodb_store/message_store.h>
#include <quasimodo_msgs/fused_world_state_object.h>
#include <std_msgs/Empty.h>
#include <sensor_msgs/PointCloud2.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/common/common.h>
#include <object_3d_benchmark/surfel_type.h>
#include <object_3d_benchmark/benchmark_visualization.h>
#include <object_3d_benchmark/surfel_renderer.h>
#include <eigen_conversions/eigen_msg.h>

const mongo::BSONObj EMPTY_BSON_OBJ;

using namespace std;
using PointT = pcl::PointXYZRGB;
using CloudT = pcl::PointCloud<PointT>;
using SurfelT = SurfelType;
using SurfelCloudT = pcl::PointCloud<SurfelT>;

class visualization_server {
public:

    ros::NodeHandle n;
    mongodb_store::MessageStoreProxy db_client;
    ros::Subscriber sub;
    size_t counter;

    visualization_server() : db_client(n, "quasimodo", "world_state"), counter(0)
    {

        sub = n.subscribe("/model/added_to_db", 1, &visualization_server::callback, this);

        callback(std_msgs::Empty());

        ros::spin();
    }

    void callback(const std_msgs::Empty& empty_msg)
    {
        std::vector<boost::shared_ptr<quasimodo_msgs::fused_world_state_object> > messages;
        mongo::BSONObj meta_query = BSON("removed_at" << "");
        db_client.query(messages); //, EMPTY_BSON_OBJ, meta_query);

        vector<SurfelCloudT::Ptr> results;
        vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > first_transforms;
        for (boost::shared_ptr<quasimodo_msgs::fused_world_state_object>& msg : messages) {
            if (!msg->removed_at.empty()) {
                continue;
            }
            if (msg->inserted_at.empty()) {
                cout << "No insertion for msg with ID: " << msg->object_id << endl;
            }
            cout << "Adding " << counter << ":th point cloud..." << endl;
            SurfelCloudT::Ptr cloud(new SurfelCloudT);
            pcl::fromROSMsg(msg->surfel_cloud, *cloud);

            if (cloud->size() < 100) {
                cout << "Cloud too small, skipping" << endl;
                ++counter;
                continue;
            }

            Eigen::Affine3d e;

            Eigen::Matrix4f T;
            if (msg->transforms.size() > 0) {
                tf::poseMsgToEigen(msg->transforms[0], e);
                T = e.matrix().cast<float>();
            }
            else {
                T.setIdentity();
            }

            cout << "Using transform: \n" << T << endl;

            results.push_back(cloud);
            first_transforms.push_back(T);

            ++counter;
        }

        cv::Mat visualization_image;
        vector<cv::Mat> individual_images;

        tie(visualization_image, individual_images) = benchmark_retrieval::make_image(results, first_transforms);

        //cv::imshow("Visualization image", visualization_image);
        //cv::waitKey();
        cv::imwrite("quasimodo_db_vis.png", visualization_image);
    }

};

int main(int argc, char** argv)
{
    ros::init(argc, argv, "visualize_database_node");

    visualization_server vs;

    return 0;
}
