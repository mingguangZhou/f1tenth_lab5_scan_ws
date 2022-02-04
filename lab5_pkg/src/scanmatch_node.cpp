#include <sstream>
#include <string>
#include <cmath>
#include <vector>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/laser_scan.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "std_msgs/msg/color_rgba.hpp"
#include "ackermann_msgs/msg/ackermann_drive_stamped.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "visualization_msgs/msg/marker.hpp"
/// TF: #include <tf/transform_broadcaster.h>
/// CHECK: include needed ROS msg type headers and libraries

#include "scan_matching_skeleton/correspond.h"
#include "scan_matching_skeleton/transform.h"
#include "scan_matching_skeleton/visualization.h"

using namespace std;

const string &TOPIC_SCAN = "/scan";
const string &TOPIC_POS = "/scan_match_location";
const string &TOPIC_RVIZ = "/scan_match_debug";
const string &FRAME_POINTS = "laser";

const float RANGE_LIMIT = 10.0;

const float MAX_ITER = 2.0;
const float MIN_INFO = 0.1;
const float A = (1 - MIN_INFO) / MAX_ITER / MAX_ITER;

class ScanMatch : public rclcpp::Node
{
    // Implement Scan Matching
    // This is just a template, you are free to implement your own node!

private:
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pos_pub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr marker_pub;

    vector<Point> points;
    vector<Point> transformed_points;
    vector<Point> prev_points;
    vector<Correspondence> corresponds;
    vector<vector<int>> jump_table;
    /// TF:
    // Transform prev_trans, curr_trans;
    // tf::TransformBroadcaster br;
    // tf::Transform tr;

    PointVisualizer *points_viz;
    // CorrespondenceVisualizer* corr_viz;

    geometry_msgs::msg::PoseStamped msg;
    Eigen::Matrix3f global_tf;
    std_msgs::msg::ColorRGBA col;

    // void lidar_callback(const sensor_msgs::msg::LaserScan::ConstSharedPtr scan_msg)
    // {

    // }

public:
    ScanMatch() : Node("scanmatch_node")
    {
        pos_pub = this->create_publisher<geometry_msgs::msg::PoseStamped>(TOPIC_POS, 1);
        marker_pub = this->create_publisher<visualization_msgs::msg::Marker>(TOPIC_RVIZ, 1);
        points_viz = new PointVisualizer(marker_pub, "scan_match", FRAME_POINTS);
        // corr_viz = new CorrespondenceVisualizer(marker_pub, "scan_match", FRAME_POINTS);
        global_tf = Eigen::Matrix3f::Identity(3, 3);
    }

    void handleLaserScan(const sensor_msgs::msg::LaserScan::ConstSharedPtr msg)
    {
        readScan(msg);

        //We have nothing to compare to!
        if (prev_points.empty())
        {
            RCLCPP_INFO(this->get_logger(), "%s\n", "First Scan");
            prev_points = points;
            return;
        }

        col.r = 1.0;
        col.b = 0.0;
        col.g = 0.0;
        col.a = 1.0;
        points_viz->addPoints(prev_points, col);

        int count = 0;
        computeJump(jump_table, prev_points);
        RCLCPP_INFO(this->get_logger(), "%s\n", "Starting Optimization");

        // curr_trans = Transform();

        /// TF: while (count < MAX_ITER && (curr_trans != prev_trans || count == 0))
        {
            /// TF: transformPoints(points, curr_trans, transformed_points);

            //************************************************ Find correspondence between points of the current and previous frames  *************** ////
            // **************************************************** getCorrespondence() function is the fast search function and getNaiveCorrespondence function is the naive search option **** ////

            // getCorrespondence(prev_points, transformed_points, points, jump_table, corresponds, A * count * count + MIN_INFO);

            getNaiveCorrespondence(prev_points, transformed_points, points, jump_table, corresponds, A * count * count + MIN_INFO);

            /// TF: prev_trans = curr_trans;
            ++count;

            // **************************************** We update the transforms here ******************************************* ////
            /// TF: updateTransform(corresponds, curr_trans);
        }

        col.r = 0.0;
        col.b = 0.0;
        col.g = 1.0;
        col.a = 1.0;
        points_viz->addPoints(transformed_points, col);
        points_viz->publishPoints();

        RCLCPP_INFO(this->get_logger(), "Count: %i", count);

        /// TF: this->global_tf = global_tf * curr_trans.getMatrix();

        publishPos();
        prev_points = points;
    }

    void readScan(const sensor_msgs::msg::LaserScan::ConstSharedPtr msg)
    {
        float range_min = msg->range_min;
        float range_max = msg->range_max;
        float angle_min = msg->angle_min;
        float angle_increment = msg->angle_increment;

        const vector<float> &ranges = msg->ranges;
        points.clear();

        for (int i = 0; i < ranges.size(); ++i)
        {
            float range = ranges.at(i);
            if (range > RANGE_LIMIT)
            {
                continue;
            }
            if (!isnan(range) && range >= range_min && range <= range_max)
            {
                points.push_back(Point(range, angle_min + angle_increment * i));
            }
        }
    }

    void publishPos()
    {   
        /// TF:
        // msg.pose.position.x = global_tf(0, 2);
        // msg.pose.position.y = global_tf(1, 2);
        // msg.pose.position.z = 0;
        // tf::Matrix3x3 tf3d;
        // tf3d.setValue(static_cast<double>(global_tf(0, 0)), static_cast<double>(global_tf(0, 1)), 0,
        //             static_cast<double>(global_tf(1, 0)), static_cast<double>(global_tf(1, 1)), 0, 0, 0, 1);

        // tf::Quaternion q;
        // tf3d.getRotation(q);
        // msg.pose.orientation.x = q.x();
        // msg.pose.orientation.y = q.y();
        // msg.pose.orientation.z = q.z();
        // msg.pose.orientation.w = q.w();
        // msg.header.frame_id = "laser";
        // msg.header.stamp = this->get_clock()->now();
        // pos_pub.publish(msg);
        // tr.setOrigin(tf::Vector3(global_tf(0, 2), global_tf(1, 2), 0));
        // tr.setRotation(q);
        // br.sendTransform(tf::StampedTransform(tr, this->get_clock()->now(), "map", "laser"));
    }

    ~ScanMatch() {}
};
int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<ScanMatch>());
    rclcpp::shutdown();
    return 0;
}