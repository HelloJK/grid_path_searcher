#include <iostream>
#include <fstream>
#include <math.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>

#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include "Astar_searcher.h"
#include "JPS_searcher.h"
#include "backward.hpp"
#include "obvp_tool.h"

using namespace std;
using namespace Eigen;

namespace backward {
backward::SignalHandling sh;
}

// simulation param from launch file
double _resolution, _inv_resolution, _cloud_margin;
double _x_size, _y_size, _z_size;    

// useful global variables
bool _has_map   = false;

Vector3d _start_pt, _start_velocity;
Vector3d _map_lower, _map_upper;
int _max_x_id, _max_y_id, _max_z_id;

// ros related
ros::Subscriber _map_sub, _pts_sub;
ros::Publisher  _grid_path_vis_pub, _visited_nodes_vis_pub, _grid_map_vis_pub, _path_vis_pub;

AstarPathFinder * _astar_path_finder     = new AstarPathFinder();
JPSPathFinder   * _jps_path_finder       = new JPSPathFinder();

obvp::OBVPtool * _obvp_tool     = new obvp::OBVPtool();
obvp::TrajectoryStatePtr *** TraLibrary;

// Integral parameter
double _max_input_acc     = 1.0;
int    _discretize_step   = 2;
double _time_interval     = 1.25;
int    _time_step         = 50;

void rcvWaypointsCallback(const nav_msgs::Path & wp);
void rcvWaypointsCallback2(const nav_msgs::Path & wp);
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map);


void visGridPath( vector<Vector3d> nodes, bool is_use_jps );
void visVisitedNode( vector<Vector3d> nodes );
void visTraLibrary(obvp::TrajectoryStatePtr *** TraLibrary);
void AstarPathFinding(const Vector3d start_pt, const Vector3d target_pt);
void trajectoryLibrary(const Eigen::Vector3d start_pt, const Eigen::Vector3d start_velocity, const Eigen::Vector3d target_pt);

void rcvWaypointsCallback(const nav_msgs::Path & wp)
{     
    if( wp.poses[0].pose.position.z < 0.0 || _has_map == false )
        return;

    Vector3d target_pt;
    target_pt << wp.poses[0].pose.position.x,
                 wp.poses[0].pose.position.y,
                 wp.poses[0].pose.position.z;

    ROS_INFO("[node] receive the planning target");
    ROS_INFO("[node] start: %lf %lf %lf", _start_pt(0), _start_pt(1), _start_pt(2));
    ROS_INFO("[node] end: %lf %lf %lf", target_pt(0), target_pt(1), target_pt(2));
    AstarPathFinding(_start_pt, target_pt); 
}

void rcvWaypointsCallback2(const nav_msgs::Path & wp) {
    if( wp.poses[0].pose.position.z < 0.0 || _has_map == false )
        return;

    Vector3d target_pt;
    target_pt << wp.poses[0].pose.position.x,
                 wp.poses[0].pose.position.y,
                 wp.poses[0].pose.position.z;

    ROS_INFO("[node] receive the planning target");
    trajectoryLibrary(_start_pt, _start_velocity, target_pt);
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map)
{   
    if(_has_map ) return;

    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::PointCloud<pcl::PointXYZ> cloud_vis;
    sensor_msgs::PointCloud2 map_vis;

    pcl::fromROSMsg(pointcloud_map, cloud);
    
    if( (int)cloud.points.size() == 0 ) return;

    pcl::PointXYZ pt;
    for (int idx = 0; idx < (int)cloud.points.size(); idx++)
    {    
        pt = cloud.points[idx];        

        // set obstalces into grid map for path planning
        _astar_path_finder->setObs(pt.x, pt.y, pt.z);
        _jps_path_finder->setObs(pt.x, pt.y, pt.z);
        _obvp_tool->setObs(pt.x, pt.y, pt.z);

        // for visualize only
        Vector3d cor_round = _astar_path_finder->coordRounding(Vector3d(pt.x, pt.y, pt.z));
        pt.x = cor_round(0);
        pt.y = cor_round(1);
        pt.z = cor_round(2);
        cloud_vis.points.push_back(pt);
    }

    cloud_vis.width    = cloud_vis.points.size();
    cloud_vis.height   = 1;
    cloud_vis.is_dense = true;

    pcl::toROSMsg(cloud_vis, map_vis);

    map_vis.header.frame_id = "/world";
    _grid_map_vis_pub.publish(map_vis);

    _has_map = true;
}

void AstarPathFinding(const Vector3d start_pt, const Vector3d target_pt)
{
    {
      ROS_INFO("[node] ==========================Euclidean");
      //Call A* to search for a path
      _astar_path_finder->setHeutype(0);
      _astar_path_finder->AstarGraphSearch(start_pt, target_pt);

      //Retrieve the path
      auto grid_path     = _astar_path_finder->getPath();
      auto visited_nodes = _astar_path_finder->getVisitedNodes();

      //Visualize the result
      visGridPath (grid_path, false);
      visVisitedNode(visited_nodes);

      //Reset map for next call
      _astar_path_finder->resetUsedGrids();
      ros::Duration(3.0).sleep();
    }

    {
      ROS_INFO("[node] ==========================manhattan");
      //Call A* to search for a path
      _astar_path_finder->setHeutype(1);
      _astar_path_finder->AstarGraphSearch(start_pt, target_pt);

      //Retrieve the path
      auto grid_path     = _astar_path_finder->getPath();
      auto visited_nodes = _astar_path_finder->getVisitedNodes();

      //Visualize the result
      visGridPath (grid_path, false);
      visVisitedNode(visited_nodes);

      //Reset map for next call
      _astar_path_finder->resetUsedGrids();
    }

    {
      ROS_INFO("[node] ==========================Diagonal");
      //Call A* to search for a path
      _astar_path_finder->setHeutype(2);
      _astar_path_finder->AstarGraphSearch(start_pt, target_pt);

      //Retrieve the path
      auto grid_path     = _astar_path_finder->getPath();
      auto visited_nodes = _astar_path_finder->getVisitedNodes();

      //Visualize the result
      visGridPath (grid_path, false);
      visVisitedNode(visited_nodes);

      //Reset map for next call
      _astar_path_finder->resetUsedGrids();
    }

    {
      ROS_INFO("[node] ==========================dijkstra");
      //Call A* to search for a path
      _astar_path_finder->setHeutype(3);
      _astar_path_finder->AstarGraphSearch(start_pt, target_pt);

      //Retrieve the path
      auto grid_path     = _astar_path_finder->getPath();
      auto visited_nodes = _astar_path_finder->getVisitedNodes();

      //Visualize the result
      visGridPath (grid_path, false);
      visVisitedNode(visited_nodes);

      //Reset map for next call
      _astar_path_finder->resetUsedGrids();
    }

    //_use_jps = 0 -> Do not use JPS
    //_use_jps = 1 -> Use JPS
    //you just need to change the #define value of _use_jps
#define _use_jps 0
#if _use_jps
    {
        //Call JPS to search for a path
        _jps_path_finder -> JPSGraphSearch(start_pt, target_pt);

        //Retrieve the path
        auto grid_path     = _jps_path_finder->getPath();
        auto visited_nodes = _jps_path_finder->getVisitedNodes();

        //Visualize the result
        visGridPath   (grid_path, _use_jps);
        visVisitedNode(visited_nodes);

        //Reset map for next call
        _jps_path_finder->resetUsedGrids();
    }
#endif
}


void trajectoryLibrary(const Vector3d start_pt, const Vector3d start_velocity, const Vector3d target_pt)
{
    Vector3d acc_input;
    Vector3d pos,vel;
    int a =0 ;
    int b =0 ;
    int c =0 ;

    double min_Cost = 100000.0;
    double Trajctory_Cost;
    TraLibrary  = new obvp::TrajectoryStatePtr ** [_discretize_step + 1];     //recored all trajectories after input

    for(int i=0; i <= _discretize_step; i++){           //acc_input_ax
        TraLibrary[i] = new obvp::TrajectoryStatePtr * [_discretize_step + 1];
        for(int j=0;j <= _discretize_step; j++){        //acc_input_ay
            TraLibrary[i][j] = new obvp::TrajectoryStatePtr [_discretize_step + 1];
            for(int k=0; k <= _discretize_step; k++){   //acc_input_az
                vector<Vector3d> Position;
                vector<Vector3d> Velocity;
                acc_input(0) = double(-_max_input_acc + i * (2 * _max_input_acc / double(_discretize_step)) );
                acc_input(1) = double(-_max_input_acc + j * (2 * _max_input_acc / double(_discretize_step)) );
                acc_input(2) = double( k * (2 * _max_input_acc / double(_discretize_step) ) + 0.1);                          //acc_input_az >0.1
                
                pos(0) = start_pt(0);
                pos(1) = start_pt(1);
                pos(2) = start_pt(2);
                vel(0) = start_velocity(0);
                vel(1) = start_velocity(1);
                vel(2) = start_velocity(2);
                Position.push_back(pos);
                Velocity.push_back(vel);

                bool collision = false;
                double delta_time;
                delta_time = _time_interval / double(_time_step);
                for(int step=0 ; step<=_time_step ; step ++){
                    pos = pos + vel * delta_time + 0.5 * acc_input * delta_time * delta_time;
                    vel = vel + delta_time * acc_input;
                    
                    Position.push_back(pos);
                    Velocity.push_back(vel);
                    double coord_x = pos(0);
                    double coord_y = pos(1);
                    double coord_z = pos(2);
                    //check if if the trajectory face the obstacle
                    if(_obvp_tool->isObsFree(coord_x,coord_y,coord_z) != 1){
                        collision = true;
                    }
                }
                Trajctory_Cost = _obvp_tool -> OptimalBVP(pos,vel,target_pt);
                std::vector<Eigen::Vector3d> optimal_path;
                std::vector<Eigen::Vector3d> optimal_vel;
                _obvp_tool->GetOptimalTraj(optimal_path, optimal_vel);
                visGridPath(optimal_path, false);

                //input the trajetory in the trajectory library
                TraLibrary[i][j][k] = new obvp::TrajectoryState(Position,Velocity,Trajctory_Cost);
                
                //if there is not any obstacle in the trajectory we need to set 'collision_check = true', so this trajectory is useable
                if(collision)
                    TraLibrary[i][j][k]->setCollisionfree();
                
                //record the min_cost in the trajectory Library, and this is the part pf selecting the best trajectory cloest to the planning traget
                if (Trajctory_Cost<min_Cost && TraLibrary[i][j][k]->collision_check == false) {
                    a = i;
                    b = j;
                    c = k;
                    min_Cost = Trajctory_Cost;
                }
            }
        }
    }
    TraLibrary[a][b][c] -> setOptimal();
    visTraLibrary(TraLibrary);
    return;
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "demo_node");
    ros::NodeHandle nh("~");

    _map_sub  = nh.subscribe( "map",       1, rcvPointCloudCallBack );
    // _pts_sub  = nh.subscribe( "waypoints", 1, rcvWaypointsCallback );
    _pts_sub  = nh.subscribe( "waypoints", 1, rcvWaypointsCallback2 );

    _grid_map_vis_pub             = nh.advertise<sensor_msgs::PointCloud2>("grid_map_vis", 1);
    _grid_path_vis_pub            = nh.advertise<visualization_msgs::Marker>("grid_path_vis", 1);
    _visited_nodes_vis_pub        = nh.advertise<visualization_msgs::Marker>("visited_nodes_vis",1);
    _path_vis_pub                 = nh.advertise<visualization_msgs::MarkerArray>("RRTstar_path_vis",1);

    nh.param("map/cloud_margin",  _cloud_margin, 0.0);
    nh.param("map/resolution",    _resolution,   0.2);
    
    nh.param("map/x_size",        _x_size, 50.0);
    nh.param("map/y_size",        _y_size, 50.0);
    nh.param("map/z_size",        _z_size, 5.0 );
    
    nh.param("planning/start_x",  _start_pt(0),  0.0);
    nh.param("planning/start_y",  _start_pt(1),  0.0);
    nh.param("planning/start_z",  _start_pt(2),  0.0);

    _map_lower << - _x_size/2.0, - _y_size/2.0,     0.0;
    _map_upper << + _x_size/2.0, + _y_size/2.0, _z_size;
    
    _inv_resolution = 1.0 / _resolution;
    
    _max_x_id = (int)(_x_size * _inv_resolution);
    _max_y_id = (int)(_y_size * _inv_resolution);
    _max_z_id = (int)(_z_size * _inv_resolution);
    ROS_INFO("[node] size: %d %d %d", _max_x_id, _max_y_id, _max_z_id);

    _astar_path_finder  = new AstarPathFinder();
    _astar_path_finder  -> initGridMap(_resolution, _map_lower, _map_upper, _max_x_id, _max_y_id, _max_z_id);

    _jps_path_finder    = new JPSPathFinder();
    _jps_path_finder    -> initGridMap(_resolution, _map_lower, _map_upper, _max_x_id, _max_y_id, _max_z_id);

    _obvp_tool  = new obvp::OBVPtool();
    _obvp_tool  -> initGridMap(_resolution, _map_lower, _map_upper, _max_x_id, _max_y_id, _max_z_id);
    
    ros::Rate rate(100);
    bool status = ros::ok();
    while(status) 
    {
        ros::spinOnce();      
        status = ros::ok();
        rate.sleep();
    }

    delete _astar_path_finder;
    delete _jps_path_finder;
    return 0;
}

void visGridPath( vector<Vector3d> nodes, bool is_use_jps )
{   
    visualization_msgs::Marker node_vis; 
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    
    if(is_use_jps)
        node_vis.ns = "demo_node/jps_path";
    else
        node_vis.ns = "demo_node/astar_path";

    node_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;

    if(is_use_jps){
        node_vis.color.a = 1.0;
        node_vis.color.r = 1.0;
        node_vis.color.g = 0.0;
        node_vis.color.b = 0.0;
    }
    else{
        node_vis.color.a = 1.0;
        node_vis.color.r = 0.0;
        node_vis.color.g = 1.0;
        node_vis.color.b = 0.0;
    }


    node_vis.scale.x = _resolution;
    node_vis.scale.y = _resolution;
    node_vis.scale.z = _resolution;

    geometry_msgs::Point pt;
    for(int i = 0; i < int(nodes.size()); i++)
    {
        Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    _grid_path_vis_pub.publish(node_vis);
}

void visVisitedNode( vector<Vector3d> nodes )
{   
    visualization_msgs::Marker node_vis; 
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "demo_node/expanded_nodes";
    node_vis.type = visualization_msgs::Marker::CUBE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;
    node_vis.color.a = 0.5;
    node_vis.color.r = 0.0;
    node_vis.color.g = 0.0;
    node_vis.color.b = 1.0;

    node_vis.scale.x = _resolution;
    node_vis.scale.y = _resolution;
    node_vis.scale.z = _resolution;

    geometry_msgs::Point pt;
    for(int i = 0; i < int(nodes.size()); i++)
    {
        Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    _visited_nodes_vis_pub.publish(node_vis);
}

void visTraLibrary(obvp::TrajectoryStatePtr *** TraLibrary)
{
    double _resolution = 0.2;
    visualization_msgs::MarkerArray  LineArray;
    visualization_msgs::Marker       Line;

    Line.header.frame_id = "world";
    Line.header.stamp    = ros::Time::now();
    Line.ns              = "demo_node/TraLibrary";
    Line.action          = visualization_msgs::Marker::ADD;
    Line.pose.orientation.w = 1.0;
    Line.type            = visualization_msgs::Marker::LINE_STRIP;
    Line.scale.x         = _resolution/5;

    Line.color.r         = 0.0;
    Line.color.g         = 0.0;
    Line.color.b         = 1.0;
    Line.color.a         = 1.0;

    int marker_id = 0;

    for(int i = 0; i <= _discretize_step; i++){
        for(int j = 0; j<= _discretize_step;j++){  
            for(int k = 0; k<= _discretize_step;k++){
                if(TraLibrary[i][j][k]->collision_check == false){
                    if(TraLibrary[i][j][k]->optimal_flag == true){
                        Line.color.r         = 0.0;
                        Line.color.g         = 1.0;
                        Line.color.b         = 0.0;
                        Line.color.a         = 1.0;
                    }else{
                        Line.color.r         = 0.0;
                        Line.color.g         = 0.0;
                        Line.color.b         = 1.0;
                        Line.color.a         = 1.0;
                    }
                }else{
                    Line.color.r         = 1.0;
                    Line.color.g         = 0.0;
                    Line.color.b         = 0.0;
                    Line.color.a         = 1.0;
                }
                   Line.points.clear();
                    geometry_msgs::Point pt;
                    Line.id = marker_id;
                    for(int index = 0; index < int(TraLibrary[i][j][k]->Position.size());index++){
                        Vector3d coord = TraLibrary[i][j][k]->Position[index];
                        pt.x = coord(0);
                        pt.y = coord(1);
                        pt.z = coord(2);
                        Line.points.push_back(pt);
                    }
                    LineArray.markers.push_back(Line);
                    _path_vis_pub.publish(LineArray);
                    ++marker_id; 
            }
        }
    }    
}