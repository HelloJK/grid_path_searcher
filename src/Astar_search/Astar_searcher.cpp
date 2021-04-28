#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
  for(int i=0; i < GLX_SIZE ; i++)
    for(int j=0; j < GLY_SIZE ; j++)
      for(int k=0; k < GLZ_SIZE ; k++)
        resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                //if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and close list
                if(GridNodeMap[i][j][k]->id == -1)  // visualize nodes in close list only
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d & pt) 
{
  Vector3i idx;
  idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
          min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
          min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  

  return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline bool AstarPathFinder::isInMap(const int& idx_x, const int & idx_y, const int & idx_z) const {
  return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE);
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
  neighborPtrSets.clear();
  edgeCostSets.clear();
  
  int cur_x_idx = currentPtr->index(0);
  int cur_y_idx = currentPtr->index(1);
  int cur_z_idx = currentPtr->index(2);

  for (int dx = -1; dx <= 1; dx++) {
    int new_x = cur_x_idx + dx;
    for (int dy = -1; dy <= 1; dy++) {
      int new_y = cur_y_idx + dy;
      for (int dz = -1; dz <= 1; dz++) {
        int new_z = cur_z_idx + dz;
        if (dx == 0 && dy == 0 && dz == 0) {
          continue;
        }
        
        if (isInMap(new_x, new_y, new_z)) {
          // Eigen::Vector3i idx;
          // idx << new_x, new_y, new_z;
          // Eigen::Vector3d coor;
          // coor = gridIndex2coord(idx);
          // GridNodePtr tmp_ptr = new GridNode(idx, coor);
          // tmp_ptr->cameFrom = currentPtr;
          // tmp_ptr->gScore = currentPtr->gScore + std::sqrt(dx * dx + dy * dy + dz * dz);
          // double hScore = getHeu(tmp_ptr, goalPtr_);
          // tmp_ptr->fScore = tmp_ptr->gScore + hScore;
          
          neighborPtrSets.push_back(GridNodeMap[new_x][new_y][new_z]);
        }
      }
    }
  }
  
  
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
  double dx = std::fabs(node2->index(0)- node1->index(0));
  double dy = std::fabs(node2->index(1)- node1->index(1));
  double dz = std::fabs(node2->index(2)- node1->index(2));
  double D1 = 1;
  double D2 = std::sqrt(2);
  double D3 = std::sqrt(3);

  // int type = 2;
  switch (heuType_) { 
    case 0: {// Euclidean
      return std::sqrt(dx * dx + dy * dy + dz * dz);
      break;
    }
    case 1:{  // manhattan
      // return std::max(std::max(std::fabs(node2->coord(0)- node1->coord(0)), 
      //                          std::fabs(node2->coord(1)- node1->coord(1))), 
      //                 std::fabs(node2->coord(2)- node1->coord(2)));
      return dx + dy + dz;
      break;
    }
    case 2: { // Diagonal
      // 2D
      // dx = abs(node.x - goal.x)
      // dy = abs(node.y - goal.y)
      // return D * (dx + dy) + (D2 - 2 * D) * min(dx, dy)
      
      // 3D
      // dx = absdiff(node.x, goal.x)
      // dy = absdiff(node.y, goal.y)
      // dz = absdiff(node.z, goal.z)
      // dmin = min(dx, dy, dz)
      // dmax = max(dx, dy, dz)
      // dmid = dx + dy + dz - dmin - dmax
      // return (D3 - D2) * dmin + (D2 - D1) * dmid + D1 * dmax

      double dmin = std::min(std::min(dx, dy), dz);
      double dmax = std::max(std::max(dx, dy), dz);
      double dmid = dx + dy + dz - dmin - dmax;
      return (D3 - D2) * dmin + (D2 - D1) * dmid + D1 * dmax;

      break;
    }
    case 3: { // 0 for dijkstra
      return 0;
      break;
    }
    default: {
      return 0;
      break;  
    }
  } 
  return 0;
}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt)
{
  ros::Time time_1 = ros::Time::now();    

  //index of start_point and end_point
  Vector3i start_idx = coord2gridIndex(start_pt);
  Vector3i end_idx   = coord2gridIndex(end_pt);
  goalIdx = end_idx;

  //position of start_point and end_point
  start_pt = gridIndex2coord(start_idx);
  end_pt   = gridIndex2coord(end_idx);

  // ROS_INFO("[Astar] start: %lf %lf %lf", start_pt(0), start_pt(1), start_pt(2));
  // ROS_INFO("[Astar] start_idx: %d %d %d", start_idx(0), start_idx(1), start_idx(2));

  //Initialize the pointers of struct GridNode which represent start node and goal node
  // GridNodePtr startPtr = new GridNode(start_idx, start_pt);
  // GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);
  startPtr_ = GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
  goalPtr_ = GridNodeMap[end_idx(0)][end_idx(1)][end_idx(2)];
  // ROS_INFO("[Astar] start: %lf %lf %lf", startPtr_->coord(0), startPtr_->coord(1), startPtr_->coord(2));
  // ROS_INFO("[Astar] start_idx: %d %d %d", startPtr_->index(0), startPtr_->index(1), startPtr_->index(2));
  // ROS_INFO("[Astar] end: %lf %lf %lf", goalPtr_->coord(0), goalPtr_->coord(1), goalPtr_->coord(2));
  // ROS_INFO("[Astar] end_idx: %d %d %d", goalPtr_->index(0), goalPtr_->index(1), goalPtr_->index(2));
  openSet.clear();
  GridNodePtr currentPtr  = NULL;
  GridNodePtr neighborPtr = NULL;

  //put start node in open set
  startPtr_->gScore = 0;
  startPtr_->fScore = getHeu(startPtr_, goalPtr_);   
  startPtr_->id = 1; 
  startPtr_->coord = start_pt;
  openSet.insert( make_pair(startPtr_->fScore, startPtr_) );
  
  vector<GridNodePtr> neighborPtrSets;
  vector<double> edgeCostSets;

  // this is the main loop
  terminatePtr = nullptr;
  int count = 0;
  while ( !openSet.empty()){
    count++;

    auto begin_iter = openSet.begin();
    currentPtr = begin_iter->second;
    currentPtr->id = -1;
    openSet.erase(begin_iter);

    // ROS_INFO("[Astar] current_idx: %d %d %d", currentPtr->index(0), currentPtr->index(1), currentPtr->index(2));
    
    // if the current node is the goal 
    if( currentPtr->index == goalIdx ){
      ROS_INFO("search succeed");
      ros::Time time_2 = ros::Time::now();
      terminatePtr = currentPtr;
      ROS_WARN("[A*]{sucess} count:%d Time in A*  is %f ms, path cost is %f m", count, (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );            
      return;
    }
    //get the succetion
    AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);
    // if (count < 2) {
    //   cout << "===" << endl;
    //   cout << neighborPtrSets.size() << endl;
    //   for (int i = 0; i < neighborPtrSets.size(); i++) {
    //     cout << neighborPtrSets.at(i)->index(0) << " "<< neighborPtrSets.at(i)->index(1) << " "<< neighborPtrSets.at(i)->index(2) << endl;
    //   }
    //   cout << "=====" << endl;
    // }       
    for(int i = 0; i < (int)neighborPtrSets.size(); i++) {
      neighborPtr = neighborPtrSets.at(i);
      if (!isFree(neighborPtr->index(0), neighborPtr->index(1), neighborPtr->index(2))) {
        continue;
      }
      if(neighborPtr->id == 0){ //discover a new node, which is not in the closed set and open set
        neighborPtr->id = 1;
        neighborPtr->cameFrom = currentPtr;
        double dx = neighborPtr->index(0) - currentPtr->index(0);
        double dy = neighborPtr->index(1) - currentPtr->index(1);
        double dz = neighborPtr->index(2) - currentPtr->index(2);
        neighborPtr->gScore = currentPtr->gScore + std::sqrt(dx * dx + dy * dy + dz * dz);
        neighborPtr->fScore = neighborPtr->gScore + getHeu(neighborPtr, goalPtr_);
        openSet.insert(std::make_pair(neighborPtr->fScore, neighborPtr));
        // continue;
      } else if (neighborPtr->id == 1){
        double dx = neighborPtr->index(0) - currentPtr->index(0);
        double dy = neighborPtr->index(1) - currentPtr->index(1);
        double dz = neighborPtr->index(2) - currentPtr->index(2);
        double gScore = currentPtr->gScore + std::sqrt(dx * dx + dy * dy + dz * dz);
        double hScore = getHeu(neighborPtr, goalPtr_);
        if (gScore + hScore < neighborPtr->fScore) {
          neighborPtr->cameFrom = currentPtr;
          neighborPtr->gScore = gScore;
          neighborPtr->fScore = gScore + hScore;
        }
        // continue;
      } else if (neighborPtr->id == -1) {
        continue;
      }
    }
  }
  //if search fails
  ROS_INFO("search failed, count: %d", count);
  ros::Time time_2 = ros::Time::now();
  if((time_2 - time_1).toSec() > 0.1)
      ROS_WARN("Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}


vector<Vector3d> AstarPathFinder::getPath() 
{   
  vector<Vector3d> path;
  vector<GridNodePtr> gridPath;

  if (terminatePtr == nullptr) {
    return path;
  }
  auto iter_ptr = terminatePtr;
  while (iter_ptr->cameFrom != startPtr_ ) {
    gridPath.push_back(iter_ptr);
    iter_ptr = iter_ptr->cameFrom;
  }
  gridPath.push_back(startPtr_);

  for (auto ptr: gridPath) {
    path.push_back(ptr->coord);
  }
  ROS_INFO("path size: %d", path.size());
      
  reverse(path.begin(),path.end());
  // for (int i = 0; i < path.size(); i++) {
  //   ROS_INFO("[Astar] index:%d : (%lf, %lf, %lf)", i, path.at(i)(0), path.at(i)(1), path.at(i)(2));
  // }

  return path;
}