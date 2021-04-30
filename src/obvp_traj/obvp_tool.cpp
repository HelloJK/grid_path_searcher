#include "../../include/obvp_traj/obvp_tool.h"

using namespace std;
using namespace Eigen;

namespace obvp {

void OBVPtool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
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
}

void OBVPtool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool OBVPtool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d OBVPtool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i OBVPtool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d OBVPtool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double OBVPtool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = DBL_MAX; // this just to initial the optimal_cost, you can delete it 
    
    Eigen::Vector3d target_velocity;
    target_velocity << 0.0, 0.0, 0.0;

    start_position_ = _start_position;
    start_velocity_ = _start_velocity;
    target_position_ = _target_position;
    target_velocity_ = target_velocity;
    
    double delta_x = _target_position[0] - _start_position[0];
    double delta_y = _target_position[1] - _start_position[1];
    double delta_z = _target_position[2] - _start_position[2];
    double delta_vx = _start_velocity[0] + target_velocity[0];
    double delta_vy = _start_velocity[1] + target_velocity[1];
    double delta_vz = _start_velocity[2] + target_velocity[2];
    double delta_vvx = _start_velocity[0] * _start_velocity[0] + target_velocity[0] * target_velocity[0] + 
                       _start_velocity[0] * target_velocity[0];
    double delta_vvy = _start_velocity[1] * _start_velocity[1] + target_velocity[1] * target_velocity[1] + 
                       _start_velocity[1] * target_velocity[1];
    double delta_vvz = _start_velocity[2] * _start_velocity[2] + target_velocity[2] * target_velocity[2] + 
                       _start_velocity[2] * target_velocity[2];

    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    Eigen::VectorXd coeff(5);
    coeff[0] = -36 * (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
    coeff[1] = 24 * (delta_x * delta_vx + delta_y * delta_vy, delta_z * delta_vz);
    coeff[2] = -4 * (delta_vvx + delta_vvy + delta_vvz);
    coeff[3] = 0.0;
    coeff[4] = 1.0;
    
    solver.compute(coeff);
    const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType& r = solver.roots();
    std::vector<double> real_root;
    for (int i = 0; i < r.rows(); ++i) {
      if (r[i].imag() != 0.0) {
          continue;
      }
      real_root.push_back(r[i].real());
    }
    if (real_root.empty()) {
        ROS_WARN("[OBVP] no real root");
        return 100000;
    }
    double optimal_T = 0;
    for (int i = 0; i < real_root.size(); ++i) {
        double T = real_root[i];
        double delta_p_x = _target_position[0] - _start_position[0] - _start_velocity[0] * T;
        double delta_p_y = _target_position[1] - _start_position[1] - _start_velocity[1] * T;
        double delta_p_z = _target_position[2] - _start_position[2] - _start_velocity[2] * T;
        double delta_v_x = target_velocity[0] - _start_velocity[0];
        double delta_v_y = target_velocity[1] - _start_velocity[1];
        double delta_v_z = target_velocity[2] - _start_velocity[2];
        double T_12 = 12.0 / (T * T * T);
        double T_6 = 6.0 / (T * T);
        double T_2 = 2.0 / T;
        double alpha_1 = -T_12 * delta_p_x + T_6 * delta_v_x;
        double alpha_2 = -T_12 * delta_p_y + T_6 * delta_v_y;
        double alpha_3 = -T_12 * delta_p_z + T_6 * delta_v_z;
        double beta_1 = T_6 * delta_p_x - T_2 * delta_v_x;
        double beta_2 = T_6 * delta_p_y - T_2 * delta_v_y;
        double beta_3 = T_6 * delta_p_z - T_2 * delta_v_z;
        double J = T + 1/3 * (alpha_1 * alpha_1 + alpha_2 * alpha_2 + alpha_3 * alpha_3) * T * T * T +
                  (alpha_1 * beta_1 + alpha_2 * beta_2 + alpha_3 * beta_3) * T * T +
                  (beta_1 * beta_1 + beta_2 * beta_2 + beta_3 * beta_3) * T;
        if (optimal_cost > J) {
            optimal_cost = J;
            optimal_T = T;
        }
    }

    optimal_J_ = optimal_cost;
    optimal_T_ = optimal_T;
    ROS_INFO("[OBVP] optimal_J: %lf optimal_T: %lf", optimal_J_, optimal_T_);

    return optimal_cost;
}

void OBVPtool::GetOptimalTraj(std::vector<Eigen::Vector3d>& position, std::vector<Eigen::Vector3d>& velocity) {
    double T = optimal_T_;
    double delta_p_x = target_position_[0] - start_position_[0] - start_velocity_[0] * T;
    double delta_p_y = target_position_[1] - start_position_[1] - start_velocity_[1] * T;
    double delta_p_z = target_position_[2] - start_position_[2] - start_velocity_[2] * T;
    double delta_v_x = target_velocity_[0] - start_velocity_[0];
    double delta_v_y = target_velocity_[1] - start_velocity_[1];
    double delta_v_z = target_velocity_[2] - start_velocity_[2];
    double T_12 = 12 / (T * T * T);
    double T_6 = 6 / (T * T);
    double T_2 = 2 / T;
    double alpha_1 = -T_12 * delta_p_x + T_6 * delta_v_x;
    double alpha_2 = -T_12 * delta_p_y + T_6 * delta_v_y;
    double alpha_3 = -T_12 * delta_p_z + T_6 * delta_v_z;
    double beta_1 = T_6 * delta_p_x - T_2 * delta_v_x;
    double beta_2 = T_6 * delta_p_y - T_2 * delta_v_y;
    double beta_3 = T_6 * delta_p_z - T_2 * delta_v_z;

    Eigen::Vector3d pos;
    Eigen::Vector3d vel;
    double count = 100;
    double time_step = T / count;
    for (double t = 0; t <= T; t += time_step) {
      pos[0] = 1.0/6.0 * alpha_1 * t * t * t + 1.0/2.0 * beta_1 * t * t + start_velocity_[0] * t + start_position_[0];
      pos[1] = 1.0/6.0 * alpha_2 * t * t * t + 1.0/2.0 * beta_2 * t * t + start_velocity_[1] * t + start_position_[1];
      pos[2] = 1.0/6.0 * alpha_3 * t * t * t + 1.0/2.0 * beta_3 * t * t + start_velocity_[2] * t + start_position_[2];
      vel[0] = 1.0/2.0 * alpha_1 * t * t + beta_1 * t + start_velocity_[0];
      vel[1] = 1.0/2.0 * alpha_2 * t * t + beta_2 * t + start_velocity_[1];
      vel[2] = 1.0/2.0 * alpha_3 * t * t + beta_3 * t + start_velocity_[2];
      position.push_back(pos);
      velocity.push_back(vel);
    }
    ROS_INFO("[OVBP] point size: %d", position.size());
    return;
}

} // namespace obvp