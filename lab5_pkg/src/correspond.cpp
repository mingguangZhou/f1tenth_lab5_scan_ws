#include "scan_matching_skeleton/correspond.h"
#include "cmath"
#include <iostream>

using namespace std;

const int UP_SMALL = 0;
const int UP_BIG = 1;
const int DOWN_SMALL = 2;
const int DOWN_BIG = 3;

void getNaiveCorrespondence(vector<Point> &old_points, vector<Point> &trans_points, vector<Point> &points,
                            vector<vector<int>> &jump_table, vector<Correspondence> &c, float prob)
{

  c.clear();
  int last_best = -1;
  const int n = trans_points.size();
  const int m = old_points.size();
  int min_index = 0;
  int second_min_index = 0;

  //Do for each point
  for (int ind_trans = 0; ind_trans < n; ++ind_trans)
  {
    float min_dist = 100000.00;
    for (int ind_old = 0; ind_old < m; ++ind_old)
    {
      // float dist = old_points[ind_trans].distToPoint2(&trans_points[ind_old]);
      float dist = trans_points[ind_trans].distToPoint2(&old_points[ind_old]);
      if (dist < min_dist)
      {
        min_dist = dist;
        min_index = ind_old;
        if (ind_old == 0)
        {
          second_min_index = ind_old + 1;
        }
        else
        {
          second_min_index = ind_old - 1;
        }
      }
    }
    c.push_back(Correspondence(&trans_points[ind_trans], &points[ind_trans], &old_points[min_index], &old_points[second_min_index]));
  }
}

void getCorrespondence(vector<Point> &old_points, vector<Point> &trans_points, vector<Point> &points,
                       vector<vector<int>> &jump_table, vector<Correspondence> &c, float prob)
{

  // Written with inspiration from: https://github.com/AndreaCensi/gpc/blob/master/c/gpc.c
  // use helper functions and structs in transform.h and correspond.h
  // input : old_points : vector of struct points containing the old points (points of the previous frame)
  // input : trans_points : vector of struct points containing the new points transformed to the previous frame using the current estimated transform
  // input : points : vector of struct points containing the new points
  // input : jump_table : jump table computed using the helper functions from the transformed and old points
  // input : c: vector of struct correspondences . This is a refernece which needs to be updated in place and return the new correspondences to calculate the transforms.
  // output : c; update the correspondence vector in place which is provided as a reference. you need to find the index of the best and the second best point.
  //Initializecorrespondences
  c.clear();
  int last_best = -1;
  const int trans_size = trans_points.size();
  const int old_size = old_points.size();

  //Do for each point
  for (int ind_trans = 0; ind_trans < trans_size; ++ind_trans)
  {
    // Implement Fast Correspondence Search

    // initialization:
    // current best match and its dist
    int best = -1;
    float min_dist = std::numeric_limits<float>::max();
    // Approximated index in scan corresponding to current point 
    // trans_points[ind_trans].wrapTheta();
    // old_points[0].wrapTheta();
    int we_start_at = (last_best != -1)?(last_best + 1):0;
    // int we_start_at = (last_best != -1)? last_best:0;
    // Search is conducted in two directions: up and down
    // Distance of last point examined in the up (down) direction.
    float last_dist_up = INFINITY, last_dist_down = INFINITY;
    float dist = INFINITY;
    // True if search is finished in the up (down) direction.
    bool up_stopped = false, down_stopped = false;

    // start from last best, search upwards
    int current = we_start_at;

    int loop_counter = 0;
    const int max_loops = 500; // Arbitrary large number
    
    // search in up direction
    while ((!up_stopped) && (current < old_size))
    {
      // if (++loop_counter > max_loops)
      // {
      //     std::cerr << "Error: Infinite loop detected in correspondence search!" << std::endl;
      //     break;
      // }
      
      // Distance from this trans point to the current 'up' point
      dist = trans_points[ind_trans].distToPoint2(&old_points[current]);
      // If it is less than the best point, up is our best guess so far.
      if (dist < min_dist) {best = current, min_dist = dist;}
      
      // early stop check
      float theta_delta = old_points[current].theta - trans_points[ind_trans].theta;
      float min_dist_up = sin(theta_delta) * trans_points[ind_trans].r;
      // If going up we can’t make better than best_dist, then we stop searching in the "up" direction
      if ((min_dist_up * min_dist_up) > min_dist)
      {
          up_stopped = true;
          // std::cout << "[DEBUG] Early stopping upward search." << std::endl;
          continue;
      }
        
      // Then we can implement the jump tables optimization.
      current = // Next point to examine is...
        (old_points[current].r < trans_points[ind_trans].r) ? // is current scan point longer?
        jump_table[current][UP_BIG] // then jump to a further point
        :jump_table[current][UP_SMALL]; // else, to a closer one.
      
    } // end of up

    // back to initial starting point one down and search downwards
    current = we_start_at - 1;

    // repeat in down direction
    // This is the specular part of the previous chunk of code.
    while ((!down_stopped) && (current >= 0))
    {
      dist = trans_points[ind_trans].distToPoint2(&old_points[current]);
      if (dist < min_dist) { best = current; min_dist = dist;}

      float theta_delta = trans_points[ind_trans].theta - old_points[current].theta;
      float min_dist_down = sin(theta_delta) * trans_points[ind_trans].r;

      if ((min_dist_down * min_dist_down) > min_dist)
      {
        down_stopped = true;
        // std::cout << "[DEBUG] Early stopping downward search." << std::endl;
        continue;
      }

      current = (old_points[current].r < trans_points[ind_trans].r) ? 
              jump_table[current][DOWN_BIG] : jump_table[current][DOWN_SMALL];
    }
    
    // Find the second best point
    int second_best = (best > 0) ? best - 1 : best + 1;

    c.push_back(Correspondence(&trans_points[ind_trans], &points[ind_trans], &old_points[best], &old_points[second_best]));

    // Update last best index
    last_best = best;
  }
}

void computeJump(vector<vector<int>> &table, vector<Point> &points)
{
  table.clear();
  int n = points.size();
  for (int i = 0; i < n; ++i)
  {
    vector<int> v = {n, n, -1, -1};
    for (int j = i + 1; j < n; ++j)
    {
      if (points[j].r < points[i].r)
      {
        v[UP_SMALL] = j;
        break;
      }
    }
    for (int j = i + 1; j < n; ++j)
    {
      if (points[j].r > points[i].r)
      {
        v[UP_BIG] = j;
        break;
      }
    }
    for (int j = i - 1; j >= 0; --j)
    {
      if (points[j].r < points[i].r)
      {
        v[DOWN_SMALL] = j;
        break;
      }
    }
    for (int j = i - 1; j >= 0; --j)
    {
      if (points[j].r > points[i].r)
      {
        v[DOWN_BIG] = j;
        break;
      }
    }
    table.push_back(v);
  }
}
