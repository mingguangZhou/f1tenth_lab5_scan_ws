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

// void getCorrespondence(vector<Point> &old_points, vector<Point> &trans_points, vector<Point> &points,
//                        vector<vector<int>> &jump_table, vector<Correspondence> &c, float prob,
//                        float angle_increment_rad)
// {

//   // Written with inspiration from: https://github.com/AndreaCensi/gpc/blob/master/c/gpc.c
//   // use helper functions and structs in transform.h and correspond.h
//   // input : old_points : vector of struct points containing the old points (points of the previous frame)
//   // input : trans_points : vector of struct points containing the new points transformed to the previous frame using the current estimated transform
//   // input : points : vector of struct points containing the new points
//   // input : jump_table : jump table computed using the helper functions from the transformed and old points
//   // input : c: vector of struct correspondences . This is a refernece which needs to be updated in place and return the new correspondences to calculate the transforms.
//   // output : c; update the correspondence vector in place which is provided as a reference. you need to find the index of the best and the second best point.
//   //Initializecorrespondences
//   c.clear();
//   int last_best = -1;
//   const int trans_size = trans_points.size();
//   const int old_size = old_points.size();

//   //Do for each point
//   for (int ind_trans = 0; ind_trans < min(old_size, trans_size); ++ind_trans)
//   {
//     // Implement Fast Correspondence Search

//     // initialization:
//     // current best match and its dist
//     int best = -1;
//     float min_dist = 100000.00;
//     // Approximated index in scan corresponding to current point 
//     trans_points[ind_trans].wrapTheta();
//     old_points[0].wrapTheta();
//     int start_index;
//     if (trans_points[ind_trans].theta > old_points[0].theta)
//     {
//       start_index = (trans_points[ind_trans].theta - old_points[0].theta) / angle_increment_rad;
//     } else {
//       start_index = 0;
//     }
//     // If last match was succesful, then start at that index + 1
//     int we_start_at = (last_best != -1)?(last_best + 1):start_index;
//     // Search is conducted in two directions: up and down
//     int up = we_start_at+1, down = we_start_at;
//     // Distance of last point examined in the up (down) direction.
//     float last_dist_up = INFINITY, last_dist_down = INFINITY;
//     // True if search is finished in the up (down) direction.
//     bool up_stopped = false, down_stopped = false;

//     int loop_counter = 0;
//     const int max_loops = 500; // Arbitrary large number

//     // select direction to search
//     // Until the search is stopped in both directions...
//     while (!(up_stopped && down_stopped))
//     {
//       if (++loop_counter > max_loops)
//       {
//           // std::cerr << "Error: Infinite loop detected in correspondence search!" << std::endl;
//           break;
//       }

//       // Start of loop
//       // std::cout << "[DEBUG] Start of loop: up = " << up << ", down = " << down
//       //           << ", up_stopped = " << up_stopped
//       //           << ", down_stopped = " << down_stopped 
//       //           << ", old_size = " << old_size << std::endl;
      
      
//       // Should we try to explore up or down?
//       bool now_up = false;
//       if (!up_stopped && !down_stopped)
//       {
//           now_up = (last_dist_up < last_dist_down);
//       }
//       else if (!up_stopped)
//       {
//           now_up = true;
//       }
//       else if (!down_stopped)
//       {
//           now_up = false;
//       }
//       // std::cout << "[DEBUG] Direction chosen: " << (now_up ? "up" : "down") << std::endl;
//       // std::cout << "[DEBUG] last_dist_up: " << last_dist_up 
//       //           << ", last_dist_down: " << last_dist_down 
//       //           << ", now_up = " << now_up << std::endl;
   
//       // Now two symmetric chunks of code, the now_up and the !now_up
//       // search in up direction
//       if(now_up)
//       {
//         // If we have finished the points to search, we stop.
//         if (up >= old_size) 
//         {
//           up_stopped = true;
//           // std::cout << "[DEBUG] Stopping upward search: up = " << up << std::endl; 
//           continue;
//         }
//         // Distance from current point to the 'up' point
//         last_dist_up = trans_points[ind_trans].distToPoint2(&old_points[up]);
//         // If it is less than the best point, up is our best guess so far.
//         if (last_dist_up < min_dist) {best = up, min_dist = last_dist_up;}
//         if (up > start_index)
//         {
//           // If we are moving away from start_cell we can compute a
//           // bound for early stopping. Currently our best point has distance
//           // best_dist; we can compute the minimum distance to current point for points j > up
//           float theta_delta = old_points[up].theta - trans_points[ind_trans].theta;
//           float min_dist_up = sin(theta_delta) * trans_points[ind_trans].r;
//           // If going up we can’t make better than best_dist, then we stop searching in the "up" direction
//           if ((min_dist_up * min_dist_up) > min_dist)
//           {
//               up_stopped = true;
//               // std::cout << "[DEBUG] Early stopping upward search." << std::endl;
//               continue;
//           }
          
//           // If we are moving away, then we can implement the jump tables optimization.
//           up = // Next point to examine is...
//             (old_points[up].r < trans_points[ind_trans].r) ? // is current scan point longer?
//             jump_table[up][UP_BIG] // then jump to a further point
//             :jump_table[up][UP_SMALL]; // else, to a closer one.

//           // Debug updated up
//           // std::cout << "[DEBUG] Updated up with jump table: " << up << std::endl;
          
//           // Add safety check here
//           if (up >= old_size || jump_table[up][UP_BIG] >= old_size || jump_table[up][UP_SMALL] >= old_size)
//           {
//               up_stopped = true;
//               // std::cout << "[DEBUG] Stopping upward search due to jump table. " << std::endl;
//           }
//         } else { 
//           // If we are moving towards "start_cell", we can’t do any ot the previous optimizations and we just move to the next point.
//           up++;
//           // std::cout << "[DEBUG] Update up by ++: up = " << up << std::endl; 
//         }
//       } // if(now_up)

//       // repeat in down direction
//       // This is the specular part of the previous chunk of code.
//       else
//       {
//         if (down < 0)
//         {
//           down_stopped = true;
//           // std::cout << "[DEBUG] Stopping downward search: down = " << down << std::endl;
//           continue;
//         }

//         last_dist_down = trans_points[ind_trans].distToPoint2(&old_points[down]);
//         if (last_dist_down < min_dist) { best = down; min_dist = last_dist_down;}

//         if (down < start_index)
//         {
//           float theta_delta = trans_points[ind_trans].theta - old_points[down].theta;
//           float min_dist_down = sin(theta_delta) * trans_points[ind_trans].r;

//           if ((min_dist_down * min_dist_down) > min_dist)
//           {
//             down_stopped = true;
//             // std::cout << "[DEBUG] Early stopping downward search." << std::endl;
//             continue;
//           }

//           down = (old_points[down].r < trans_points[ind_trans].r) ? 
//                   jump_table[down][DOWN_BIG] : jump_table[down][DOWN_SMALL];

//           // Debug updated down
//           // std::cout << "[DEBUG] Updated down after jump table: " << down << std::endl;
          
//           // Add safety check here
//           if (down < 0 || jump_table[down][DOWN_BIG] < 0 || jump_table[down][DOWN_SMALL] < 0)
//           {
//               down_stopped = true;
//               // std::cout << "[DEBUG] Stopping downward search due to jump table." << std::endl;
//           }
//         }
//         else
//         {
//           down--;
//           // std::cout << "[DEBUG] Update down by ++: down = " << down << std::endl; 
//         }
//       }

//       // End of loop
//       // std::cout << "[DEBUG] End of loop: up = " << up << ", down = " << down
//       //           << ", up_stopped = " << up_stopped
//       //           << ", down_stopped = " << down_stopped << std::endl;
//       // std::cout << "[DEBUG] Best index of this loop: " << best 
//       //           << ", with distance: " << min_dist << std::endl;
//     }

//     // std::cout << "[DEBUG] Best old point: " << best 
//     //           << ", for index of trans point: " << ind_trans
//     //           << ", with distance: " << min_dist << std::endl;
//     // Find the second best point
//     int second_best = (best > 0) ? best - 1 : best + 1;

//     c.push_back(Correspondence(&trans_points[ind_trans], &points[ind_trans], &old_points[best], &old_points[second_best]));

//     // Update last best index
//     last_best = best;
//   }
//   // std::cout << "Fast correspond search done for scan group" << std::endl;
// }

void getCorrespondence(vector<Point>& old_points, vector<Point>& trans_points, vector<Point>& points,
                        vector< vector<int> >& jump_table, vector<Correspondence>& c, float prob)
{
  c.clear();
  int last_best = -1;
  const int n = trans_points.size();
  const int m = old_points.size();
  // cout<<"size n  "<< n << endl;

  //Do for each point
  for(int i = 0; i<n; ++i)
  {
    int start_at = (last_best<0)?(last_best+1):last_best;
    int min_index = start_at;
    int second_min_index;
    float min_dist = 10000.0;

    bool up_stop = false;
    bool down_stop = false;

    // start from last best, search upwards
    int current = start_at;

    while((!up_stop) && (current < m))
    {
      float dist = trans_points[i].distToPoint2(&old_points[current]); 
      if (dist < min_dist)
      {
        min_dist = dist;
    		min_index = current;

    		if(min_index==0)
        {
          second_min_index = min_index+1;
        } else {
          second_min_index = min_index-1;
        }
      }

      //check for early termination
      float theta_diff_up = old_points[current].theta - trans_points[i].theta;
      if (trans_points[i].r * sin(theta_diff_up) > min_dist)
      { 
    		up_stop = true;
    		//cout<<"up stopped:"<< endl;
    	}

      // advance based on jump table
    	if(old_points[current].r < trans_points[i].r)
      {
    		current = jump_table[current][UP_BIG];
    	} else {
    		current = jump_table[current][UP_SMALL];
    	}
      //	cout<<"current up:  "<< current << endl;
    }

    // back to initial starting point and search downwards
    current = start_at;

    while ((!down_stop) && (current >= 0))
    {
      float dist = trans_points[i].distToPoint2(&old_points[current]);
      if (dist < min_dist)
      {
    		min_dist = dist;
    		min_index = current;

    		if(min_index==0)
        {
          second_min_index = min_index+1;
        } else {
          second_min_index = min_index-1;
        }
      }

      //check for early termination
      float theta_diff_down = trans_points[i].theta - old_points[current].theta;
    	if (trans_points[i].r * sin(theta_diff_down) > min_dist)
      { 
    		down_stop = true;
    		//cout<<"down stopped:"<< endl;
    	}

      // advance based on jump table
    	if(old_points[current].r < trans_points[i].r)
      {
    		current = jump_table[current][DOWN_BIG];
    	} else {
    		current = jump_table[current][DOWN_SMALL];
    	}
    	//cout<<"current down:  "<< current << endl;
    }

    c.push_back(Correspondence(&trans_points[i], &points[i], &old_points[min_index], &old_points[second_min_index]));
    last_best = min_index;
    // cout<<"min_index:  "<< min_index << endl;
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
