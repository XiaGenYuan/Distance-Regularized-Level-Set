//
//  level_set.hpp
//  DRLevelSet
//
//  Created by Genyuan Xia on 2016/12/24.
//  Copyright © 2016年 Genyuan Xia. All rights reserved.
//

#ifndef level_set_hpp
#define level_set_hpp
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using std::vector;
/// LevelSet2D: 2D level set model
/// The class is implemented based on "Distrance Regularied Level Set Evolution and Its
/// Application to Image Segmentation" by Chunming Li.

enum PotentialFunction {
    SINGLE_WELL = 0,
    DOUBLE_WELL = 1
};

class LevelSet2D {
public:
    // constructor with parameters
    LevelSet2D(std::vector<double> &image_data, int image_width, int image_height,
                     double mu, double timestep, double lambda, double alfa, double epsilon,
                     PotentialFunction potential_function, int iters);
    // destructor
    ~LevelSet2D();
    // initialization
    void Initialization();
    
public:
    // init level set function
    void InitPhi(std::vector<bool>& mask, const int& image_width, const int& image_height);
    // compute gradient in x and y direction
    void Gradient(std::vector<double>& g, std::vector<double>& vx, std::vector<double>& vy,
                  const int& image_width, const int& image_height);
    // level set evolution
    void Evolution();
    // make a function satisfy Neumann boundary condition
    void NeumannBoundCond(vector<double>& phi, const int& image_width, const int& image_height);
    // compute del2 function
    void Del2(vector<double>& phi, vector<double>& del2, const int& image_width,
              const int& image_height);
    // compute Dirac function
    void Dirac(vector<double>& phi, vector<double>& diracPhi, const int& image_width,
               const int& image_height);
    // get phi function
    std::vector<double>& GetPhi();
    // Gauss low-pass filter
    void GaussLowpassFilter(vector<double>& image_data, const int& image_width,
                            const int& image_height);
    // compute edge indicator function
    void EdgeIndicator(vector<double>& image_data, const int& image_width, const int& image_height,
                       vector<double>& edge_indicator);
    
private:
    // image data
    vector<double> image_data_;
    // image width
    int image_width_;
    // image height
    int image_height_;
    // level set function to be updated by level set evolution
    vector<double> phi_;
    // edge indicator function
    vector<double> g_;
    // distance regularization term
    vector<double> dist_reg_term_;
    // Dirac phi function
    vector<double> dirac_phi_;
    // area term
    vector<double> area_term_;
    // edge term
    vector<double> edge_term_;
    // widget to distance regularization term
    double mu_;
    // time step for updating
    double timestep_;
    // weight of the weighted length term
    double lambda_;
    // weight of the weighted area term
    double alfa_;
    // width of Dirac Delta function
    double epsilon_;
    // choice of potential function in distance regularization term.
    // single-well is good for region-based model
    // double-well is good for both edge and region based models.
    PotentialFunction potential_function_;
    // number of iteration in one evolution
    int iters_;
};

#endif /* level_set_hpp */
