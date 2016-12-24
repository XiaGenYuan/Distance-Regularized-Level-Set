//
//  main.cpp
//  DRLevelSet
//
//  Created by Genyuan Xia on 2016/12/24.
//  Copyright © 2016年 Genyuan Xia. All rights reserved.
//

#include <iostream>
#include "level_set.hpp"

int main(int argc, const char * argv[]) {
    /* 
     *  This is one demo how to use DeformLevelSet2D class
     */
    int image_width = 128;
    int image_height = 128;
    int image_size = image_width * image_height;
    // TODO: input your image data
    vector<double> image_data(image_size, 0);
    double timestep = 5;
    double mu = 0.2 / timestep;
    double lambda = 20;
    double alfa = 4;
    double epsilon = 1.5;
    int iters = 40;
    PotentialFunction potential_function = DOUBLE_WELL;
    LevelSet2D level_set2D(image_data, image_width, image_height, mu, timestep, lambda, alfa,
                                 epsilon, potential_function, iters);
    // TODO: init your segmentation mask
    vector<bool> mask(image_size, false);
    level_set2D.InitPhi(mask, image_width, image_height);
    for(int i = 0; i < 5; ++ i) {
        level_set2D.Evolution();
    }
    alfa = 0;
    level_set2D.Evolution();
    vector<double>& phi = level_set2D.GetPhi();
    // use updated phi function to update the segmentation result
    return 0;
}
