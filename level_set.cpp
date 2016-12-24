//
//  level_set.cpp
//  DRLevelSet
//
//  Created by Genyuan Xia on 2016/12/24.
//  Copyright © 2016年 Genyuan Xia. All rights reserved.
//

#include "level_set.hpp"

// constructor with parameters
LevelSet2D::LevelSet2D(std::vector<double> &image_data, int image_width,
                                   int image_height, double mu, double timestep, double lambda, double alfa, double epsilon,
                                   PotentialFunction potential_function, int iters) {
    this->image_data_ = image_data;
    this->image_width_ = image_width;
    this->image_height_ = image_height;
    this->mu_ = mu;
    this->timestep_ = timestep;
    this->lambda_ = lambda;
    this->alfa_ = alfa;
    this->epsilon_ = epsilon;
    this->potential_function_ = potential_function;
    this->iters_ = iters;
    this->Initialization();
}

// destructor
LevelSet2D::~LevelSet2D(){
    
}

// initialization
void LevelSet2D::Initialization() {
    int image_size = this->image_width_ * this->image_height_;
    this->phi_.resize(image_size, 0.0);
    this->g_.resize(image_size, 0.0);
    this->dist_reg_term_.resize(image_size, 0.0);
    this->dirac_phi_.resize(image_size, 0.0);
    this->area_term_.resize(image_size, 0.0);
    this->edge_term_.resize(image_size, 0.0);
    
    this->GaussLowpassFilter(this->image_data_, image_width_, image_height_);
    this->EdgeIndicator(this->image_data_, image_width_, image_height_, this->g_);
}

// init phi function
void LevelSet2D::InitPhi(std::vector<bool>& mask, const int& image_width,
                               const int& image_height){
    double c0 = 2;
    int image_size = image_width * image_height;
    for(int j = 0; j < image_size; ++ j) {
        if(mask[j]) {
            this->phi_[j] = -c0;
        } else {
            this->phi_[j] = c0;
        }
    }
}

// Gauss low-pass filter
void LevelSet2D::GaussLowpassFilter(vector<double>& image_data, const int& image_width,
                                          const int& image_height){
    // define gauss kernel
    const int kernel_size = 9;
    const int kernel_center = kernel_size / 2;
    const double sigma = 0.8;
    double gauss_kernel[kernel_size][kernel_size];
    double sumg = 0.0;
    for(int i = 0; i < kernel_size; ++ i) {
        for(int j = 0; j < kernel_size; ++ j) {
            gauss_kernel[i][j] = exp(-(std::pow((i - kernel_center), 2.0) +
                                       std::pow((j - kernel_center), 2.0)) / (2 * std::pow(sigma, 2.0)));
            sumg += gauss_kernel[i][j];
        }
    }
    for(int i = 0; i < kernel_size; ++ i) {
        for(int j = 0; j < kernel_size; ++ j) {
            gauss_kernel[i][j] /= sumg;
        }
    }
    int image_index, center_image_index;
    double kernel_result = 0.0;
    for(int x = 0; x < image_width - kernel_size; ++ x) {
        for(int y = 0; y < image_height - kernel_size; ++ y) {
            kernel_result = 0.0;
            center_image_index = (y + kernel_center) * image_width + x + kernel_center;
            for(int i = 0; i < kernel_size; ++ i) {
                for(int j = 0; j < kernel_size; ++ j) {
                    image_index = (y + j) * image_width + x + i;
                    kernel_result += image_data[image_index] * gauss_kernel[j][i];
                }
            }
            image_data[center_image_index] = kernel_result;
        }
    }
}

// compute edge indicator function
void LevelSet2D::EdgeIndicator(vector<double>& image_data, const int& image_width,
                                     const int& image_height, vector<double>& edge_indicator) {
    int image_size = image_width * image_height;
    vector<double> ix(image_size), iy(image_size);
    this->Gradient(image_data, ix, iy, image_width, image_height);
    for (int i = 0; i < image_size; ++i){
        edge_indicator[i] = 1.0 / (1.0 + ix[i] * ix[i] + iy[i] * iy[i]);
    }
}

// compute gradient in x and y direction
void LevelSet2D::Gradient(std::vector<double>& g, std::vector<double>& vx,
                                std::vector<double>& vy, const int& image_width, const int& image_height) {
    int image_size = image_width * image_height;
    if (vx.size() != image_size) {
        vx.resize(image_size, 0.0);
    }
    if (vy.size() != image_size) {
        vy.resize(image_size, 0.0);
    }
    int cur_index, left_index, right_index, up_index, down_index;
    for(int x = 1; x < image_width - 1; ++ x) {
        for(int y = 1; y < image_height - 1; ++ y) {
            cur_index = y * image_width + x;
            left_index = y * image_width + x - 1;
            right_index = y * image_width + x + 1;
            up_index = (y - 1) * image_width + x;
            down_index = (y + 1) * image_width + x;
            vx[cur_index] = (g[right_index] - g[left_index]) / 2.0;
            vy[cur_index] = (g[down_index] - g[up_index]) / 2.0;
        }
    }
}

// make a function satisfy Neumann boundary condition
void LevelSet2D::NeumannBoundCond(vector<double>& phi, const int& image_width,
                                        const int& image_height) {
    int nrow = image_height;
    int ncol = image_width;
    vector<double> g(phi);
    phi[0 * ncol + 0] = g[2 * ncol + 2];
    phi[0 * ncol + ncol - 1] = g[2 * ncol + ncol - 3];
    phi[(nrow - 1) * ncol + 0] = g[(nrow - 3) * ncol + 2];
    phi[(nrow - 1) * ncol + ncol - 1] = g[(nrow - 3) * ncol + ncol - 3];
    g = phi;
    for(int y = 1; y < ncol - 1; ++ y) {
        phi[0 * ncol + y] = g[2 * ncol + y];
        phi[(nrow - 1) * ncol + y] = g[(nrow - 3) * ncol + y];
    }
    g = phi;
    for(int x = 1; x < nrow - 1; ++ x) {
        phi[x * ncol + 0] = g[x * ncol + 2];
        phi[x * ncol + ncol - 1] = g[x * ncol + ncol - 3];
    }
}

// compute del2 fuction
void LevelSet2D::Del2(vector<double>& phi, vector<double>& del2, const int& image_width,
                            const int& image_height) {
    del2.resize(image_width * image_height, 0);
    int cur_index, left_index, right_index, up_index, down_index;
    for(int x = 1; x < image_width - 1; ++ x) {
        for(int y = 1; y < image_height - 1; ++ y) {
            cur_index = y * image_width + x;
            left_index = y * image_width + x - 1;
            right_index = y * image_width + x + 1;
            up_index = (y - 1) * image_width + x;
            down_index = (y + 1) * image_width + x;
            del2[cur_index] = (phi[up_index] + phi[down_index] + phi[left_index] +
                               phi[right_index]) / 4.0 - phi[cur_index];
        }
    }
}

// compute Dirac function
void LevelSet2D::Dirac(vector<double>& phi, vector<double>& dirac_phi, const int& image_width,
                             const int& image_height){
    int image_size = image_width * image_height;
    const double PI = 3.1415926;
    for(int j = 0; j < image_size; ++ j) {
        dirac_phi[j] = (1.0 / 2.0 / this->epsilon_) * (1 + std::cos(PI * phi[j] / this->epsilon_));
        if(std::fabs(phi[j]) > std::fabs(this->epsilon_)){
            dirac_phi[j] = 0;
        }
    }
}

// level set evolution
void LevelSet2D::Evolution(){
    const double PI = 3.1415926;
    const double small_number = 1e-10;
    int image_size = this->image_width_ * this->image_height_;
    vector<double> vx(image_size), vy(image_size);
    vector<double> phi_x(image_size), phi_y(image_size), s(image_size);
    vector<double> nx(image_size), ny(image_size);
    vector<double> nxx(image_size), nyy(image_size), junk(image_size), curvature(image_size);
    vector<double> del2(image_size);
    vector<double> ps(image_size), dps(image_size);
    vector<double> nxx_p2(image_size), nyy_p2(image_size);
    
    this->Gradient(this->g_, vx, vy, this->image_width_, this->image_height_);
    for(int i = 0; i < this->iters_; ++ i) {
        this->NeumannBoundCond(this->phi_, this->image_width_, this->image_height_);
        this->Gradient(this->phi_, phi_x, phi_y, this->image_width_, this->image_height_);
        
        for(int j = 0; j < image_size; ++ j) {
            s[j] = std::sqrt(phi_x[j] * phi_x[j] + phi_y[j] * phi_y[j]);
            nx[j] = phi_x[j] / (small_number + s[j]);
            ny[j] = phi_y[j] / (small_number + s[j]);
        }
        
        this->Gradient(nx, nxx, junk, this->image_width_, this->image_height_);
        this->Gradient(ny, junk, nyy, this->image_width_, this->image_height_);
        for(int j = 0; j < image_size; ++ j) {
            curvature[j] = nxx[j] + nyy[j];
        }
        this->Del2(this->phi_, del2, this->image_width_, this->image_height_);
        if(SINGLE_WELL == this->potential_function_) {
            for(int j = 0; j < image_size; ++ j) {
                this->dist_reg_term_[j] = 4 * del2[j] - curvature[j];
            }
        } else {
            for(int j = 0; j < image_size; ++ j) {
                if(s[j] <= 1) {
                    ps[j] = std::sin(2 * PI * s[j]) / (2 * PI);
                } else {
                    ps[j] = s[j] - 1;
                }
            }
            for(int j  = 0; j < image_size; ++ j) {
                if(0 == ps[j] && 0 == s[j]) {
                    dps[j] = 1.0;
                } else if(0 == ps[j]){
                    dps[j] = 1.0 / s[j];
                } else if (0 == s[j]) {
                    dps[j] = ps[j] / 1.0;
                } else {
                    dps[j] = ps[j] / s[j];
                }
            }
            for(int j = 0; j < image_size; ++ j) {
                nxx_p2[j] = dps[j] * phi_x[j] - phi_x[j];
                nyy_p2[j] = dps[j] * phi_y[j] - phi_y[j];
            }
            this->Gradient(nxx_p2, nxx, junk, this->image_width_, this->image_height_);
            this->Gradient(nyy_p2, junk, nyy, this->image_width_, this->image_height_);
            for(int j = 0; j < image_size; ++ j) {
                this->dist_reg_term_[j] = nxx[j] + nyy[j] + 4 * del2[j];
            }
        }
        
        this->Dirac(this->phi_, this->dirac_phi_, this->image_width_, this->image_height_);
        for(int j = 0; j < image_size; ++ j) {
            this->area_term_[j] = this->dirac_phi_[j] * this->g_[j];
        }
        
        for(int j = 0; j < image_size; ++ j) {
            this->edge_term_[j] = this->dirac_phi_[j] * (vx[j] * nx[j] + vy[j] * ny[j]) +
            this->dirac_phi_[j] * this->g_[j] * curvature[j];
        }
        
        for(int j = 0; j < image_size; ++ j) {
            this->phi_[j] += this->timestep_ * (this->mu_ * this->dist_reg_term_[j] +
                                                this->lambda_ * this->edge_term_[j] + this->alfa_ * this->area_term_[j]);
        }
    }
}

// get phi function
std::vector<double>& LevelSet2D::GetPhi(){
    return this->phi_;
}