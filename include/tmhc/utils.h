//
//  find_holes.hpp
//  tri-mesh-hole-closer
//
//  Created by Paul R on 03.10.21.
//

#ifndef tmhc_utils_hpp
#define tmhc_utils_hpp

#include <Eigen/Dense>
#include <vector>
#include <iostream>

namespace tmhc {
namespace utils {

inline bool allEqual(Eigen::MatrixXi inp1, Eigen::MatrixXi inp2) {
    assert(inp1.rows() == inp2.rows());
    assert(inp1.cols() == inp2.cols());

    return (inp1 - inp2).norm() == 0;
}

const float getTriangleArea(Eigen::MatrixXi triangle, Eigen::MatrixXf V) {
    // get vertices of the triangle
    Eigen::Vector3f v0 = V.row(triangle(0));
    Eigen::Vector3f v1 = V.row(triangle(1));
    Eigen::Vector3f v2 = V.row(triangle(2));

    // extract two edges
    Eigen::Vector3f e1 = v0 - v1;
    Eigen::Vector3f e2 = v0 - v2;

    // triangle area via cross product
    return 0.5 * e1.cross(e2).norm();
}

Eigen::MatrixXf getTriangleAreas(Eigen::MatrixXi F, Eigen::MatrixXf V) {
    Eigen::MatrixXf areas(F.rows(), 1);
    for (int i = 0; i < F.rows(); i++) {
        areas(i) = getTriangleArea(F(i, Eigen::all), V);
    }
    return areas;
}

Eigen::MatrixXi linspaced(int start, int end, int step) {
    assert(step > 0);
    assert(end > start);
    int length = (end - start)/step;
    Eigen::MatrixXi A(length, 1);
    for (int i = 0; i < length; i++) {
        A(i, 0) = start + i * step;
    }
    return A;
}

Eigen::MatrixXi linspaced(int start, int end) {
    return linspaced(start, end, 1);
}

} // namespace utils
} // namespace tmhc

#endif /* tmhc_utils_hpp */
