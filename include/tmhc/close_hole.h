//
//  close_hole.hpp
//  tri-mesh-hole-closer
//
//  Created by Paul R on 03.10.21.
//

#ifndef close_hole_hpp
#define close_hole_hpp

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <igl/triangle/triangulate.h>

#include "extract_edges.h"

namespace tmhc {

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

/* close_hole
    closes a hole in a triangle mesh using libigl triangle function and svd
    Input:
        - V: contains vertex positions of open mesh
        - F: contains triangulation information of open mesh
        - E: contains the edges of the mesh
        - Vclosed: contains vertex positions of closed mesh
        - Fclosed: contains triangulation information of closed mesh
        - surfaceFaring: do surface fairing
    returns true if the hole was completely closed
 */
bool close_hole(Eigen::MatrixXd V,
                Eigen::MatrixXi F,
                Eigen::MatrixXi E,
                std::vector<int> holeBoundaryEdges,
                Eigen::MatrixXd &Vclosed,
                Eigen::MatrixXi &Fclosed) {


    // find plane which best fits all boundary points via svd
    Eigen::MatrixXf BoundaryPoints(holeBoundaryEdges.size(), 3);
    Eigen::MatrixXi boundEdgeToVertMap(holeBoundaryEdges.size(), 1);
    for (int i = 0; i < holeBoundaryEdges.size(); i++) {
        const int vertexIdx = E(holeBoundaryEdges[i], 0);
        boundEdgeToVertMap(i) = vertexIdx;
        BoundaryPoints.row(i) = V.row(vertexIdx).cast<float>();
    }
    Eigen::Matrix<float, 1, 3> meanBoundaryPoints = BoundaryPoints.colwise().mean();
    const Eigen::MatrixX3f BoundaryPointsCentered = BoundaryPoints.rowwise() - meanBoundaryPoints;

    //TODO: use other eigen SVD to be able to handle more points
    Eigen::JacobiSVD<Eigen::MatrixX3f> svd = BoundaryPointsCentered.jacobiSvd(Eigen::ComputeFullV);

    Eigen::MatrixXf PlaneProjection = svd.matrixV().block(0, 0, 3, 2);


    // Project Boundary points to best fitting 2D Plane we found before
    Eigen::MatrixXd VBound2D = (BoundaryPointsCentered * PlaneProjection).cast<double>();
    // Connect all boundary points in EBound
    Eigen::MatrixXi EBound(VBound2D.rows(), 2);
    for (int i = 0; i < VBound2D.rows(); i++) {
        EBound(i, 0) = i;
        EBound(i, 1) = i == VBound2D.rows()-1 ? 0 : i+1;
    }

    // specify a point that is inside a closed shape
    // where we do not want triangulation to happen
    Eigen::MatrixXd H;

    // Patch matrices
    Eigen::MatrixXd Vpatch;
    Eigen::MatrixXi Fpatch;


    const float meanTriangleArea = getTriangleAreas(F, V.cast<float>()).mean();

    // Triangulate the interior of Vbound2D connected with Ebound
    //  - a0.005 means that the area of each triangle should
    //    not be greater than 0.005
    //  - q15 means that no angles will be smaller than 15 degrees
    //  - Y means no new points at the boundary
    // for a detailed set of commands please refer to:
    // https://www.cs.cmu.edu/~quake/triangle.switch.html
    std::string triCmd = "q15a" + std::to_string(2*meanTriangleArea) + "Y";
    igl::triangle::triangulate(VBound2D, EBound, H, triCmd, Vpatch, Fpatch);
    Eigen::MatrixXf Vpatchf = Vpatch.cast<float>();

    // project the points up in 3D space (we will replace the boundary with their originals later)
    Vpatchf = Vpatchf * PlaneProjection.transpose();
    Vpatchf = Vpatchf.rowwise() + meanBoundaryPoints;

    // translate entries in Fpatch for the full triangulation matrix
    Eigen::MatrixXi boundEToVbigMap(Vpatchf.rows(), 1);
    const int numBoundaryPoints = BoundaryPoints.rows();
    const int numNewPoints = Vpatchf.rows() - numBoundaryPoints;
    const int numOldPoints = V.rows();
    boundEToVbigMap.block(0, 0, numBoundaryPoints, 1) = boundEdgeToVertMap;
    if (numNewPoints) {
        boundEToVbigMap.block(numBoundaryPoints, 0, numNewPoints, 1) = utils::linspaced(numOldPoints, numOldPoints + numNewPoints);
    }
    for (int i = 0; i < Fpatch.rows(); i++) {
        for (int j = 0; j < Fpatch.cols(); j++) {
            Fpatch(i, j) = boundEToVbigMap(Fpatch(i, j));
        }
    }
    // assemble Fbig (fuse mesh with the patch)
    const int oldNumFaces = F.rows();
    const int numNewFaces = Fpatch.rows();
    Eigen::MatrixXi Fbig(oldNumFaces + numNewFaces, 3);
    Fbig.block(0, 0, oldNumFaces, 3) = F;
    Fbig.block(oldNumFaces, 0, numNewFaces, 3) = Fpatch;

    // assemble Vbig (fuse mesh with the patch)
    Eigen::MatrixXf Vbig(numOldPoints + numNewPoints, 3);
    Vbig.block(0, 0,  numOldPoints, 3) = V.cast<float>();
    Vbig.block(numOldPoints, 0, numNewPoints, 3) = Vpatchf.block(numBoundaryPoints, 0, numNewPoints, 3);


    // check if we closed the hole successfully and if the triangles are correctly oriented
    Eigen::MatrixXi Epatch, Lpatch;
    extract_edges(Fpatch, Epatch, Lpatch);
    bool completelyClosed = true;
    bool correctOrientation = true;
    for (int e = 0; e < holeBoundaryEdges.size(); e++) {
        const int edgeIdx = holeBoundaryEdges[e];
        bool foundEdge = false;
        for (int i = 0; i < Epatch.rows(); i++) {
            // for a closed hole we need to find a counterpart for each boundary edge within the patch
            if (E(edgeIdx, 0) - Epatch(i, 0) == 0 &&
                 E(edgeIdx, 1) - Epatch(i, 1) == 0) {
                correctOrientation = false;
                foundEdge = true;
                break;
            }
            if (E(edgeIdx, 1) - Epatch(i, 0) == 0 &&
                E(edgeIdx, 0) - Epatch(i, 1) == 0) {
                foundEdge = true;
                break;
            }
        }
        if (!foundEdge) {
            completelyClosed = false;
        }
    }
    if (!correctOrientation) {
        // if the triangles are not correctly oriented we have to flip their orientation
        Fbig.block(oldNumFaces, 1, Fpatch.rows(), 1) = Fpatch.col(2);
        Fbig.block(oldNumFaces, 2, Fpatch.rows(), 1) = Fpatch.col(1);
    }

    Fclosed = Fbig;
    Vclosed = Vbig.cast<double>();
    return completelyClosed;
}

} // namespace tmhc

#endif /* close_hole_hpp */
