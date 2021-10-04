//
//  close_holes.hpp
//  tri-mesh-hole-closer
//
//  Created by Paul R on 03.10.21.
//

#ifndef close_holes_hpp
#define close_holes_hpp

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include "find_holes.h"
#include "extract_edges.h"
#include "close_hole.h"

namespace tmhc {

/* close_holes
    Identifies all holes in a given mesh and closes them
    Input:
        - V: contains vertex positions of open mesh
        - F: contains triangulation information of open mesh
        - Vclosed: contains vertex positions of closed mesh
        - Fclosed: contains triangulation information of closed mesh
        - surfaceFaring: do surface fairing
    returns number of holes
 */
int close_holes(Eigen::MatrixXd V,
                Eigen::MatrixXi F,
                Eigen::MatrixXd &Vclosed,
                Eigen::MatrixXi &Fclosed,
                bool surfaceFaring) {
    assert(V.cols() == 3);
    assert(F.cols() == 3);
    assert(F.maxCoeff() < V.rows());

    using namespace tmhc;

    // extract edges to check watertightness
    Eigen::MatrixXi E, L;
    const bool isWaterTight = extract_edges(F, E, L);
    if (isWaterTight) {
        std::cout << "Mesh does not contain any holes" << std::endl;
        Vclosed = V;
        Fclosed = F;
        return 0;
    }

    std::vector<std::vector<int>> holeEdgeIndices;
    const int numHoles = find_holes(F, E, L, holeEdgeIndices);

    Eigen::MatrixXd Vwork = V;
    Eigen::MatrixXi Fwork = F;
    for (int hole = 0; hole < numHoles; hole++) {
        close_hole(Vwork, Fwork, E, holeEdgeIndices[hole], Vclosed, Fclosed);
        // we dont need an edge matrix to work one since we will find new holes in this step
        Vwork = Vclosed;
        Fwork = Fclosed;
    }

    if (surfaceFaring && Vclosed.rows() > V.rows()) {
        std::cout << "SurfaceFaring not yet implemented" << std::endl;
    }
    return numHoles;
}

int close_holes(Eigen::MatrixXd V,
                Eigen::MatrixXi F,
                Eigen::MatrixXd &Vclosed,
                Eigen::MatrixXi &Fclosed) {
    return close_holes(V, F, Vclosed, Fclosed, false);
}

} // namespace tmhc

#endif /* close_holes_hpp */
