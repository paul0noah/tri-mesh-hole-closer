//
//  surface_fairing.hpp
//  tri-mesh-hole-closer
//
//  Created by Paul R on 03.10.21.
//

#ifndef surface_fairing_hpp
#define surface_fairing_hpp

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <igl/harmonic.h>

namespace tmhc {

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
inline void surface_fairing(Eigen::MatrixXd &Vclosed,
                            Eigen::MatrixXi Fclosed,
                            int numOldPoints,
                            std::vector<std::vector<int>> holeEdgeIndices) {

    if (numOldPoints > 2000) {
        std::cout << "surface_fairing not yet implemented efficiently" << std::endl;
        return;
    }
    Eigen::MatrixXd fairedV = Vclosed;
    Eigen::MatrixXi fairedF = Fclosed;
    // now we shall do surface fairing on the mesh, to ensure
    // that the patch conforms to the surrounding curvature.
    {

        Eigen::VectorXi b(numOldPoints);
        Eigen::MatrixXd bc(numOldPoints, 3);
        // setup the boundary conditions. This is simply the vertex positions of the vertices not part of the patch.
        // since we add all patch indices to the end of our matrix V we just add the old points to the boundary condition
        for (int i = 0; i < numOldPoints; i++) {
            b(i) = i;

            bc(i, 0) = fairedV(i, 0);
            bc(i, 1) = fairedV(i, 1);
            bc(i, 2) = fairedV(i, 2);
        }


        Eigen::MatrixXd Z;
        int k = 2;
        // surface fairing simply means that we solve the equation
        // Delta^2 f = 0
        // with appropriate boundary conditions.
        // this function igl::harmonic from libigl takes care of that.

        // note that this is pretty inefficient thought.
        // the only boundary conditions necessary are the 2-ring of vertices around the patch.
        // the rest of the mesh vertices need not be specified
        igl::harmonic(fairedV, fairedF, b, bc, k, Z);
        if (Z.rows() == Vclosed.rows())
            Vclosed = Z;
        else
            std::cout << "Surface fairing was not successfull" << std::endl;
    }
}

} // namespace tmhc

#endif /* surface_fairing_hpp */
