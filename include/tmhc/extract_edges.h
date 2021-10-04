//
//  extract_edges.hpp
//  tri-mesh-hole-closer
//
//  Created by Paul R on 03.10.21.
//

#ifndef extract_edges_hpp
#define extract_edges_hpp

#include <Eigen/Dense>
#include <vector>
#include <iostream>

namespace tmhc {

inline bool allEqual(Eigen::MatrixXi inp1, Eigen::MatrixXi inp2) {
    assert(inp1.rows() == inp2.rows());
    assert(inp1.cols() == inp2.cols());

    return (inp1 - inp2).norm() == 0;
}

/* extract_edges
    Finds all holes in a given mesh and extract the indices of the boundaries
    Input:
        - F: contains triangulation information of open mesh
        - E: contains edges of mesh
        - L: maps edges and inverse orientation of edge to faces
             => any "-1" int the second column of L means that we did not find both orientations
    returns if mesh is watertight
 */
bool extract_edges(Eigen::MatrixXi F,
                   Eigen::MatrixXi &E,
                   Eigen::MatrixXi &L) {
    using namespace tmhc;

    const int numFaces = F.rows();
    const int numEdgesClosedMesh = numFaces * 3 / 2;
    const int maxNumEdges = numFaces * 3;
    E.resize(maxNumEdges, 2);
    L.resize(maxNumEdges, 2);
    L = - L.setOnes();

    Eigen::Vector2i idxEdge0; idxEdge0 << 0, 1;
    Eigen::Vector2i idxEdge1; idxEdge1 << 1, 2;
    Eigen::Vector2i idxEdge2; idxEdge2 << 2, 0;

    // we add the first edges before the loop => otherwise first iteration does not work
    E(0, Eigen::all) = F(0, idxEdge0);
    L(0, 0) = 0;
    E(1, Eigen::all) = F(0, idxEdge1);
    L(1, 0) = 0;
    E(2, Eigen::all) = F(0, idxEdge2);
    L(2, 0) = 0;
    int numEdgesAdded = 3;

    Eigen::Vector2i edge0, edge1, edge2;
    bool found0, found1, found2;
    Eigen::Vector2i minusEdge; minusEdge << 1, 0;

    // will be overwritten later
    bool watertight = true;

    for (int f = 1; f < numFaces; f++) {
        edge0 = F(f, idxEdge0); found0 = false;
        edge1 = F(f, idxEdge1); found1 = false;
        edge2 = F(f, idxEdge2); found2 = false;

        for (int e = 0; e < numEdgesAdded; e++) {
            // edge0
            if (allEqual( E(e, Eigen::all), edge0.transpose()) ) {
                found0 = true;
            }
            else if (allEqual( E(e, Eigen::all), edge0(minusEdge).transpose()) ) {
                L(e, 1) = f;
                found0 = true;
            }

            // edge1
            if (allEqual( E(e, Eigen::all), edge1.transpose()) ) {
                found1 = true;
            }
            else if (allEqual( E(e, Eigen::all), edge1(minusEdge).transpose()) ) {
                L(e, 1) = f;
                found1 = true;
            }

            // edge2
            if (allEqual( E(e, Eigen::all), edge2.transpose()) ) {
                found2 = true;
            }
            else if (allEqual( E(e, Eigen::all), edge2(minusEdge).transpose()) ) {
                L(e, 1) = f;
                found2 = true;
            }

        }

        try {
            // Add edges if any new were found
            if (!found0) {
                E(numEdgesAdded, Eigen::all) = edge0;
                L(numEdgesAdded, 0) = f;
                numEdgesAdded++;
            }
            if (!found1) {
                E(numEdgesAdded, Eigen::all) = edge1;
                L(numEdgesAdded, 0) = f;
                numEdgesAdded++;
            }
            if (!found2) {
                E(numEdgesAdded, Eigen::all) = edge2;
                L(numEdgesAdded, 0) = f;
                numEdgesAdded++;
            }
        } catch (const std::exception& e) {
            watertight = false;
        }
    }

    L = L.block(0, 0, numEdgesAdded, 2);
    E = E.block(0, 0, numEdgesAdded, 2);

    // Watertightness check
    if (numEdgesAdded != numEdgesClosedMesh) {
        watertight = false;
    }
    if ((L.array() == -1).any()) {
        watertight = false;
    }
    return watertight;
}

} // namespace tmhc

#endif /* extract_edges_hpp */
