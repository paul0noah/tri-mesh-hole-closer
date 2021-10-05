//
//  find_holes.hpp
//  tri-mesh-hole-closer
//
//  Created by Paul R on 03.10.21.
//

#ifndef find_holes_hpp
#define find_holes_hpp

#include <Eigen/Dense>
#include <vector>
#include <iostream>

namespace tmhc {

/* find_holes
    Finds all holes in a given mesh and extract the indices of the boundaries
    Input:
        - F: contains triangulation information of open mesh
        - E: contains the edges of the mesh
        - L: contains the mapping between both orientations of E and F
        - holeEdgeIndices: indices of edges of each hole => each hole has its own std::vector
    returns number of holes
 */
inline int find_holes(Eigen::MatrixXi F,
               Eigen::MatrixXi E,
               Eigen::MatrixXi L,
               std::vector<std::vector<int>> &holeEdgeIndices) {
    using namespace tmhc;
    const int numEdges = E.rows();

    // extract edges without inverse counterpart
    std::vector<int> edgesWithoutNeighbor;
    edgesWithoutNeighbor.reserve((L.array() == -1).count());
    for (int i = 0; i < numEdges; i++) {
        if (L(i, 1) == -1) {
            edgesWithoutNeighbor.push_back(i);
        }
    }

    // find all holes by traversing the border of each hole
    while (edgesWithoutNeighbor.size() > 0) {
        // always take the first element of edgesWithoutNeighbor
        std::vector<int> holeBoundary;
        holeBoundary.push_back(edgesWithoutNeighbor[0]);
        const int firstVertex = E(edgesWithoutNeighbor[0], 0);
        int currentVertex = E(edgesWithoutNeighbor[0], 1);
        edgesWithoutNeighbor.erase(edgesWithoutNeighbor.begin());
        bool loopClosed = false;
        while (!loopClosed) {
            for (int e = 0; e < edgesWithoutNeighbor.size(); e++) {
                const int nextVertex0 = E(edgesWithoutNeighbor[e], 0);
                const int nextVertex1 = E(edgesWithoutNeighbor[e], 1);
                if (currentVertex == nextVertex0) {
                    holeBoundary.push_back(edgesWithoutNeighbor[e]);
                    edgesWithoutNeighbor.erase(edgesWithoutNeighbor.begin() + e);
                    currentVertex = nextVertex1; // the other vertex of this edge
                    break;
                }
                if (currentVertex == nextVertex1) {
                    holeBoundary.push_back(edgesWithoutNeighbor[e]);
                    edgesWithoutNeighbor.erase(edgesWithoutNeighbor.begin() + e);
                    currentVertex = nextVertex0; // the other vertex of this edge
                    break;
                }
            }
            loopClosed = currentVertex == firstVertex;
        }
        holeEdgeIndices.push_back(holeBoundary);
    }

    return holeEdgeIndices.size();
}

} // namespace tmhc

#endif /* find_holes_hpp */
