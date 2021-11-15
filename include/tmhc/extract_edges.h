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
#include <tsl/robin_set.h>

#include "utils.h"


struct TMHC_EDG {
    int idx0;
    int idx1;
    int e; // index of edge
    TMHC_EDG () {}
    TMHC_EDG (int iidx0, int iidx1, int ie) {
        idx0 = iidx0;
        idx1 = iidx1;
        e = ie;
    }
    TMHC_EDG (Eigen::MatrixXi edge, int ie) {
        idx0 = edge(0);
        idx1 = edge(1);
        e = ie;
    }
    TMHC_EDG (Eigen::MatrixXi edge) {
        idx0 = edge(0);
        idx1 = edge(1);
        e = -1;
    }
    TMHC_EDG operator-() const {
        TMHC_EDG minusEDG;
        minusEDG.idx0 = idx1;
        minusEDG.idx1 = idx0;
        return minusEDG;
    }
    bool operator==(const TMHC_EDG& edg) const {
        return (idx0 == edg.idx0) && (idx1 == edg.idx1);
    }
};
namespace std {
    template<> struct hash<TMHC_EDG> {
        std::size_t operator()(TMHC_EDG const& edg) const noexcept {
            //size_t idx0hash = std::hash<int>()(edg.idx0);
            //size_t idx1hash = std::hash<int>()(edg.idx1) << 1;
            //return idx0hash ^ idx1hash;
            int k1 = edg.idx0;
            int k2 = edg.idx1;
            return (k1 + k2 ) * (k1 + k2 + 1) / 2 + k2;
        }
    };
    template<> struct equal_to<TMHC_EDG>{
        constexpr bool operator()(const TMHC_EDG &lhs, const TMHC_EDG &rhs) const {
            return (lhs.idx0 == rhs.idx0) && (lhs.idx1 == rhs.idx1);
        }
    };
}

namespace tmhc {

int findEdge(const tsl::robin_set<TMHC_EDG> &ELookup, const TMHC_EDG &edg) {
    auto it = ELookup.find(edg);
    if(it != ELookup.end()) {
        TMHC_EDG foundEdg = *it;
        return foundEdg.e;
    }
    return -1;
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
inline bool extract_edges(Eigen::MatrixXi F,
                   Eigen::MatrixXi &E,
                   Eigen::MatrixXi &L) {
    using namespace tmhc;

    const int numFaces = F.rows();
    const int numEdgesClosedMesh = numFaces * 3 / 2;
    const int maxNumEdges = numFaces * 3;
    Eigen::MatrixXi Ework(maxNumEdges, 2);
    Ework = -Ework.setOnes();
    Eigen::MatrixXi Lwork(maxNumEdges, 2);
    Lwork.setOnes();
    Lwork = - Lwork;
    tsl::robin_set<TMHC_EDG> ELookup; ELookup.reserve(maxNumEdges);

    Eigen::Vector2i idxEdge0; idxEdge0 << 0, 1;
    Eigen::Vector2i idxEdge1; idxEdge1 << 1, 2;
    Eigen::Vector2i idxEdge2; idxEdge2 << 2, 0;

    // we add the first edges before the loop => otherwise first iteration does not work
    Ework(0, Eigen::all) = F(0, idxEdge0);
    ELookup.insert(TMHC_EDG(Ework.row(0), 0));
    Lwork(0, 0) = 0;
    Ework(1, Eigen::all) = F(0, idxEdge1);
    ELookup.insert(TMHC_EDG(Ework.row(1), 1));
    Lwork(1, 0) = 0;
    Ework(2, Eigen::all) = F(0, idxEdge2);
    ELookup.insert(TMHC_EDG(Ework.row(2), 2));
    Lwork(2, 0) = 0;
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

        if (findEdge(ELookup, TMHC_EDG(edge0)) != -1) {
            found0 = true;
        }
        else {
            int e = findEdge(ELookup, -TMHC_EDG(edge0));
            if (e != -1) {
                found0 = true;
                Lwork(e, 1) = f;
            }
        }

        if (findEdge(ELookup, TMHC_EDG(edge1)) != -1) {
            found1 = true;
        }
        else {
            int e = findEdge(ELookup, -TMHC_EDG(edge1));
            if (e != -1) {
                found1 = true;
                Lwork(e, 1) = f;
            }
        }

        if (findEdge(ELookup, TMHC_EDG(edge2)) != -1) {
            found2 = true;
        }
        else {
            int e = findEdge(ELookup, -TMHC_EDG(edge2));
            if (e != -1) {
                found2 = true;
                Lwork(e, 1) = f;
            }
        }

        // Add edges if any new were found
        if (!found0) {
            Ework(numEdgesAdded, Eigen::all) = edge0;
            ELookup.insert(TMHC_EDG(edge0, numEdgesAdded));
            Lwork(numEdgesAdded, 0) = f;
            numEdgesAdded++;

        }
        if (!found1) {
            Ework(numEdgesAdded, Eigen::all) = edge1;
            ELookup.insert(TMHC_EDG(edge1, numEdgesAdded));
            Lwork(numEdgesAdded, 0) = f;
            numEdgesAdded++;

        }
        if (!found2) {
            Ework(numEdgesAdded, Eigen::all) = edge2;
            ELookup.insert(TMHC_EDG(edge2, numEdgesAdded));
            Lwork(numEdgesAdded, 0) = f;
            numEdgesAdded++;
        }
    }

    L = Lwork.block(0, 0, numEdgesAdded, 2);
    E = Ework.block(0, 0, numEdgesAdded, 2);

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
