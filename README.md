# tri-mesh-hole-closer
Are you searching for a library which is able to close holes of triangular meshes with nice triangulation while only using `libigl` and `eigen`? Here you are!

This library closes holes of triangular meshes as follows
1) Finds holes in your mesh
2) Projects each boundary of a hole to a best fitting 2D plane with help of [SVD or rather PCA]( https://stats.stackexchange.com/questions/172300/low-rank-svd-reconstruction-and-linear-projection-of-the-data)
3) Uses `libigl::triangulate` function to obtain new points and nice triangulation in the inside of the hole
4) Projects the data back up in 3D space
5) Performs surface fairing inspired by [this article](https://erkaman.github.io/posts/hole_filling.html) to smooth out the new triangulation within the holes (this is not yet implemented)

## Usage

1) Add the library to your C++ project by adding the following lines to your CMake
    ```
    add_subdirectory(path/to/tri-mesh-hole-closer)
    target_link_libraries(your_project tri-mesh-hole-closer)
    ```
    if you use libigl already within your project add the following line to your CMake
    ```
    option(LIBIGL_WITH_TRIANGLE "Use Triangle" ON)
    ```
2) Include the following header
    ```
    #include <tmhc/close_holes.h>
    ```
3) Call as follows
    ```
    Eigen::MatrixXi F = ...; // should contain triangulation of your mesh with holes
    Eigen::MatrixXd V = ...; // should contain vertices of your mesh with holes
    Eigen::MatrixXi Fclosed; // will contain triangulation of closed mesh
    Eigen::MatrixXd Vclosed; // will contain vertices of closed mesh
    bool do_surface_fairing = true;
    int num_holes = tmhc::close_holes(V, F, Vclosed, Fclosed, do_surface_fairing);
    ```

## Dependencies

We need the following libraries
| Lib | URL |
| ------ | ------ |
| libigl | `https://github.com/libigl/libigl` | 
| Eigen3 | `https://gitlab.com/libeigen/eigen` |

## License

This project can be used with the MIT License

