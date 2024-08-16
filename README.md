# 3D Mesh Volume
**Goal**: to measure the cell volume from 3D segmented cells in two different methods, and compare the differences. One is the **exact cell volume**, and the other is an approximation where cells are assumed to be straight, columnar (with the same area across the height of the cell) and volume is measured as the **product of the area and height**.


**Input**: This `python script` takes as input a `binarymeshformat` .bmf file comprised of triangular meshes that are created to segment cells in a tissue. This file is generated from ImageJ plug-in  `DM3D-https://github.com/FrancisCrickInstitute/dm3d-pages/blob/gh-pages/index.md` 

**Python VEDO Library**: It uses `VEDO-https://vedo.embl.es` library to extract the points and triangles forming the mesh object, and measures the volume of the mesh by counting pixels inside the mesh. The function `Fillmeshes()` measures the volume inside a mesh, as well as creates a binary image with the pixels inside marked with value 1. It then creates a height map by projecting the 3D stack on a single z-plane and adding them, such that the value at each pixel gives the height of the mesh at that pixel. The function  `Heightmap() ` takes this heightmap as input, selects one 2D plane where the cell mesh is present, to represent the apical plane of the cell, and calculates the value of the height at each pixel inside the apical boundary of the cell on that plane.

**Output**: It outputs two `csv` files with the coordinates of the centroid of the meshes, and the mesh volumes calculated using the two techniques respectively.
