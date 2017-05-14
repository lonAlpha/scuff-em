vim ./args_mesh
rm -rf ./cube_geo_bgmesh.hdf5
rm -rf ./cube_geo_bgmesh.pos
rm -rf ./Cube.msh
cp ./Cube_orig.msh Cube.msh 
../../src/applications/generate-bgmesh/generate-bgmesh <args_mesh
gmsh -2 Cube.geo -bgm ./cube_geo_bgmesh.pos
gmsh Cube.msh  ./cube_geo_MeshIndicator.pp
