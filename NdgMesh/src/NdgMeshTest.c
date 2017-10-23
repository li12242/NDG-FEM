#include "NdgMesh.h"
#include "NdgMeshIO.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    NdgUMeshUnion *mesh = mxGetNdgUMeshUnion(prhs[0]);
    //mxPrintUMeshUnion(mesh, "UMesh");
    NdgNcFile *ncfile = setupBasicUMeshUnionOutputFile(mesh[0], "TestMesh", 1);
    closeNetcdfFile(ncfile[0]);
    freeNdgUMeshUnion(mesh);
    return;
}