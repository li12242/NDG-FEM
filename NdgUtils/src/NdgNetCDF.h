#ifndef NDGNETCDF_H
#define NDGNETCDF_H

#include "netcdf.h"
#define NC_MAX_NAME_LEN 120
#define NC_VAR_MAX_DIM 5

/** structure of NetCDF dimension */
typedef struct NdgNcDim{
    char name[NC_MAX_NAME_LEN]; ///< name of dimension
    int length;        ///< length of dimension
    int id;      ///< dimension id in NetCDF file
} NdgNcDim;

/** structure of NetCDF variable */
typedef struct NdgNcVar{
    char name[NC_MAX_NAME_LEN]; ///< name of variables
    int Ndim;       ///< number of dimension in the variable
    int dimId[NC_VAR_MAX_DIM];  ///< pointer to the array of nc_dim
    nc_type type;       ///< variable type
    int id;      ///< variable id in NetCDF file
} NdgNcVar;

/** structure of NetCDF file */
typedef struct NdgNcFile
{
    int IOStepInd;
    char name[NC_MAX_NAME_LEN]; ///< name of files;
    int Ndim;       ///< number of dimensions;
    int Nvar;       ///< number of variable;
    NdgNcDim *dim;  ///< pointer to the array of nc_dim;
    NdgNcVar *var;  ///< pointer to the array of nc_var;
    int ncid;       ///< file id of NetCDF file;
} NdgNcFile;

/***/
void setNcDim(NdgNcDim *dim,
                 const char *name,
                 int len);

/***/
void setNcVar(NdgNcVar *var,
              const char *name,
              int Ndim,
              int *dimId,
              nc_type type);

/***/
void setNcFile(NdgNcFile *file,
               const char *casename,
               int Ndim,
               NdgNcDim *dim,
               int Nvar,
               NdgNcVar *var);

/** print NdgNcFile variable to Mablab command window */
void mexPrintNdgNcFile(NdgNcFile file);

/** read the NetCDF file information to the NdgNcFile structure */
NdgNcFile *readNetcdfFile(const char *filename);

/** close NetCDF file */
void closeNetcdfFile(NdgNcFile file);

/** create the NetCDF file and define its dimensions and variables */
void makeNetcdfFile(NdgNcFile *file);

#endif //NDGNETCDF_H