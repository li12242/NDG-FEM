#include "NdgNetCDF.h"
#include "mex.h"
#include <string.h>

#define ndgPrintf mexPrintf

void setNcDim(NdgNcDim *dim,
              const char *name,
              int len)
{
    
    strcpy(dim->name, name);
    if (len > 0){
        dim->length = len;
    }else if (len == 0){
        dim->length = NC_UNLIMITED;
    }
    return;
}

void setNcVar(NdgNcVar *var,
              const char *name,
              int Ndim,
              int *dimId,
              nc_type type)
{
    
    strcpy(var->name, name);
    if(Ndim > NC_VAR_MAX_DIM){
        fprintf(stderr,
                "%s(%d): The NDim = %d excessed the max dimension number NC_VAR_MAX_DIM = %d.\n",
                __FILE__, __LINE__, Ndim, NC_VAR_MAX_DIM);
    }
    var->Ndim = Ndim; /* number of dimensions */
    for (int n=0; n<Ndim; n++){
        var->dimId[n] = dimId[n];
    }
    switch (type){
        case NC_INT: var->type = NC_INT; break;
        case NC_FLOAT: var->type = NC_FLOAT; break;
        case NC_DOUBLE: var->type = NC_DOUBLE; break;
        case NC_SHORT: var->type = NC_SHORT; break;
        default:
            fprintf(stderr,
                    "%s(%d): Unknown type of NetCDF variable type: %d!\n",
                    __FILE__, __LINE__, (int)type);
    }
}

void setNcFile(NdgNcFile *file,
               const char *casename,
               int Ndim,
               NdgNcDim *dim,
               int Nvar,
               NdgNcVar *var)
{
    if (snprintf(file->name, NC_MAX_NAME_LEN, "%s.nc", casename) < 0){
        fprintf(stderr, "%ps(%d): The casename of the NetCDF file is too long!\n", __FILE__, __LINE__);
    }
    
    file->Ndim = Ndim;
    if ( (dim == NULL) || (var == NULL) ) {
        fprintf(stderr, "%s(%d): The pointer of dimensions or variables is invaild!\n",
                __FILE__, __LINE__);
    }
    file->IOStepInd = 0; // start from 0.
    file->dim = dim;
    file->Nvar = Nvar;
    file->var = var;
    return;
}

/** print NdgNcDim variable to Matlab command window */
static void mexPrintNdgNcDim(NdgNcDim dim){
    mexPrintf("NetCDF dimension: %s, { length = %d }\n",
              dim.name, dim.length);
    return;
}

/** print NdgNcVar variable to Matlab command window */
static void mexPrintNcVar(NdgNcVar var){
    mexPrintf("NetCDF Variable: %s,", var.name);
    mexPrintf("\t{ Ndim = %d, type = %d, dimId = [", var.Ndim, var.type);
    for (int i = 0; i<var.Ndim; i++){
        mexPrintf(" %d", var.dimId[i]);
    }
    mexPrintf("] }\n");
}

/** print NdgNcFile variable to Mablab command window */
void mexPrintNdgNcFile(NdgNcFile file){
    mexPrintf("NetCDF file: %s,", file.name);
    mexPrintf("\t{ Nvar: %d, Ndim: %d}\n", file.Nvar, file.Ndim);

    mexPrintf("Including dimensions:\n");
    for (int n=0; n<file.Ndim; n++){
        mexPrintNdgNcDim(file.dim[n]);
    }
    
    mexPrintf("Including variables:\n");
    for (int i=0; i<file.Nvar; i++){
        mexPrintNcVar(file.var[i]);
        for (int n=0; n<file.var[i].Ndim; n++) {
            int dimId = file.var[i].dimId[n];
            mexPrintf("\tdim[%d]: dimId = %d, name = %s\n",
                      n, file.dim[dimId].id, file.dim[dimId].name);
        }
    }
    return;
}

#define nc_error(t)                                            \
    do                                                         \
    {                                                          \
        if (t != NC_NOERR)                                     \
        {                                                      \
            ndgPrintf("%s: error at line %d: %s\n",            \
                      __FUNCTION__, __LINE__, nc_strerror(t)); \
        }                                                      \
    } while (0)


NdgNcFile *readNetcdfFile(const char *filename){
    NdgNcFile *ncfile = (NdgNcFile*)calloc(1, sizeof(NdgNcFile));
    // read dimensions
    // read variables
    
    return ncfile;
}

/** close NetCDF file */
void closeNetcdfFile(NdgNcFile file){
    nc_error(nc_close(file.ncid));
}

/** create the NetCDF file and define its dimensions and variables */
void makeNetcdfFile(NdgNcFile *file){
    
    const int Ndim = file->Ndim;
    const int Nvar = file->Nvar;
    /* create output file */
    nc_error(nc_create(file->name, NC_CLOBBER | NC_64BIT_OFFSET, &(file->ncid)));
    /* define dimensions */
    for (int i=0; i<Ndim; i++){
        nc_error( nc_def_dim(file->ncid,
                             file->dim[i].name,
                             file->dim[i].length,
                             &(file->dim[i].id)
                             )
                 );
    }
    /* define variables */
    for (int i=0; i < Nvar; i++){
        NdgNcVar *var = file->var+i; /* pointer of var */
        int Ndim = var->Ndim;
        int dimids[Ndim];
        for (int j = 0; j < Ndim; j++){
            int ind = var->dimId[j];
            dimids[j] = file->dim[ind].id;
        }
        nc_error(nc_def_var(file->ncid,
                            var->name,
                            var->type,
                            Ndim,
                            dimids,
                            &(var->id)
                            )
                 );
    }
    /* finish definition */
    nc_error(nc_enddef(file->ncid));
}
