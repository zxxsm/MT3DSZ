// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ post_prepocess.h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <ctype.h>
//------------constant---------------
extern const double pi;
extern const double mu0;

extern double omega;

struct doublecomplex_z
{    double r, i;    };
//-------------function---------------
void post_prepocess(int ifre);

void H_compute(double *real, double *imag);

struct doublecomplex_z Zxy_compute(struct doublecomplex_z ExA, struct doublecomplex_z EyA,\
                           struct doublecomplex_z HxA, struct doublecomplex_z HyA,\
                           struct doublecomplex_z ExB, struct doublecomplex_z EyB,\
                           struct doublecomplex_z HxB, struct doublecomplex_z HyB);

struct doublecomplex_z Zyx_compute(struct doublecomplex_z ExA, struct doublecomplex_z EyA,\
                           struct doublecomplex_z HxA, struct doublecomplex_z HyA,\
                           struct doublecomplex_z ExB, struct doublecomplex_z EyB,\
                           struct doublecomplex_z HxB, struct doublecomplex_z HyB);

int itoa(int n, char s[]);
//----------extern_variables----------
//----------MPI_variables-------------
extern int numprocs;
extern int myrank;
extern int mynewrank;
extern int MPI_FRE;
extern int MPI_FEM;
extern int LU_npcol;
extern int LU_nprow;
//----------time_variables------------
extern double time_s,time_e;
extern double time1,time2;
extern double Prepocess_time;
extern double Fem_compute_time;
extern double Solver_time;
extern double Total_time;

//----------read_variables------------
extern int     bound_mode;

extern int     phy_num;
extern double *phy_cond;

extern int     fre_num;
extern double *frequency;

extern int     layer;
extern double *layer_depth;
extern double *layer_cond;

extern double *boundary_E0_r_loc;
extern double *boundary_E0_i_loc;
extern double *boundary_E0_r;
extern double *boundary_E0_i;

//----------measure_point-------------
extern int     collect_num;
extern double *collect_x;
extern double *collect_y;
extern double *collect_z;

extern int     col_elem_num;
extern int    *col_elem;
extern int    *col_num;

//----------mesh_cordinate------------
extern int     node_number_total;
extern double *node_x;
extern double *node_y;
extern double *node_z;

//-------mesh_boundary_xyz------------
extern double  max_x;
extern double  max_y;
extern double  max_z;
extern double  min_x;
extern double  min_y;
extern double  min_z;

//----------elem_node-----------------
extern int  elem_number_total;

extern int *elem_node1;
extern int *elem_node2;
extern int *elem_node3;
extern int *elem_node4;
extern int *elem_pythsical;

//----------elem_edge-----------------
extern int *elem_edge1;
extern int *elem_edge2;
extern int *elem_edge3;
extern int *elem_edge4;
extern int *elem_edge5;
extern int *elem_edge6;

extern int  edge_number_total;
extern int *edge_node1;
extern int *edge_node2;


//-------boundary_edge_1--------------
extern int  bound_edge_total;
extern int *bound_edge;

//-------boundary_face_2--------------
extern struct boundary_face face_up;
extern struct boundary_face face_down;
extern struct boundary_face face_left;
extern struct boundary_face face_right;
extern struct boundary_face face_forward;
extern struct boundary_face face_behind;

//-------matrix_assemble--------------
extern int     matrix_number;
extern int    *matrix_i;
extern int    *matrix_j;
extern double *matrix_a_r;
extern double *matrix_a_i;

//-----matrix_assemble_deform---------
extern int     num_out_row;
extern int     num_out_column;
extern int     matrix_out_row;
extern int     matrix_out_rc;

extern int    *matrix_i_rh;
extern int    *matrix_j_rh;
extern int    *matrix_index_rh;
extern double *matrix_a_r_rh;
extern double *matrix_a_i_rh;

//---------Matrix_SupLU_solver--------
//-----------SLU_NR_forma-------------
extern int    nnz_loc;
extern int    m_loc;
extern int    n_loc;
extern int    fst_row;
extern int    *matrix_asub;
extern int    *matrix_xa ;

//---------Gather result-------------- 
extern struct doublecomplex_z *b_x_total; 
extern struct doublecomplex_z *b_y_total; 

