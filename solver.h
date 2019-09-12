// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ sovler.h
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "superlu_zdefs.h"

//------------constant---------------
extern const double pi;
extern const double mu0;

extern double omega;

struct boundary_face
{
    int bound_face_number;
    int bound_face_mode;
    int *bound_face_edge1;
    int *bound_face_edge2;
    int *bound_face_edge3;
    int *bound_face_elem;
    int *bound_face_point;
    int *bound_elem_edge1;
    int *bound_elem_edge2;
    int *bound_elem_edge3;
    int *bound_elem_edge4;
    int *bound_elem_edge5;
    int *bound_elem_edge6;

};

struct doublecomplex_z
{    double r, i;    };


//-------------function---------------
void   Main_Solver();

void   b_generate_Dirichlet();
void   b_generate_Newman();

void   generate_condition(struct boundary_face face);

void   p1_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b2,double b3,double b4,double c2,double c3 ,double c4,double d2,double d3,double d4,\
int i,struct boundary_face);
void   p2_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b1,double b3,double b4,double c1,double c3 ,double c4,double d1,double d3,double d4,\
int i,struct boundary_face);
void   p3_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b1,double b2,double b4,double c1,double c2 ,double c4,double d1,double d2,double d4,\
int i,struct boundary_face);
void   p4_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b1,double b2,double b3,double c1,double c2 ,double c3,double d1,double d2,double d3,\
int i,struct boundary_face);

void   Boundary_E0(int edge_ind,complex double *E0);
void   Boundary_H0(int edge_ind,complex double *H0);
void   field_projection(int edge_ind,complex double *E0,complex double *Ex,complex double *Ey);

extern void   gather_b();
extern void   post_prepocess(int ifre);


//----------extern_variables----------
//----------MPI_variables-------------
extern int numprocs;
extern int myrank;
extern int mynewrank;
extern int MPI_FRE;
extern int MPI_FEM;
extern int LU_npcol;
extern int LU_nprow;

extern int new_rank_tran;
extern int new_rank_comp;
extern MPI_Comm NewWorld_tran;
extern MPI_Comm NewWorld_comp;

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

extern doublecomplex *b_x;
extern doublecomplex *b_y;

//---------Gather result-------------- 
extern struct doublecomplex_z *b_x_total;
extern struct doublecomplex_z *b_y_total;

