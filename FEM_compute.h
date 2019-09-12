// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ Femcompute.h
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <unistd.h>
//------------constant---------------
extern const double pi;
extern const double mu0;

extern double omega;
//-------------function---------------
void   Elem_compute();
void   Elem_analyse();
void   Elem_assemble();
void   Process_divide();;
void   Matrix_Deform();
void   Generate_SLU_NR();

void   quick_sort_2D_double(int *e1, int *e2, double *n1, double *n2, int low, int high);
int    partition_double(int *e1, int *e2, double *n1, double *n2,int low, int high);
void   swap_double(double *a, double *b);

extern void swap(int *a, int *b);
extern void Create_mpi_newcomm();
extern void Matrix_merge_FRE();

//----------Variables----------------


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

//-------matrix_item------------------
extern int row_start;
extern int row_end;
extern int item_start;
extern int item_end;
extern int item_num;

//-------matrix_before_assemble-------
extern int    *matrix_item_loc;
extern int    *matrix_i_loc;
extern int    *matrix_j_loc;
extern double *matrix_complex_real;
extern double *matrix_complex_imag;

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
//-----------SLU_NR_format------------
int    nnz_loc;
int    m_loc;
int    n_loc;
int    fst_row;
int    *matrix_asub;
int    *matrix_xa ;

