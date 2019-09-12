// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ global_varible.c
#include <mpi.h>
#include <math.h>
#include <complex.h>
#include "superlu_zdefs.h" 
//----------constant------------------
const double pi=3.141592653589793;
const double mu0=3.141592653589793*4.0*0.0000001;

double omega;

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
{   double r, i;  };

//----------MPI_variables-------------
int numprocs;
int myrank;
int mynewrank;
int MPI_FRE;
int MPI_FEM;
int LU_npcol;
int LU_nprow;

int new_rank_tran;
int new_rank_comp;
MPI_Comm NewWorld_tran;
MPI_Comm NewWorld_comp;
//----------time_variables------------
double time_s,time_e;
double time1,time2;
double Prepocess_time;
double Fem_compute_time;
double Solver_time;
double Total_time;

//----------read_variables------------
int     bound_mode;

int     phy_num; 
double *phy_cond;

int     fre_num;
double *frequency;

int     layer;
double *layer_depth;
double *layer_cond;


//----------measure_point-------------
int     collect_num;
double *collect_x;
double *collect_y;
double *collect_z;

int     col_elem_num;
int    *col_elem;
int    *col_num;

//----------mesh_cordinate------------
int     node_number_total;
double *node_x;
double *node_y;
double *node_z;

//-------mesh_boundary_xyz------------
double  max_x;
double  max_y;
double  max_z;
double  min_x;
double  min_y;
double  min_z;

//----------elem_node-----------------
int  elem_number_total;

int *elem_node1;
int *elem_node2;
int *elem_node3;
int *elem_node4;
int *elem_pythsical;

//----------elem_edge-----------------
int *elem_edge1;
int *elem_edge2;
int *elem_edge3;
int *elem_edge4;
int *elem_edge5;
int *elem_edge6;

int  edge_number_total;
int *edge_node1;
int *edge_node2;

//-------boundary_edge_1--------------
int  bound_edge_total;
int *bound_edge;

//-------boundary_face_2--------------
struct boundary_face face_up;
struct boundary_face face_down;
struct boundary_face face_left;
struct boundary_face face_right;
struct boundary_face face_forward;
struct boundary_face face_behind;

//-------matrix_item------------------
int row_start;
int row_end;
int item_start;
int item_end;
int item_num;

//-------matrix_before_assemble-------
int    *matrix_item_loc;
int    *matrix_i_loc;
int    *matrix_j_loc;
double *matrix_complex_real;
double *matrix_complex_imag;

//-------matrix_assemble--------------
int     matrix_number;
int    *matrix_i;
int    *matrix_j;
double *matrix_a_r;
double *matrix_a_i;

//-----matrix_assemble_deform---------
int     num_out_row;
int     num_out_column;
int     matrix_out_row;
int     matrix_out_rc;

int    *matrix_i_rh;
int    *matrix_j_rh;
int    *matrix_index_rh;
double *matrix_a_r_rh;
double *matrix_a_i_rh;

//---------Matrix_SupLU_solver--------
//-----------SLU_NR_format------------
int    nnz_loc;
int    m_loc;
int    n_loc;
int    fst_row;
int    *matrix_asub;
int    *matrix_xa ;

doublecomplex *b_x;
doublecomplex *b_y;

//---------Gather result--------------
struct doublecomplex_z *b_x_total;
struct doublecomplex_z *b_y_total;

