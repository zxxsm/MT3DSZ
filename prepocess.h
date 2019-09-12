// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ prepocess.h
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//------------constant---------------
extern const double pi;
extern const double mu0;

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
//-------------function---------------
void    Prepocess();
void    Read_parameters_file();
void    Read_mesh_file();
void    Compute_edges();
void    Compute_boundary_edge();
void    Compute_boundary_face();
void    Task_partition();
void    Find_collect();

struct  boundary_face malloc_face_struct(struct boundary_face face);
struct  boundary_face generate_face(struct boundary_face face, double min_max,double *node_c);
void    verify_face(struct boundary_face face, double min_max, double *node_c);

void    quick_sort_2D(int *e1, int *e2, int *n1 ,int low, int high);
int     partition_int(int *e1, int *e2, int *n1 ,int low, int high);
void    swap(int *a, int *b);
int     min(int a,int b);
int     max(int a,int b);
double  min_arry(double *arry, int length);
double  max_arry(double *arry, int length);

//---------extern_function------------
extern void    Data_distribute();
extern void    matrix_item_distribute(int* sendbuf1, int* sendbuf2, int* sendbuf3,\
               int* sendbuf4,int* sendbuf5, int* sendbuf6,int* sendbuf7, int* sendcounts);

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

