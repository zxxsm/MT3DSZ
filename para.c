// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ para.c 
#include "para.h"

//-------------------------------------------------------------------------------------------
// MPI_Init
//-------------------------------------------------------------------------------------------
void MTSZ_task_init(int *argc, char **argv[])
{

     MPI_Init(argc, argv);
     MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
     MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
     
     printf("MPI_INIT_SUCESS,NUM_PRO=%d,RANK=%d\n",numprocs,myrank);
     MPI_Barrier(MPI_COMM_WORLD); 

     if (myrank==0) printf("MPI_INIT_DOWN!\n");
    
}
//-------------------------------------------------------------------------------------------
// Data_distribute
//-------------------------------------------------------------------------------------------
void Data_distribute()
{
     int par_list[20];
     int i;
     double *temp_d;
     int    *temp_i;
     int    *temp_i2;
     int    l_temp_d,l_temp_i;
     int    temp,temp_s;


     if (myrank==0)
     {
          par_list[0]=MPI_FRE;
          par_list[1]=MPI_FEM;
          par_list[2]=LU_npcol;
          par_list[3]=LU_nprow;
          par_list[4]=bound_mode;
          par_list[5]=phy_num;
          par_list[6]=fre_num;
          par_list[7]=collect_num;
          par_list[8]=layer;          
          par_list[9]=node_number_total;
          par_list[10]=elem_number_total;
          par_list[11]=edge_number_total;
          par_list[12]=max_x;
          par_list[13]=max_y;
          par_list[14]=max_z;
          par_list[15]=min_x;
          par_list[16]=min_y;
          par_list[17]=min_z;
          par_list[18]=col_elem_num;
     }

     MPI_Bcast(par_list,19,MPI_INT,0,MPI_COMM_WORLD);

     if (myrank!=0)
     {
          MPI_FRE=par_list[0]; 
          MPI_FEM=par_list[1];
          LU_npcol=par_list[2];
          LU_nprow=par_list[3];
          bound_mode=par_list[4];
          phy_num=par_list[5];
          fre_num=par_list[6];
          collect_num=par_list[7];
          layer=par_list[8];          
          node_number_total  =par_list[9];
          elem_number_total  =par_list[10];
          edge_number_total  =par_list[11];

          max_x=par_list[12];
          max_y=par_list[13];
          max_z=par_list[14];
          min_x=par_list[15];
          min_y=par_list[16];
          min_z=par_list[17];
          col_elem_num=par_list[18];;

          phy_cond   =(double*) malloc( (phy_num+1)*sizeof(double) );
          frequency  =(double*) malloc( (fre_num+1)*sizeof(double) );
          layer_depth=(double*) malloc( (layer+1)*sizeof(double) );
          layer_cond =(double*) malloc( (layer+1)*sizeof(double) );

          collect_x=(double*) malloc( (collect_num+1)*sizeof(double) );
          collect_y=(double*) malloc( (collect_num+1)*sizeof(double) );
          collect_z=(double*) malloc( (collect_num+1)*sizeof(double) );

          node_x=(double*) malloc( (node_number_total+1)*sizeof(double) );
          node_y=(double*) malloc( (node_number_total+1)*sizeof(double) );
          node_z=(double*) malloc( (node_number_total+1)*sizeof(double) );

          elem_node1=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_node2=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_node3=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_node4=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_pythsical =(int*) malloc( (elem_number_total+1)*sizeof(int) );
          edge_node1=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          edge_node2=(int*) malloc( (edge_number_total+1)*sizeof(int) );

          
     }



     l_temp_d= node_number_total*3+phy_num+fre_num+collect_num*3+layer*2+1;
     l_temp_i= elem_number_total*5+1;

     temp_d=(double*) malloc( l_temp_d*sizeof(double) );
     temp_i=(int   *) malloc( l_temp_i*sizeof(int));

     if (myrank==0)
     {
          //double 
          temp=0;
          for (i=1;i<=node_number_total;i++)
               temp_d[temp+i]=node_x[i];
          temp=temp+node_number_total;
          for (i=1;i<=node_number_total;i++)
               temp_d[temp+i]=node_y[i];
          temp=temp+node_number_total;
          for (i=1;i<=node_number_total;i++)
               temp_d[temp+i]=node_z[i];
          temp=temp+node_number_total;

          for (i=1;i<=phy_num;i++)
               temp_d[temp+i]=phy_cond[i];
          temp=temp+phy_num;
          for (i=1;i<=fre_num;i++)
               temp_d[temp+i]=frequency[i];
          temp=temp+fre_num;
          for (i=1;i<=collect_num;i++)
               temp_d[temp+i]=collect_x[i];
          temp=temp+collect_num;
          for (i=1;i<=collect_num;i++)
               temp_d[temp+i]=collect_y[i];
          temp=temp+collect_num;
          for (i=1;i<=collect_num;i++)
               temp_d[temp+i]=collect_z[i];
          temp=temp+collect_num;
          for (i=1;i<=layer;i++)
               temp_d[temp+i]=layer_depth[i];
          temp=temp+layer;
          for (i=1;i<=layer;i++)
               temp_d[temp+i]=layer_cond[i];

          //int 
          temp=0;
          for (i=1;i<=elem_number_total;i++) temp_i[i+temp]=elem_node1[i];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) temp_i[i+temp]=elem_node2[i];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) temp_i[i+temp]=elem_node3[i];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) temp_i[i+temp]=elem_node4[i];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) temp_i[i+temp]=elem_pythsical[i];

     }
     
    
     MPI_Bcast(temp_d,l_temp_d,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(temp_i,l_temp_i,MPI_INT,0,MPI_COMM_WORLD);

     if (myrank!=0)
     {
          //double
          temp=0;
          for (i=1;i<=node_number_total;i++)
               node_x[i]=temp_d[temp+i];
          temp=temp+node_number_total;
          for (i=1;i<=node_number_total;i++)
               node_y[i]=temp_d[temp+i];
          temp=temp+node_number_total;
          for (i=1;i<=node_number_total;i++)
               node_z[i]=temp_d[temp+i];
          temp=temp+node_number_total;
          for (i=1;i<=phy_num;i++)
               phy_cond[i]=temp_d[temp+i];
          temp=temp+phy_num;
          for (i=1;i<=fre_num;i++)
               frequency[i]=temp_d[temp+i];
          temp=temp+fre_num;
          for (i=1;i<=collect_num;i++)
               collect_x[i]=temp_d[temp+i];
          temp=temp+collect_num;
          for (i=1;i<=collect_num;i++)
               collect_y[i]=temp_d[temp+i];
          temp=temp+collect_num;
          for (i=1;i<=collect_num;i++)
               collect_z[i]=temp_d[temp+i];
          temp=temp+collect_num;
          for (i=1;i<=layer;i++)
               layer_depth[i]=temp_d[temp+i];
          temp=temp+layer;
          for (i=1;i<=layer;i++)
               layer_cond[i]=temp_d[temp+i];

          //int
          temp=0;
          for (i=1;i<=elem_number_total;i++) elem_node1[i]=temp_i[i+temp];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) elem_node2[i]=temp_i[i+temp];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) elem_node3[i]=temp_i[i+temp];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) elem_node4[i]=temp_i[i+temp];
          temp=temp+elem_number_total;
          for (i=1;i<=elem_number_total;i++) elem_pythsical[i]=temp_i[i+temp];

     }


     MPI_Bcast(edge_node1,edge_number_total+1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(edge_node2,edge_number_total+1,MPI_INT,0,MPI_COMM_WORLD);


     /**************verification*****************/

     printf("myrank=%d,layer_depth=%lf,layer_cond=%lf,elem_node4=%d,elem_py=%d\n",\
     myrank,layer_depth[layer],layer_cond[layer],elem_node4[elem_number_total],\
     elem_pythsical[elem_number_total]);

     MPI_Barrier(MPI_COMM_WORLD);
     
     if (bound_mode==1)
     {
          MPI_Bcast(&bound_edge_total,1,MPI_INT,0,MPI_COMM_WORLD);

          if (myrank!=0)
          {
               bound_edge=(int*) malloc((bound_edge_total+1)*sizeof(int) );
          }                     
          MPI_Bcast(bound_edge,bound_edge_total+1 ,MPI_INT,0,MPI_COMM_WORLD);
 
          if (myrank==0) printf("--------------------------------------------------\n");

          printf("myrank=%d,bedge_num=%d,last boundedge=%d\n",myrank,\
          bound_edge_total,bound_edge[bound_edge_total]);
          MPI_Barrier(MPI_COMM_WORLD);
     }
     else if (bound_mode==2)    
     {

          MPI_Bcast(&face_up.bound_face_number,1,MPI_INT,0,MPI_COMM_WORLD);
          MPI_Bcast(&face_down.bound_face_number,1,MPI_INT,0,MPI_COMM_WORLD);
          MPI_Bcast(&face_left.bound_face_number,1,MPI_INT,0,MPI_COMM_WORLD);
          MPI_Bcast(&face_right.bound_face_number,1,MPI_INT,0,MPI_COMM_WORLD);
          MPI_Bcast(&face_forward.bound_face_number,1,MPI_INT,0,MPI_COMM_WORLD);
          MPI_Bcast(&face_behind.bound_face_number,1,MPI_INT,0,MPI_COMM_WORLD);

          if (myrank!=0)
          {
               face_up.bound_face_mode=1;
               face_down.bound_face_mode=2;
               face_left.bound_face_mode=3;
               face_right.bound_face_mode=4;
               face_forward.bound_face_mode=5;
               face_behind.bound_face_mode=6;

               face_up=malloc_face_struct(face_up);
               face_down=malloc_face_struct(face_down);
               face_left=malloc_face_struct(face_left);
               face_right=malloc_face_struct(face_right);
               face_forward=malloc_face_struct(face_forward);
               face_behind=malloc_face_struct(face_behind);
          }

          bound_face_bcast(face_up);
          bound_face_bcast(face_down);
          bound_face_bcast(face_left);
          bound_face_bcast(face_right);
          bound_face_bcast(face_forward);
          bound_face_bcast(face_behind);

          if (myrank==0) printf("--------------------------------------------------\n");
          i=face_behind.bound_face_number;
          printf("myrank=%d,face_num=%d,last e1=%d,e2=%d,e3=%d,el=%d,m=%d,p=%d\n",\
          myrank,i,face_behind.bound_face_edge1[i],face_behind.bound_face_edge2[i],\
          face_behind.bound_face_edge3[i],face_behind.bound_face_elem[i],\
          face_behind.bound_face_mode,face_behind.bound_face_point[i]);
          
     }


     
     MPI_Barrier(MPI_COMM_WORLD);
     if (myrank==0) printf("-------Parameter distribute done.-------\n");
}
//-------------------------------------------------------------------------------------------
// bound_face_bcast
//-------------------------------------------------------------------------------------------
void bound_face_bcast(struct boundary_face face)
{
     int i,j;
     int *temp;

     temp =(int*) malloc((face.bound_face_number+1)*11*sizeof(int) );
     j=face.bound_face_number+1;

     if (myrank==0)
     {
          for (i=0;i<=face.bound_face_number;i++)
               temp[i]=face.bound_face_edge1[i];
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+j]=face.bound_face_edge2[i];
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+2*j]=face.bound_face_edge3[i];
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+3*j]=face.bound_face_elem [i];
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+4*j]=face.bound_face_point[i];          
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+5*j]=face.bound_elem_edge1[i];          
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+6*j]=face.bound_elem_edge2[i];          
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+7*j]=face.bound_elem_edge3[i];          
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+8*j]=face.bound_elem_edge4[i];          
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+9*j]=face.bound_elem_edge5[i];          
          for (i=0;i<=face.bound_face_number;i++)
               temp[i+10*j]=face.bound_elem_edge6[i];          
     }

     MPI_Bcast(temp,11*(face.bound_face_number+1),MPI_INT,0,MPI_COMM_WORLD);

     if (myrank!=0)
     {
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_face_edge1[i]=temp[i];
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_face_edge2[i]=temp[i+j];
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_face_edge3[i]=temp[i+2*j];
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_face_elem [i]=temp[i+3*j];
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_face_point[i]=temp[i+4*j];
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_elem_edge1[i]=temp[i+5*j];          
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_elem_edge2[i]=temp[i+6*j];          
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_elem_edge3[i]=temp[i+7*j];          
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_elem_edge4[i]=temp[i+8*j];          
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_elem_edge5[i]=temp[i+9*j];          
          for (i=0;i<=face.bound_face_number;i++)
               face.bound_elem_edge6[i]=temp[i+10*j];          
     }

     free(temp);     
}
//-------------------------------------------------------------------------------------------
// matrix_item_distribute
//-------------------------------------------------------------------------------------------
void matrix_item_distribute(int* sendbuf1, int* sendbuf2, int* sendbuf3,int* sendbuf4,\
int* sendbuf5, int* sendbuf6,int* sendbuf7, int* sendcounts)
{
     int i;

     /***************matrix_item_distribute******************/
     MPI_Scatter(sendbuf1,1,MPI_INT,&row_start,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Scatter(sendbuf2,1,MPI_INT,&row_end  ,1,MPI_INT,0,MPI_COMM_WORLD);

     MPI_Scatter(sendbuf3,1,MPI_INT,&item_start,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Scatter(sendbuf4,1,MPI_INT,&item_end  ,1,MPI_INT,0,MPI_COMM_WORLD);

     item_num  = item_end-item_start+1;

     matrix_i_loc   = (int*) malloc( item_num*sizeof(int) );
     matrix_j_loc   = (int*) malloc( item_num*sizeof(int) );
     matrix_item_loc= (int*) malloc( item_num*sizeof(int) );

     MPI_Scatterv(sendbuf5,sendcounts,sendbuf3,MPI_INT,matrix_i_loc,item_num,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Scatterv(sendbuf6,sendcounts,sendbuf3,MPI_INT,matrix_j_loc,item_num,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Scatterv(sendbuf7,sendcounts,sendbuf3,MPI_INT,matrix_item_loc,item_num,MPI_INT,0,MPI_COMM_WORLD);

      /*******************verification************************/
     printf("myrank=%d,row_start=%d,row_end=%d,item_start=%d,item_end=%d\n",\
     myrank,row_start,row_end,item_start,item_end);

     MPI_Barrier(MPI_COMM_WORLD);

}
//-------------------------------------------------------------------------------------------
// gather_b
//-------------------------------------------------------------------------------------------
void gather_b()
{    
     int    i,j,k;
     int    *recvbuf;
     int    *displs;
     int    *tasks;
     int     num_task;
     int     bound_num;
     double *temp_b_x;
     double *temp_b_y;

     complex double E0;
     complex double Ex;
     complex double Ey;

     if (mynewrank==0)
     {             
          tasks   = (int*) malloc( MPI_FEM*sizeof(int) );
          displs  = (int*) malloc( MPI_FEM*sizeof(int) );
     }

     MPI_Gather (&m_loc, 1, MPI_INT, tasks, 1 , MPI_INT, 0, NewWorld_comp);
        
     if (mynewrank==0)
     {             
          for (i=0;i<=edge_number_total-1;i++)
          {
               b_x_total[i].r=0;
               b_x_total[i].i=0;
               b_y_total[i].r=0;
               b_y_total[i].i=0;
          }

          num_task=0;
          for (i=0;i<=MPI_FEM-1;i++)
          {
               tasks[i]=tasks[i]*2;
               num_task=num_task+tasks[i];
          }

          
          temp_b_x= (double*) malloc( num_task*sizeof(double) );
          temp_b_y= (double*) malloc( num_task*sizeof(double) );

          displs[0]=0; 
          for (i=1;i<=MPI_FEM-1;i++)
               displs[i]=displs[i-1]+tasks[i-1];
     }


     MPI_Gatherv(b_x, 2*m_loc ,MPI_DOUBLE, temp_b_x, tasks, displs, MPI_DOUBLE,0, NewWorld_comp);
     MPI_Gatherv(b_y, 2*m_loc ,MPI_DOUBLE, temp_b_y, tasks, displs, MPI_DOUBLE,0, NewWorld_comp);

     /*if (myrank==0)
     {
          for (i=0;i<=num_task-1;i++)
               printf("i=%d,temp_b_x=%lf\n",i,temp_b_x[i]);
          for (i=0;i<=num_task-1;i++)
               printf("i=%d,temp_b_y=%lf\n",i,temp_b_y[i]);
     }*/
    
     if (mynewrank==0)
     {
          if (bound_mode==1)
          {                   
               if ( num_task!=(2*edge_number_total-2*bound_edge_total) )
               {   
                    printf( "Mode1 b_gather compute error !\n");
                    exit(0);
               }

               j=1;
               k=1;

               for (i=0;i<=edge_number_total-1;i++)
               {

                    if ((i+1) != bound_edge[j] )
                    {
                         b_x_total[i].r=temp_b_x[(k-1)*2];
                         b_x_total[i].i=temp_b_x[(k-1)*2+1];
                         b_y_total[i].r=temp_b_y[(k-1)*2];
                         b_y_total[i].i=temp_b_y[(k-1)*2+1];
                         k++;
                    }
                    else if ((i+1) ==  bound_edge[j])
                    {
                         bound_num=bound_edge[j];
                         Boundary_E0(bound_num,&E0);
                         field_projection(bound_num,&E0,&Ex,&Ey);
 
                         b_x_total[i].r=creal(Ex);
                         b_x_total[i].i=cimag(Ex);
                         b_y_total[i].r=creal(Ey);
                         b_y_total[i].i=cimag(Ey);

                         if (j<bound_edge_total)
                             j++;
                    }                    
               }
          }
          else if (bound_mode==2)
          {
               if ( num_task!=(2*edge_number_total) )
               {   
                    printf( "Mode2 b_gather compute error !\n");
                    exit(0);
               }
          
               for (i=0;i<=edge_number_total-1;i++) 
               {
               b_x_total[i].r=temp_b_x[i*2];
               b_x_total[i].i=temp_b_x[i*2+1];
               b_y_total[i].r=temp_b_y[i*2];
               b_y_total[i].i=temp_b_y[i*2+1];
               }
          }
     
     }
}
//-------------------------------------------------------------------------------------------
// Create_mpi_newcomm
//-------------------------------------------------------------------------------------------
void Create_mpi_newcomm()
{
     int *rank_tran;
     int *rank_comp;
     int  i;
     int  numpro_tran;
     int  numpro_comp;
     int  color1,color2;

     /*MPI_Group world_group;
     MPI_Group new_group_tran;
     MPI_Group new_group_comp;

     rank_tran=(int*) malloc( MPI_FRE*sizeof(int));
     rank_comp=(int*) malloc( MPI_FEM*sizeof(int));

     rank_tran[0]=(myrank/MPI_FRE)*MPI_FRE;

     if (MPI_FRE>1)
     {
         for (i=1;i<=MPI_FRE-1;i++)
         {
              rank_tran[i]=rank_tran[i-1]+1; 
         }
     }

     rank_comp[0]=myrank%MPI_FRE;
     if (MPI_FEM>1)
     {
         for (i=1;i<=MPI_FEM-1;i++)
         {
              rank_comp[i]=rank_comp[i-1]+MPI_FRE; 
         }
     }
          
     MPI_Comm_group (MPI_COMM_WORLD,&world_group);

     MPI_Group_incl (world_group,MPI_FRE,rank_tran,&new_group_tran);
     MPI_Group_incl (world_group,MPI_FEM,rank_comp,&new_group_comp);

     MPI_Comm_create(MPI_COMM_WORLD,new_group_tran,&NewWorld_tran);
     MPI_Comm_create(MPI_COMM_WORLD,new_group_comp,&NewWorld_comp);
     */
     color1=myrank/MPI_FRE;
     color2=myrank%MPI_FRE;
     
     MPI_Comm_split(MPI_COMM_WORLD,color1,myrank,&NewWorld_tran);
     MPI_Comm_split(MPI_COMM_WORLD,color2,myrank,&NewWorld_comp);

     MPI_Comm_rank(NewWorld_tran,&new_rank_tran);
     MPI_Comm_rank(NewWorld_comp,&new_rank_comp);

     MPI_Comm_size(NewWorld_tran,&numpro_tran);
     MPI_Comm_size(NewWorld_comp,&numpro_comp);

     printf("Tran:myrank=%d,newrank=%d,ranks=%d\n",myrank,new_rank_tran,numpro_tran);
     MPI_Barrier(MPI_COMM_WORLD);
     printf("Comp:myrank=%d,newrank=%d,ranks=%d\n",myrank,new_rank_comp,numpro_comp);
     MPI_Barrier(MPI_COMM_WORLD);
}
//-------------------------------------------------------------------------------------------
//  Matrix_merge_FRE 
//-------------------------------------------------------------------------------------------
void Matrix_merge_FRE()
{
     int     i;
     int    *matrix_num_arry;    
     int     matrix_number_new;
     int    *temp_i;
     int    *temp_j;
     double *temp_xr;
     double *temp_xi;

     int    *displs;
     int    *tasks;
     int     num_task;
     
     matrix_num_arry =(int*) malloc( MPI_FRE*sizeof(int));
     displs          =(int*) malloc( MPI_FRE*sizeof(int));

     MPI_Allgather(&matrix_number,1,MPI_INT,matrix_num_arry,1,MPI_INT,NewWorld_tran);

     matrix_number_new=0;

     for (i=0;i<=(MPI_FRE-1);i++)
          matrix_number_new=matrix_number_new+matrix_num_arry[i];

     printf("myrank=%d,total=%d\n",myrank,matrix_number_new);
     displs[0]=0;

     if (MPI_FRE > 1)
     {
          for (i=0;i<=(MPI_FRE-2);i++)
              displs[i+1]=displs[i]+matrix_num_arry[i];
     }     

     //printf("myrank=%d,mtx_num=%d,mtx_num_new=%d\n",myrank,matrix_number,matrix_number_new);
    
     temp_i =(int*) malloc(matrix_number_new *sizeof(int));
     temp_j =(int*) malloc(matrix_number_new *sizeof(int));
     temp_xr=(double*) malloc(matrix_number_new *sizeof(double));
     temp_xi=(double*) malloc(matrix_number_new *sizeof(double));

     MPI_Allgatherv(matrix_i,matrix_number,MPI_INT,temp_i,matrix_num_arry,displs,MPI_INT,NewWorld_tran);
     MPI_Allgatherv(matrix_j,matrix_number,MPI_INT,temp_j,matrix_num_arry,displs,MPI_INT,NewWorld_tran);

     MPI_Allgatherv(matrix_a_r,matrix_number,MPI_DOUBLE,temp_xr,matrix_num_arry,displs,MPI_DOUBLE,NewWorld_tran);
     MPI_Allgatherv(matrix_a_i,matrix_number,MPI_DOUBLE,temp_xi,matrix_num_arry,displs,MPI_DOUBLE,NewWorld_tran);

     /*if (myrank==15)
     {
          printf("displs1=%d,2=%d,3=%d,4=%d \n",displs[0],displs[1],displs[2],displs[3] );    
          for (i=0;i<matrix_number;i++)
               printf("beforei=%d,matrix_a_r=%lf\n",i,matrix_a_r[i]);
     }*/
     

     free(matrix_i);
     free(matrix_j);
     free(matrix_a_r);
     free(matrix_a_i);

     matrix_number=matrix_number_new;     

         
     matrix_i =(int*) malloc(matrix_number *sizeof(int));
     matrix_j =(int*) malloc(matrix_number *sizeof(int));
     matrix_a_r=(double*) malloc(matrix_number *sizeof(double));
     matrix_a_i=(double*) malloc(matrix_number *sizeof(double));

     memcpy (matrix_i,temp_i,matrix_number *sizeof(int));
     memcpy (matrix_j,temp_j,matrix_number *sizeof(int));
     memcpy (matrix_a_r,temp_xr,matrix_number *sizeof(double));
     memcpy (matrix_a_i,temp_xi,matrix_number *sizeof(double));

     /*if (myrank==15)
     {   
          for (i=0;i<matrix_number;i++)
               printf("after1i=%d,matrix_a_r=%lf\n",i,matrix_a_r[i]);
     }*/

     free(temp_i);
     free(temp_j);
     free(temp_xr);
     free(temp_xi);
     
     if (new_rank_comp==0&& myrank!=0)
     {
          elem_edge1=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          elem_edge2=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          elem_edge3=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          elem_edge4=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          elem_edge5=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          elem_edge6=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          col_elem=(int*)malloc((col_elem_num)*sizeof(int)); 
          col_num =(int*)malloc((col_elem_num)*sizeof(int));
     }

     
     if (new_rank_comp==0)
     {
          MPI_Bcast(elem_edge1,elem_number_total+1,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(elem_edge2,elem_number_total+1,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(elem_edge3,elem_number_total+1,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(elem_edge4,elem_number_total+1,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(elem_edge5,elem_number_total+1,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(elem_edge6,elem_number_total+1,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(col_elem,col_elem_num,MPI_INT,0,NewWorld_tran);
          MPI_Bcast(col_num,col_elem_num,MPI_INT,0,NewWorld_tran);
     }

     MPI_Comm_free(&NewWorld_tran);
	 
     MPI_Comm_rank(NewWorld_comp,&mynewrank);
     //MPI_Barrier(MPI_COMM_WORLD);
     //exit(0);

}
