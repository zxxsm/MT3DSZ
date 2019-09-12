// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ prepocess.c
#include "prepocess.h"

//-------------------------------------------------------------------------------------------
// Prepocess
// 1.read_param
// 2.read_mesh
// 3.compute_edge
// 4.compute_boundary_edge 
//   compute_boundary_face 
// 5.data_distribute
// 6.Task_partition
// 7.Find_collect
//-------------------------------------------------------------------------------------------
void Prepocess()
{
     if (myrank==0)
     {
          printf("****************************************\n");
          printf("           Prepocess_start              \n"); 
          printf("****************************************\n"); 
     }

     if (myrank==0)
     {
          time1=MPI_Wtime();
          time_s=time1;
     }
     // ----------------Read parameters file------------//
     Read_parameters_file();

     // ----------------Read mesh file------------------//
     Read_mesh_file();

     // ----------------Compute edges ------------------//
     Compute_edges();

     // ----------------Compute boundary edge-----------//
     if (bound_mode==1)
     {
          Compute_boundary_edge();
     }
     // ----------------Compute boundary face-----------//
     else if (bound_mode==2)
     {
          Compute_boundary_face();
     }
     // ----------------Find_collect-------------------//
     Find_collect();

     // ----------------Data distribute-----------------//
     Data_distribute();

     // ----------------Task_partition------------------//
     Task_partition();

     MPI_Barrier(MPI_COMM_WORLD);
     if (myrank==0) 
     {
          time2=MPI_Wtime();
          Prepocess_time=time2-time1;
     }

     if (myrank==0)
     {
          printf("****************************************\n");
          printf("           Prepocess_done               \n"); 
          printf("****************************************\n"); 
     }
}
//-------------------------------------------------------------------------------------------
// Read_parameters_from_file
//-------------------------------------------------------------------------------------------
void Read_parameters_file()
{
     char  para_infile[50];
     char  message[20];
     int   temp;
     int   i;

     FILE *prb_file;

     if (myrank==0)
     {	

          printf("----------Input Parameter Set.----------\n"); 
          strcpy(para_infile,"MT-inp");
          if((prb_file=fopen(para_infile, "r"))==NULL)
          {
               printf("\nCannot open: %s\n", para_infile);
               printf("Please ensure the file is within the directory\n");
               exit(-1);
          }           
          else
          {
               while(!feof(prb_file))
               {
                    if (fscanf(prb_file,"%s",message)!=1)
                      printf("message read start!\n");
 

                    /***************$MPI_FREQUENCY*****************/
                    if (strcmp(message,"$MPI_FEM")==0)
                    {
                         if (fscanf(prb_file,"%d",&MPI_FEM)!=1)
                              printf("MPI_FEM read error!\n");
                    }
                    /***************$MPI_FREQUENCY*****************/
                    if (strcmp(message,"$MPI_FRE")==0)
                    {
                         if (fscanf(prb_file,"%d",&MPI_FRE)!=1)
                              printf("MPI_FRE read error!\n");
                    }

                    /***************$MPI_nprow_npcol***************/
                    if (strcmp(message,"$LU_NPROW")==0)
                    {
                         if (fscanf(prb_file,"%d",&LU_nprow)!=1);
                    }

                    if (strcmp(message,"$LU_NPCOL")==0)
                    {
                         fscanf(prb_file,"%d",&LU_npcol);
                    }

                    /***************$BOUNDARY_CONDITIO*************/
                    if (strcmp(message,"$BOUNDARY_CONDITION")==0)
                    {
                         fscanf(prb_file,"%d",&bound_mode);
                    }

                    /***************$PHYSICAL_AREA*****************/
                    if (strcmp(message,"$PHYSICAL_AREA_COND")==0)
                    {
                         fscanf(prb_file,"%d",&phy_num);
                         phy_cond=(double*) malloc( (phy_num+1)*sizeof(double) );
          
                         for (i=1;i<=phy_num;i++)
                         {
                              fscanf(prb_file,"%d%lf",&temp,&phy_cond[i]);
                         }
          
                         fscanf(prb_file,"%s",message);
                         if (strcmp(message,"$END_PHY")!=0)
                         {
                              printf("PHY_COND data read error!\n");
                              exit(-1);
                         }
                    }
                    /*******************$FREQUENCY*****************/
                    if (strcmp(message,"$FREQUENCY")==0)
                    {
                         fscanf(prb_file,"%d",&fre_num);
                         frequency=(double*) malloc( (fre_num+1)*sizeof(double) );
                         for (i=1;i<=fre_num;i++)
                         {
                              fscanf(prb_file,"%d%lf",&temp,&frequency[i]);
                         }
          
                         fscanf(prb_file,"%s",message);
                         if (strcmp(message,"$END_FREQUENCY")!=0)
                         {
                              printf("FREQUENCY data read error!\n");
                              exit(-1);
                         }
                    }
                     
                    /***************$LOCATION********************/
                    if (strcmp(message,"$COLLECT_POINT_XYZ")==0)
                    {    
                         fscanf(prb_file,"%d",&collect_num);
                         collect_x=(double*) malloc( (collect_num+1)*sizeof(double) );
                         collect_y=(double*) malloc( (collect_num+1)*sizeof(double) );
                         collect_z=(double*) malloc( (collect_num+1)*sizeof(double) );
          
                         for (i=1;i<=collect_num;i++)
                         {    
                              fscanf(prb_file,"%d%lf%lf%lf",&temp,&collect_x[i],\
   	                      &collect_y[i],&collect_z[i]);
                         }
                         
                         fscanf(prb_file,"%s",message);
                         if (strcmp(message,"$END_COLLECT")!=0)
                         {    
                              printf("$LOCATION data read error!\n");
                              exit(-1);
                         }
                    }
                    /***************$BACKGROUND********************/
                    if (strcmp(message,"$LAYER_DEPTH_CONDUCTION")==0)
                    {    
                        fscanf(prb_file,"%d",&layer);
                        layer_depth=(double*) malloc( (layer+1)*sizeof(double) );
                        layer_cond =(double*) malloc( (layer+1)*sizeof(double) );
                        for (i=1;i<=layer;i++)
                        {
   	                     fscanf(prb_file,"%d%lf%lf",&temp,&layer_depth[i],\
                             &layer_cond[i]);
                        }
                        fscanf(prb_file,"%s",message);
                        if (strcmp(message,"$END_BACK_GROUND")!=0)
                        {
                             printf("$LOCATION data read error!\n");
                             exit(-1);
                        }
          
                    }
          
               }//end_while
          }
          fclose(prb_file);

          /********************verification*******************/
                     
          if(LU_nprow*LU_npcol!=MPI_FEM)
          {
               printf("MPI_PROCESS_LU set error!\n");
               exit(-1);
          }
          if(numprocs!=MPI_FEM*MPI_FRE)
          {
               printf("MPI_PROCESS_NUM set error!\n");
               exit(-1);
          }
          if(MPI_FRE > fre_num)
          {
               printf("MPI_FRE set error\n");
               exit(-1);
          }
          if(fre_num%MPI_FRE!=0)
          {
               printf("MPI_FRE best set to be divisible by read_par.fre_num");
          }
          
          printf("MPI_FEM=%d\n",MPI_FEM);
          printf("MPI_FRE=%d\n",MPI_FRE);
          printf("MPI_LU: LU_nprow=%d  LU_npcol=%d\n",LU_nprow,LU_npcol);
         
          printf("Bound_mod=%d\n",bound_mode);
          if (bound_mode==1) 
          {   printf("Load Dirichlet Boundary Condition.\n"); }
          else if (bound_mode==2)
          {   printf("Load Newman Boundary Condition.\n"); }
          else
          {   printf("Boundary Condition set error !.\n"); exit(0); }


          printf("phy_num=%d\n",phy_num);
          for (i=1;i<=phy_num;i++) printf("phy_cond[%d]=%.16lf\n",i,phy_cond[i]);
    
          printf("layer_num=%d\n",layer);
          for (i=1;i<=layer;i++) printf("layer_depth[%d]=%lf,layer_cond[%d]=%.16lf\n",\
          i,layer_depth[i],i,layer_cond[i]);
          
          for (i=1;i<=layer;i++)
          {
              if (layer_cond[i]!=phy_cond[i])
              {printf("Background conduction set error !.\n"); exit(0);}
          }

          printf("fre_num=%d\n",fre_num); 
          for (i=1;i<=fre_num;i++) printf("frequency[%d]=%lf\n",i,frequency[i]);
          
          printf("collect_num=%d\n",collect_num);
          for (i=1;i<=collect_num;i++) printf("col_x[%d]=%lf,col_y[%d]=%lf,col_z[%d]=%lf\n",\
          i,collect_x[i],i,collect_y[i],i,collect_z[i]);
          
          
          printf("--------Input file read done.-----------\n"); 
     }     
}
//-------------------------------------------------------------------------------------------
// Read_mesh
//-------------------------------------------------------------------------------------------
void Read_mesh_file()
{
     int     i,j;
     int     physical_number;
     int     type;
     int     temp,temp1;
     int     elementary;
     char    message[20];
     char    mesh_infile[20];
     FILE   *mesh_file;

     
     if (myrank==0)
     {
          printf("---------Mesh file read start.----------\n"); 

          strcpy(mesh_infile,"MT.msh");
          if((mesh_file=fopen(mesh_infile, "r"))==NULL)
          {
               printf("\nCannot open: %s\n", mesh_infile);
               printf("Please ensure the file is within the directory.\n");
               exit(-1);
          }
          else
          {
               while(!feof(mesh_file))
               {
                    fscanf(mesh_file,"%s",message);
          
                    /***********$PhysicalNames*************/
                    if (strcmp(message,"$PhysicalNames")==0)
                    {
                         fscanf(mesh_file,"%d",&physical_number); 
                    }
                    /***************$Nodes*****************/
                    if (strcmp(message,"$Nodes")==0)
                    {
                         fscanf(mesh_file,"%d",&node_number_total);
                         node_x=(double*) malloc( (node_number_total+1)*sizeof(double) );
                         node_y=(double*) malloc( (node_number_total+1)*sizeof(double) );
                         node_z=(double*) malloc( (node_number_total+1)*sizeof(double) );
                         node_x[0]=0;
                         node_y[0]=0;
                         node_z[0]=0;
          
                         for (i=1;i<=node_number_total;i++)
                         {
                              fscanf(mesh_file,"%d%lf%lf%lf",&j,&node_x[i],&node_y[i],&node_z[i]);
                         }
          
                         fscanf(mesh_file,"%s",message);
                         if (strcmp(message,"$EndNodes")!=0)
                         {
                              printf("Nodes data read error!\n");
                              exit(-1);
                         }
                    }
                    /*************$Elemets******************/
                    if (strcmp(message,"$Elements")==0)
                    {
                         fscanf(mesh_file,"%d",&elem_number_total);
                         elem_node1=(int*) malloc( (elem_number_total+1)*sizeof(int) );
                         elem_node2=(int*) malloc( (elem_number_total+1)*sizeof(int) );
                         elem_node3=(int*) malloc( (elem_number_total+1)*sizeof(int) );
                         elem_node4=(int*) malloc( (elem_number_total+1)*sizeof(int) );
                         elem_pythsical=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          
                         elem_node1[0]=0;
                         elem_node2[0]=0;
                         elem_node3[0]=0;
                         elem_node4[0]=0;
                         elem_pythsical[0]=0;
          
                         for (i=1;i<=elem_number_total;i++)
                         {
                              fscanf(mesh_file,"%d%d%d%d%d%d%d%d%d",&j,&type,&elementary,&elem_pythsical[i],
                              &elementary,&elem_node1[i],&elem_node2[i],&elem_node3[i],&elem_node4[i]);
                         } 
          
                         fscanf(mesh_file,"%s",message);
                         if (strcmp(message,"$EndElements")!=0)
                         {
                              printf("Elements data read error!\n");
                              exit(-1);
                         }
          
                    }
          
               }//end_while
          
          fclose(mesh_file);
          /***********************verification**************************/
          if (phy_num!=physical_number)
          {     
               printf("Phy_num are not the same in mesh and inp file.!\n");
               exit(-1);
          }
          
          printf("pysical_number=%d\n",physical_number);

          i=node_number_total;          
          printf("Node_number_total=%d\n",node_number_total);
          printf("The last node xyz:\n");
          printf("node_x[i]=%.16lf\nnode_y[i]=%.16lf\nnode_z[i]=%.16lf\n",\
          node_x[i],node_y[i],node_z[i]);
          
          i=elem_number_total;           
          printf("Elem_number_total=%d\n",elem_number_total);

          printf("The last elem node:\n");          
          printf("node_1[i]=%d\nnode_2[i]=%d\nnode_3[i]=%d\nnode_4[i]=%d\nnode_pyth[i]=%d\n",\
          elem_node1[i],elem_node2[i],elem_node3[i],elem_node4[i],elem_pythsical[i]);
          
          printf("---------Read mesh file done.-----------\n"); 
          
          
          }
     }      
}
//------------------------------------------------------------------------------------------
// Compute_edges_from_nodes
//------------------------------------------------------------------------------------------
void Compute_edges()
{

     int i,j,k;
     int temp_arr[5];
     int temp;

     int *edges1;
     int *edges2;
     int *edge_num;
     int *edge_num_index;

     if (myrank==0)
     {	 
          printf("--------Edges and Boundary compute.-----\n"); 

          /**********element_rearrange_small_to_large******/
          edges1=(int*) malloc( (elem_number_total*6+1)*sizeof(int) );
          edges2=(int*) malloc( (elem_number_total*6+1)*sizeof(int) );
          
          edge_num      =(int*) malloc( (elem_number_total*6+1)*sizeof(int) );
          edge_num_index=(int*) malloc( (elem_number_total*6+1)*sizeof(int) );
          
          for (i=0;i<=elem_number_total*6;i++)
          {
               edge_num[i]=i;
          }
          
          edges1[0]=0;
          edges2[0]=0;
          
          for (i=1;i<=elem_number_total;i++)
          {
               edges1[(i-1)*6+1]=elem_node1[i];
               edges1[(i-1)*6+2]=elem_node1[i];
               edges1[(i-1)*6+3]=elem_node1[i];
               edges1[(i-1)*6+4]=elem_node2[i];
               edges1[(i-1)*6+5]=elem_node2[i];
               edges1[(i-1)*6+6]=elem_node3[i];
          
               edges2[(i-1)*6+1]=elem_node2[i];
               edges2[(i-1)*6+2]=elem_node3[i];
               edges2[(i-1)*6+3]=elem_node4[i];
               edges2[(i-1)*6+4]=elem_node3[i];
               edges2[(i-1)*6+5]=elem_node4[i];
               edges2[(i-1)*6+6]=elem_node4[i];
          }
          
          for (i=1;i<=elem_number_total*6;i++)
          {
               if (edges1[i]>edges2[i])
               {
                   temp=edges1[i];
                   edges1[i]=edges2[i];
                   edges2[i]=temp;
               }
          } 
          
          quick_sort_2D(edges1,edges2,edge_num,1,elem_number_total*6);
          
          temp=1;
          edge_num_index[1]=temp;
          /**************Remove_duplicates*************/
          for (i=1;i<=elem_number_total*6-1;i++)
          { 
               if (edges1[i] != edges1[i+1] || edges2[i] != edges2[i+1])
               { 
                    temp++;
               }
               edge_num_index[i+1]=temp;
          }
                    
          edge_number_total=temp;
          edge_node1=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          edge_node2=(int*) malloc( (edge_number_total+1)*sizeof(int) );
          
          temp=1;
          edge_node1[1]=edges1[1];
          edge_node2[1]=edges2[1];
          
          for (i=1;i<=elem_number_total*6-1;i++)
          { 
               if (edges1[i] != edges1[i+1] || edges2[i] != edges2[i+1])
               { 
                    edge_node1[temp+1]=edges1[i+1];
                    edge_node2[temp+1]=edges2[i+1];
                    temp++;
               }
          }
          
          elem_edge1=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_edge2=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_edge3=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_edge4=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_edge5=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          elem_edge6=(int*) malloc( (elem_number_total+1)*sizeof(int) );
          
          for (i=1;i<=elem_number_total;i++)
          {
               elem_edge1[i]=-1;
               elem_edge2[i]=-1;
               elem_edge3[i]=-1;
               elem_edge4[i]=-1;
               elem_edge5[i]=-1;
               elem_edge6[i]=-1;
          }
          
          for (i=1;i<=elem_number_total*6;i++)
          {
               j   =edge_num[i]/6 +1;
               temp=edge_num[i]%6   ;
               if (temp==0) j--;
          
               switch(temp)
               {
                    case 1:elem_edge1[j]=edge_num_index[i];break;
                    case 2:elem_edge2[j]=edge_num_index[i];break;
                    case 3:elem_edge3[j]=edge_num_index[i];break;
                    case 4:elem_edge4[j]=edge_num_index[i];break;
                    case 5:elem_edge5[j]=edge_num_index[i];break;
                    case 0:elem_edge6[j]=edge_num_index[i];break;
          
               }             
          }
          /**************verification*****************/
          for (i=1;i<=elem_number_total;i++)
          {    
               if (elem_edge1[i]==-1 || elem_edge2[i]==-1 ||\
                   elem_edge3[i]==-1 || elem_edge4[i]==-1 ||\
                   elem_edge5[i]==-1 || elem_edge6[i]==-1)
               { 
                    printf("Elements edge compute error!\n");
                    exit(-1);
               }
          }
                    
          free(edges1);
          free(edges2);
          free(edge_num);
          free(edge_num_index);
          
          
          printf("EDGES_number_total=%d\n",edge_number_total);
          for (i=1;i<=elem_number_total;i++)
          {
               
               if ((edge_node1[elem_edge1[i]]+edge_node2[elem_edge1[i]])!=\
                   (elem_node1[i]+elem_node2[i])||\
                   (edge_node1[elem_edge6[i]]+edge_node2[elem_edge6[i]])!=\
                   (elem_node3[i]+elem_node4[i]))
               { 
                    printf("Edge node compute error!\n");
                    exit(-1); 
               }
          }
     }            
}
//------------------------------------------------------------------------------------------
// Compute_Boundary_edge
//------------------------------------------------------------------------------------------
void Compute_boundary_edge()
{
     int temp;
     int i,j;
     int ind1,ind2;     
     
     if (myrank==0)
     { 
          /**************find min of Z*****************/
          min_x=min_arry(node_x, node_number_total );
          min_y=min_arry(node_y, node_number_total );
          min_z=min_arry(node_z, node_number_total );
          
          max_x=max_arry(node_x, node_number_total );
          max_y=max_arry(node_y, node_number_total );
          max_z=max_arry(node_z, node_number_total );

          printf("model_bound range:\n");
          printf("min_x=%lf,max_x=%lf\n",min_x,max_x); 
          printf("min_y=%lf,max_y=%lf\n",min_y,max_y); 
          printf("min_z=%lf,max_z=%lf\n",min_z,max_z); 
          
          /********************************************/
          temp=0;
          for (i=1;i<=edge_number_total;i++)
          {
               ind1=edge_node1[i];
               ind2=edge_node2[i];
          
               if      (node_x[ind1]==min_x && node_x[ind2]==min_x) temp++ ;
               else if (node_x[ind1]==max_x && node_x[ind2]==max_x) temp++ ;
               else if (node_y[ind1]==min_y && node_y[ind2]==min_y) temp++ ;
               else if (node_y[ind1]==max_y && node_y[ind2]==max_y) temp++ ;
               else if (node_z[ind1]==min_z && node_z[ind2]==min_z) temp++ ;
               else if (node_z[ind1]==max_z && node_z[ind2]==max_z) temp++ ;
          
          }

          bound_edge_total=temp;
          bound_edge=(int*) malloc((bound_edge_total+1)*sizeof(int) );
          
          temp=0; 
          for (i=1;i<=edge_number_total;i++)
          {
               ind1=edge_node1[i];
               ind2=edge_node2[i];
          
               if (node_x[ind1]==min_x && node_x[ind2]==min_x)
               {   
                    temp++;
                    bound_edge[temp]=i;
               }
               else if (node_x[ind1]==max_x && node_x[ind2]==max_x) 
               {   
                    temp++;
                    bound_edge[temp]=i;
               }
               else if (node_y[ind1]==min_y && node_y[ind2]==min_y) 
               {
                    temp++;
                    bound_edge[temp]=i;
               }
               else if (node_y[ind1]==max_y && node_y[ind2]==max_y) 
               {
                    temp++;
                    bound_edge[temp]=i;
               }
               else if (node_z[ind1]==min_z && node_z[ind2]==min_z) 
               {
                    temp++;
                    bound_edge[temp]=i;
               }
               else if (node_z[ind1]==max_z && node_z[ind2]==max_z) 
               {
                    temp++;
                    bound_edge[temp]=i;
               }
               
          }
          
          /**************verification*****************/
          printf("Boundary_edge_num=%d\n",bound_edge_total);


          for (i=1;i<=bound_edge_total;i++)
          {
               ind1=edge_node1[bound_edge[i]];
               ind2=edge_node2[bound_edge[i]];
               if ((node_x[ind1]==max_x && node_x[ind2]==max_x)||\
                   (node_x[ind1]==min_x && node_x[ind2]==min_x)||\
                   (node_y[ind1]==max_y && node_y[ind2]==max_y)||\
                   (node_y[ind1]==min_y && node_y[ind2]==min_y)||\
                   (node_z[ind1]==max_z && node_z[ind2]==max_z)||\
                   (node_z[ind1]==min_z && node_z[ind2]==min_z)
                  )
               {
                   ;                 
               }
               else 
               {
                    printf("Bounday edge compute error!\n");
                    exit(-1);
               }
          }
          printf("----Edges and Boundary compute done.----\n"); 
          
     }     
}
//------------------------------------------------------------------------------------------
// Compute_Boundary_face
//------------------------------------------------------------------------------------------
void Compute_boundary_face()
{
     int temp;
     int i,j;
     int tempu,tempd,templ,tempr,tempf,tempb;
     int tempud,templr,tempfb;
     int fu,fd,fl,fr,ff,fb;
     int node1,node2,node3,node4;
     int n1,n2,n3,n4;

     if (myrank==0)
     {
          
          /**************find min of Z*****************/
          min_x=min_arry(node_x, node_number_total );
          min_y=min_arry(node_y, node_number_total );
          min_z=min_arry(node_z, node_number_total );
          
          max_x=max_arry(node_x, node_number_total );
          max_y=max_arry(node_y, node_number_total );
          max_z=max_arry(node_z, node_number_total );

          printf("model_bound range:\n");
          printf("min_x=%lf,max_x=%lf\n",min_x,max_x); 
          printf("min_y=%lf,max_y=%lf\n",min_y,max_y); 
          printf("min_z=%lf,max_z=%lf\n",min_z,max_z); 

          /********************************************/
          tempu=0; tempd=0;
          templ=0; tempr=0;
          tempf=0; tempb=0;
          
          for (i=1;i<=elem_number_total;i++)
          {
               fu=0;fd=0;
               fl=0;fr=0;
               ff=0;fb=0;
          
               node1=elem_node1[i];
               node2=elem_node2[i];
               node3=elem_node3[i];
               node4=elem_node4[i];
          
               if (node_z[node1]==max_z) fu++ ;
               if (node_z[node2]==max_z) fu++ ;
               if (node_z[node3]==max_z) fu++ ;
               if (node_z[node4]==max_z) fu++ ;
               if (fu==3) tempu++;
          
               if (node_z[node1]==min_z) fd++ ;
               if (node_z[node2]==min_z) fd++ ;
               if (node_z[node3]==min_z) fd++ ;
               if (node_z[node4]==min_z) fd++ ;
               if (fd==3) tempd++;
          
               if (node_x[node1]==max_x) fr++ ;
               if (node_x[node2]==max_x) fr++ ;
               if (node_x[node3]==max_x) fr++ ;
               if (node_x[node4]==max_x) fr++ ;
               if (fr==3) tempr++;
          
               if (node_x[node1]==min_x) fl++ ;
               if (node_x[node2]==min_x) fl++ ;
               if (node_x[node3]==min_x) fl++ ;
               if (node_x[node4]==min_x) fl++ ;
               if (fl==3) templ++;
          
               if (node_y[node1]==max_y) ff++ ;
               if (node_y[node2]==max_y) ff++ ;
               if (node_y[node3]==max_y) ff++ ;
               if (node_y[node4]==max_y) ff++ ;
               if (ff==3) tempf++;
          
               if (node_y[node1]==min_y) fb++ ;
               if (node_y[node2]==min_y) fb++ ;
               if (node_y[node3]==min_y) fb++ ;
               if (node_y[node4]==min_y) fb++ ;
               if (fb==3) tempb++;
          
          }

          face_up.bound_face_number=tempu;
          face_down.bound_face_number=tempd;
          face_left.bound_face_number=templ;
          face_right.bound_face_number=tempr;
          face_forward.bound_face_number=tempf;
          face_behind.bound_face_number=tempb;

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
          
          face_up=generate_face(face_up,max_z,node_z);
          face_down=generate_face(face_down,min_z,node_z);
          face_left=generate_face(face_left,min_x,node_x);
          face_right=generate_face(face_right,max_x,node_x);
          face_forward=generate_face(face_forward,max_y,node_y);
          face_behind=generate_face(face_behind,min_y,node_y);
          
          /**************verification*****************/
          verify_face(face_up,max_z,node_z);
          verify_face(face_down,min_z,node_z);
          verify_face(face_left,min_x,node_x);
          verify_face(face_right,max_x,node_x);
          verify_face(face_forward,max_y,node_y);
          verify_face(face_behind,min_y,node_y);

          printf("----Edges and Boundary compute done.----\n"); 
     }
}
//------------------------------------------------------------------------------------------
// Task_partition
//------------------------------------------------------------------------------------------
void Task_partition()
{
     int  i,j;
     int  temp;
     int  my_task;
     int  length;

     int *temp_i;
     int *temp_j;
     int *temp_num;

     int *arr_row_s;
     int *arr_row_e;
     int *arr_row_item_s;
     int *arr_row_item_e;

     int *arr_task_num;
     int *arr_calcu_s;
     int *arr_calcu_e;
     int *arr_dis_s;
     int *arr_dis_e;


     if (myrank==0)
     {

          printf("---------Task partition start.----------\n");

          length=elem_number_total*36;
          temp_i=(int*) malloc( (length)*sizeof(int) );
          temp_j=(int*) malloc( (length)*sizeof(int) );

          temp_num  =(int*) malloc( (length)*sizeof(int) );

          for (i=0;i<=elem_number_total*36-1;i++)
          {
               temp_num[i]=i;
          }

          temp=0;
          for (i=1;i<=elem_number_total;i++)
          {
               j=(i-1)*36;
               temp_i[j+0] =elem_edge1[i];temp_j[j+0] =elem_edge1[i];
               temp_i[j+1] =elem_edge1[i];temp_j[j+1] =elem_edge2[i];
               temp_i[j+2] =elem_edge1[i];temp_j[j+2] =elem_edge3[i];
               temp_i[j+3] =elem_edge1[i];temp_j[j+3] =elem_edge4[i];
               temp_i[j+4] =elem_edge1[i];temp_j[j+4] =elem_edge5[i];
               temp_i[j+5] =elem_edge1[i];temp_j[j+5] =elem_edge6[i];

               temp_i[j+6] =elem_edge2[i];temp_j[j+6] =elem_edge1[i];
               temp_i[j+7] =elem_edge2[i];temp_j[j+7] =elem_edge2[i];
               temp_i[j+8] =elem_edge2[i];temp_j[j+8] =elem_edge3[i];
               temp_i[j+9] =elem_edge2[i];temp_j[j+9] =elem_edge4[i];
               temp_i[j+10]=elem_edge2[i];temp_j[j+10]=elem_edge5[i];
               temp_i[j+11]=elem_edge2[i];temp_j[j+11]=elem_edge6[i];

               temp_i[j+12]=elem_edge3[i];temp_j[j+12]=elem_edge1[i];
               temp_i[j+13]=elem_edge3[i];temp_j[j+13]=elem_edge2[i];
               temp_i[j+14]=elem_edge3[i];temp_j[j+14]=elem_edge3[i];
               temp_i[j+15]=elem_edge3[i];temp_j[j+15]=elem_edge4[i];
               temp_i[j+16]=elem_edge3[i];temp_j[j+16]=elem_edge5[i];
               temp_i[j+17]=elem_edge3[i];temp_j[j+17]=elem_edge6[i];

               temp_i[j+18]=elem_edge4[i];temp_j[j+18]=elem_edge1[i];
               temp_i[j+19]=elem_edge4[i];temp_j[j+19]=elem_edge2[i];
               temp_i[j+20]=elem_edge4[i];temp_j[j+20]=elem_edge3[i];
               temp_i[j+21]=elem_edge4[i];temp_j[j+21]=elem_edge4[i];
               temp_i[j+22]=elem_edge4[i];temp_j[j+22]=elem_edge5[i];
               temp_i[j+23]=elem_edge4[i];temp_j[j+23]=elem_edge6[i];

               temp_i[j+24]=elem_edge5[i];temp_j[j+24]=elem_edge1[i];
               temp_i[j+25]=elem_edge5[i];temp_j[j+25]=elem_edge2[i];
               temp_i[j+26]=elem_edge5[i];temp_j[j+26]=elem_edge3[i];
               temp_i[j+27]=elem_edge5[i];temp_j[j+27]=elem_edge4[i];
               temp_i[j+28]=elem_edge5[i];temp_j[j+28]=elem_edge5[i];
               temp_i[j+29]=elem_edge5[i];temp_j[j+29]=elem_edge6[i];

               temp_i[j+30]=elem_edge6[i];temp_j[j+30]=elem_edge1[i];
               temp_i[j+31]=elem_edge6[i];temp_j[j+31]=elem_edge2[i];
               temp_i[j+32]=elem_edge6[i];temp_j[j+32]=elem_edge3[i];
               temp_i[j+33]=elem_edge6[i];temp_j[j+33]=elem_edge4[i];
               temp_i[j+34]=elem_edge6[i];temp_j[j+34]=elem_edge5[i];
               temp_i[j+35]=elem_edge6[i];temp_j[j+35]=elem_edge6[i];

          }
          quick_sort_2D(temp_i,temp_j,temp_num,0,length-1);

          /***************matrix_item_average_assign_task*************/
          arr_task_num   =(int*) malloc( numprocs*sizeof(int) );
          arr_calcu_s    =(int*) malloc( numprocs*sizeof(int) );
          arr_calcu_e    =(int*) malloc( numprocs*sizeof(int) );

          if (length%numprocs==0)
          {
               arr_calcu_s[0]=0;
               arr_calcu_e[0]=length/numprocs-1;

               for (i=0;i<=numprocs-1;i++)
                    arr_task_num[i]=length/numprocs;

               for (i=1;i<=numprocs-1;i++)
               {
                    arr_calcu_s[i]=arr_calcu_s[i-1]+arr_task_num[i-1];
                    arr_calcu_e[i]=arr_calcu_e[i-1]+arr_task_num[i];
               }
          }
          else
          {
               arr_calcu_s[0]=0;
               arr_calcu_e[0]=length/numprocs;

               for (i=0;i<=numprocs-1;i++)
               {
                    arr_task_num[i]=length/numprocs;
                    if (i<length%numprocs) arr_task_num[i]=arr_task_num[i]+1;
               }

               for (i=1;i<=numprocs-1;i++)
               {
                    arr_calcu_s[i]=arr_calcu_s[i-1]+arr_task_num[i-1];
                    arr_calcu_e[i]=arr_calcu_e[i-1]+arr_task_num[i];
               }
          }

          /***************matrix_row_assign_task********************/
          arr_row_s =(int*) malloc( numprocs*sizeof(int) );
          arr_row_e =(int*) malloc( numprocs*sizeof(int) );
          arr_row_item_s =(int*) malloc( numprocs*sizeof(int) );
          arr_row_item_e =(int*) malloc( numprocs*sizeof(int) );

          arr_row_s[0]=1;
          arr_row_e[numprocs-1]=temp_i[length-1];

          for (i=0;i<=numprocs-2;i++)
          {
               temp=arr_calcu_e[i];
               if  (temp_i[temp]==temp_i[temp+1])
               {
                    arr_row_e[i]  =temp_i[temp]-1;
                    arr_row_s[i+1]=temp_i[temp];
               }
               if  (temp_i[temp]!=temp_i[temp+1])
               {
                    arr_row_e[i]  =temp_i[temp];
                    arr_row_s[i+1]=temp_i[temp+1];
               }

          }

          arr_row_item_s[0]=0;
          arr_row_item_e[numprocs-1]=length-1;

          temp=0;
          j=0;
          for (i=0;i<=length-2;i++)
          {
               if (numprocs>1)
               {
                    if  (temp_i[i]==arr_row_e[j] && temp_i[i+1]==arr_row_s[j+1])
                    {
                         arr_row_item_e[j]  =temp  ;
                         arr_row_item_s[j+1]=temp+1;
                         j++;
                         if (j==numprocs-1) break;
                    }
               }
               temp++;
          }

          for (i=0;i<=numprocs-1;i++)
          {
               arr_task_num[i]=arr_row_item_e[i]-arr_row_item_s[i]+1;
          }

          /*******************verification************************/
          for (i=0;i<=numprocs-1;i++)
               printf("row_s=%d,row_e=%d\n",arr_row_s[i],arr_row_e[i]);

          for (i=0;i<=numprocs-1;i++)
               printf("row_it_s=%d,row_it_e=%d\n",arr_row_item_s[i],arr_row_item_e[i]);

           
          for (i=0;i<=numprocs-2;i++)
          {     
               if (temp_i[arr_row_item_e[i]]==temp_i[arr_row_item_s[i+1]])
               {
                    printf("Task partition 0 process compute error!\n");
                    exit(-1);
               }
          }  
     }

     /*********************data_distribute******************/

     matrix_item_distribute(arr_row_s,arr_row_e,arr_row_item_s,arr_row_item_e,\
     temp_i,temp_j,temp_num,arr_task_num);

     if (myrank==0)
     {
             free(temp_i);
             free(temp_j);
             free(temp_num);
             free(arr_row_s);
             free(arr_row_e);
             free(arr_calcu_s);
             free(arr_calcu_e);
             free(arr_row_item_s);
             free(arr_row_item_e);
             //free(elem_edge1);
             //free(elem_edge2);
             //free(elem_edge3);
             //free(elem_edge4);
             //free(elem_edge5);
             //free(elem_edge6);
     }

     if (myrank==0)printf("---------Task partion done.-------------\n");
}

//------------------------------------------------------------------------------------------
// malloc_face_struct
//------------------------------------------------------------------------------------------
struct boundary_face malloc_face_struct(struct boundary_face face)
{
     face.bound_face_edge1=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_face_edge2=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_face_edge3=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_face_elem =(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_face_point=(int*)malloc((face.bound_face_number+1)*sizeof(int));
     face.bound_elem_edge1=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_elem_edge2=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_elem_edge3=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_elem_edge4=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_elem_edge5=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     face.bound_elem_edge6=(int*)malloc((face.bound_face_number+1)*sizeof(int)); 
     return (face); 
}
//------------------------------------------------------------------------------------------
// generate_face
//------------------------------------------------------------------------------------------
struct boundary_face generate_face(struct boundary_face face, double min_max, double *node_c)
{
     int i;
     int temp;
     int tempc;
     int n1,n2,n3,n4;
     int node1,node2,node3,node4;


     temp=0;
     for (i=1;i<=elem_number_total;i++)
     {
          node1=elem_node1[i];
          node2=elem_node2[i];
          node3=elem_node3[i];
          node4=elem_node4[i];

          n1=0;n2=0;n3=0;n4=0;
          tempc=0;
          if (node_c[node1]==min_max) {tempc++ ;n1=1;}
          if (node_c[node2]==min_max) {tempc++ ;n2=1;}
          if (node_c[node3]==min_max) {tempc++ ;n3=1;}
          if (node_c[node4]==min_max) {tempc++ ;n4=1;}

          if (tempc==3)
          {
               temp++;
               face.bound_face_elem [temp]=i;
               face.bound_elem_edge1[temp]=elem_edge1[i];
               face.bound_elem_edge2[temp]=elem_edge2[i];
               face.bound_elem_edge3[temp]=elem_edge3[i];
               face.bound_elem_edge4[temp]=elem_edge4[i];
               face.bound_elem_edge5[temp]=elem_edge5[i];
               face.bound_elem_edge6[temp]=elem_edge6[i];
               //printf("bound_edge1=%d\n",face.bound_elem_edge1[temp]);

               if (n1 ==0)
               {
                    face.bound_face_edge1[temp]=elem_edge4[i];
                    face.bound_face_edge2[temp]=elem_edge5[i];
                    face.bound_face_edge3[temp]=elem_edge6[i];
                    face.bound_face_point[temp]=1;
               }
               if (n2 ==0)
               {
                    face.bound_face_edge1[temp]=elem_edge2[i];
                    face.bound_face_edge2[temp]=elem_edge3[i];
                    face.bound_face_edge3[temp]=elem_edge6[i];
                    face.bound_face_point[temp]=2;
               }
               if (n3 ==0)
               {
                    face.bound_face_edge1[temp]=elem_edge1[i];
                    face.bound_face_edge2[temp]=elem_edge3[i];
                    face.bound_face_edge3[temp]=elem_edge5[i];
                    face.bound_face_point[temp]=3;
               }
               if (n4 ==0)
               {
                    face.bound_face_edge1[temp]=elem_edge1[i];
                    face.bound_face_edge2[temp]=elem_edge2[i];
                    face.bound_face_edge3[temp]=elem_edge4[i];
                    face.bound_face_point[temp]=4;
               }

          }

     }

     return (face); 
}
//------------------------------------------------------------------------------------------
// Find_collect
//------------------------------------------------------------------------------------------
void Find_collect()
{
     int temp;
     int i,j,k;
     int node1,node2,node3,node4;
     double x0,x1,x2,x3,x4;
     double y0,y1,y2,y3,y4;
     double z0,z1,z2,z3,z4;
     double V0,V1,V2,V3,V4;
     int surface_num;
     int *surface_elem;
     int **l;
     int *surface_start;
     int *surface_end;
     int *recv;

     if (myrank==0)
     {

          /********************************************/
          temp=0;
          l =(int**)malloc((elem_number_total+1)*sizeof(int*));

          for (i=0;i<=elem_number_total;i++)
               l[i] =(int*)malloc((collect_num+1)*sizeof(int));

          for (i=0;i<=elem_number_total;i++)
              for (j=0;j<=collect_num;j++)
                   l[i][j]=0;

          for (i=1;i<=elem_number_total;i++)
          {
               node1=elem_node1[i];
               node2=elem_node2[i];
               node3=elem_node3[i];
               node4=elem_node4[i];
                
               x1=node_x[node1];
               x2=node_x[node2];
               x3=node_x[node3];
               x4=node_x[node4];
 
               y1=node_y[node1];
               y2=node_y[node2];
               y3=node_y[node3];
               y4=node_y[node4];

               z1=node_z[node1];
               z2=node_z[node2];
               z3=node_z[node3];
               z4=node_z[node4];

               if (z1>0 ||z2>0 ||z3>0 ||z4 >0) continue;

               V0= x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2\
                 - x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1\
                 + x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1\
                 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1;

               for (j=1;j<=collect_num;j++)
               {
                    x0=collect_x[j];
                    y0=collect_y[j];
                    z0=collect_z[j];

                    V1= x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2\
                      - x0*y3*z4 + x0*y4*z3 + x3*y0*z4 - x3*y4*z0 - x4*y0*z3 + x4*y3*z0\
                      + x0*y2*z4 - x0*y4*z2 - x2*y0*z4 + x2*y4*z0 + x4*y0*z2 - x4*y2*z0\
                      - x0*y2*z3 + x0*y3*z2 + x2*y0*z3 - x2*y3*z0 - x3*y0*z2 + x3*y2*z0;
 
                    V2= x0*y3*z4 - x0*y4*z3 - x3*y0*z4 + x3*y4*z0 + x4*y0*z3 - x4*y3*z0\
                      - x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1\
                      + x1*y0*z4 - x1*y4*z0 - x0*y1*z4 + x0*y4*z1 + x4*y1*z0 - x4*y0*z1\
                      - x1*y0*z3 + x1*y3*z0 + x0*y1*z3 - x0*y3*z1 - x3*y1*z0 + x3*y0*z1;

                    V3= x2*y0*z4 - x2*y4*z0 - x0*y2*z4 + x0*y4*z2 + x4*y2*z0 - x4*y0*z2\
                      - x1*y0*z4 + x1*y4*z0 + x0*y1*z4 - x0*y4*z1 - x4*y1*z0 + x4*y0*z1\
                      + x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1\
                      - x1*y2*z0 + x1*y0*z2 + x2*y1*z0 - x2*y0*z1 - x0*y1*z2 + x0*y2*z1;

                    V4= x2*y3*z0 - x2*y0*z3 - x3*y2*z0 + x3*y0*z2 + x0*y2*z3 - x0*y3*z2\
                      - x1*y3*z0 + x1*y0*z3 + x3*y1*z0 - x3*y0*z1 - x0*y1*z3 + x0*y3*z1\
                      + x1*y2*z0 - x1*y0*z2 - x2*y1*z0 + x2*y0*z1 + x0*y1*z2 - x0*y2*z1\
                      - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1;
                    
                    if (V0*V1>=0 && V0*V2>=0 && V0*V3>=0 && V0*V4>=0)
                    {
                         l[i][j]=i;
                         temp++; 
                    }
               }
          }

          col_elem_num=temp;
          col_elem=(int*)malloc((col_elem_num)*sizeof(int)); 
          col_num =(int*)malloc((col_elem_num)*sizeof(int));

          k=0;
          for (i=1;i<=elem_number_total;i++)
          {
               for (j=1;j<=collect_num;j++)
               {
                    if (l[i][j]!=0)
                    {
                        col_elem[k]=l[i][j];
                        col_num [k]=j;
                        k++;
                    }
               }
          }
          
          for (i=0;i<=elem_number_total;i++)
               free(l[i]);
          free(l);

          /*******************verification************************/

          for (i=0;i<=col_elem_num-1;i++)
          {
               printf("col_elem=%d,col_num=%d\n",col_elem[i],col_num[i]);
               j=col_num[i];
               printf("col_x=%6.6lf,col_y=%6.6lf,col_z=%6.6lf\n",\
               collect_x[j],collect_y[j],collect_z[j]);
               j=elem_node1[col_elem[i]];
               printf("elem_node1=%d,x=%6.6lf,y=%6.6lf,z=%6.6lf\n",j,\
               node_x[j],node_y[j],node_z[j]);
               j=elem_node2[col_elem[i]];
               printf("elem_node2=%d,x=%6.6lf,y=%6.6lf,z=%6.6lf\n",j,\
               node_x[j],node_y[j],node_z[j]);
               j=elem_node3[col_elem[i]];
               printf("elem_node3=%d,x=%6.6lf,y=%6.6lf,z=%6.6lf\n",j,\
               node_x[j],node_y[j],node_z[j]);
               j=elem_node4[col_elem[i]];
               printf("elem_node4=%d,x=%6.6lf,y=%6.6lf,z=%6.6lf\n",j,\
               node_x[j],node_y[j],node_z[j]);

          }
                  
          printf("-------Find collect point done.---------\n");
          
          
     }
        
}
//------------------------------------------------------------------------------------------
// verify_face
//------------------------------------------------------------------------------------------
void verify_face(struct boundary_face face, double min_max, double *node_c)
{
     int i,j;
     int temp;
     int tempc;
     int n1,n2,n3,n4;
     int node1,node2,node3,node4;

     for (i=1;i<=face.bound_face_number;i++)
     {
          j=face.bound_face_elem[i];
          node1=elem_node1[j];
          node2=elem_node2[j];
          node3=elem_node3[j];
          node4=elem_node4[j];

          tempc=0;
          if (node_c[node1]==min_max) {tempc++ ;}
          if (node_c[node2]==min_max) {tempc++ ;}
          if (node_c[node3]==min_max) {tempc++ ;}
          if (node_c[node4]==min_max) {tempc++ ;}

          if (tempc!=3)
          {
               printf("Bounday face compute error!\n");
               exit(-1);
          }

     }

}
//------------------------------------------------------------------------------------------
// quick_sort_2D
//------------------------------------------------------------------------------------------
void quick_sort_2D(int *e1, int *e2, int *n1, int low, int high)
{

     int privotLoc;
     if(low < high)
     {  
          privotLoc = partition_int(e1,e2,n1,low,high);

	  quick_sort_2D(e1,e2,n1,low,privotLoc -1); 
	  quick_sort_2D(e1,e2,n1,privotLoc + 1, high);       
     }  
}

int partition_int(int *e1, int *e2, int *n1, int low, int high) 
{  
     int privotKey1 = e1[low];
     int privotKey2 = e2[low];

     while(low < high)
     {

          while(low < high  && \
	  (e1[high] >  privotKey1 || (e1[high] == privotKey1 && e2[high] >= privotKey2)))\
	  --high;
	  {     
	       swap(&e1[low], &e1[high]);
	       swap(&e2[low], &e2[high]);
	       swap(&n1[low], &n1[high]);             
	  }
	  while(low < high  && \
	  (e1[low]  <  privotKey1 || (e1[low]  == privotKey1 && e2[low] <= privotKey2)))\
	  ++low;
	  {      
	       swap(&e1[low], &e1[high]);
	       swap(&e2[low], &e2[high]);
	       swap(&n1[low], &n1[high]);             
          }        
     }  
     return low; 
}  

void swap(int *a, int *b)  
{  
     int tmp = *a;  
     *a = *b;  
     *b = tmp;  
}  
//------------------------------------------------------------------------------------------
// max_arry,min_arry
//------------------------------------------------------------------------------------------
int max(int a, int b)
{
     if(a>b) return a;
     else    return b;
}

int min(int a, int b)
{
     if(a<b) return a;
     else    return b;
}

double min_arry(double *arry, int length)
{
     int    i;
     double temp;
     temp = arry[1];

     for (i=2;i<=length;i++)
     {   
          //printf("arry=%f",arry[i]);
	  if(temp>arry[i]) temp=arry[i];
     }  
     return temp;
}

double max_arry(double *arry, int length)
{
     int    i;
     double temp;
     temp = arry[1];

     for (i=2;i<=length;i++)
     {   
          //printf("arry=%f",arry[i]);
          if(temp<arry[i]) temp=arry[i];
     }  
     return temp;
}
//------------------------------------------------------------------------------------------
