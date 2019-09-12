// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ Solver.c
#include "solver.h"

//-------------------------------------------------------------------------------------------
// 1.Main_Solver();
//-------------------------------------------------------------------------------------------
void Main_Solver()
{
     superlu_dist_options_t options;
     SuperLUStat_t stat;
     SuperMatrix A;
     ScalePermstruct_t ScalePermstruct;
     LUstruct_t LUstruct;
     SOLVEstruct_t SOLVEstruct;
     gridinfo_t grid;
     double   *berr;
     doublecomplex   *xtrue;
     doublecomplex   *nzval_loc;

     int v_major, v_minor, v_bugfix;
     int iam,info,nrhs,ldb;
     int m,n;

     int ifre; 
     int ifres,ifree,ifretemp;    
     int i,j,k;
     int boundary_start;
     int boundary_end;
     double E0_r_x;
     double E0_i_x;
     double E0_r_y;
     double E0_i_y;

     int    *temp_matrix_asub;
     int    *temp_matrix_xa ;   
     double *time_fre;
     double time_fres,time_free;

     MPI_Barrier(MPI_COMM_WORLD);

     if (myrank==0) time1=MPI_Wtime();


     /*********calculate frequency range each process*************/
     ifretemp=fre_num/MPI_FRE;
     if (fre_num%MPI_FRE==0)
     {
          ifres=myrank%MPI_FRE*ifretemp+1;
          ifree=myrank%MPI_FRE*ifretemp+ifretemp;
     }
     else if (fre_num%MPI_FRE!=0) 
     {
          if(myrank%MPI_FRE<fre_num%MPI_FRE)
          {
               ifres=myrank%MPI_FRE*(ifretemp+1)+1;
               ifree=myrank%MPI_FRE*(ifretemp+1)+ifretemp+1; 
          }
          else if(myrank%MPI_FRE>=fre_num%MPI_FRE)
          {
               ifres=myrank%MPI_FRE*ifretemp+1+fre_num%MPI_FRE;
               ifree=myrank%MPI_FRE*ifretemp+ifretemp+fre_num%MPI_FRE;
          }
     }
     printf ("myrank=%d,ifres=%d,ifree=%d\n",myrank,ifres,ifree);

     /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID. 
     ------------------------------------------------------------*/

     superlu_gridinit(NewWorld_comp, LU_nprow, LU_npcol, &grid);

     iam = grid.iam;

     printf ("myrank=%d,iam=%d\n",myrank,iam);

     if (myrank==0)
     {
          printf("****************************************\n");
          printf("             Solver_start               \n");
          printf("****************************************\n");
           
          printf("SuperLU process rows    %d\n", LU_nprow);
          printf("SuperLU process columns %d\n", LU_npcol);
          superlu_dist_GetVersionNumber(&v_major, &v_minor, &v_bugfix);
          printf("Library version:\t%d.%d.%d\n", v_major, v_minor, v_bugfix);
     }

     b_x  = (doublecomplex *) doublecomplexMalloc_dist(m_loc);
     b_y  = (doublecomplex *) doublecomplexMalloc_dist(m_loc);

     if (iam==0)
     {
          b_x_total =(struct doublecomplex_z*)malloc(edge_number_total*sizeof(struct doublecomplex_z)); 
          b_y_total =(struct doublecomplex_z*)malloc(edge_number_total*sizeof(struct doublecomplex_z));
     }

     /* Number of right-hand side. */ 
     ldb  = m_loc;
     nrhs = 1;

     berr = doubleMalloc_dist(nrhs); 

     /************************frequence_loop**************************/

     //ifre=1;
 
     MPI_Barrier(MPI_COMM_WORLD); 
     time_fre         =(double*)  malloc( (ifree-ifres+1)*sizeof(double) );         
  
    
     for (ifre=ifres;ifre<=ifree;ifre++)
     {
          MPI_Barrier(NewWorld_comp); 
          if (mynewrank==0) 
          {
                time_fres=MPI_Wtime();
                printf("--------Solver_Frequence=%lf------------\n",frequency[ifre]);
          }
          omega=2*pi*frequency[ifre];

          nzval_loc = (doublecomplex *) doublecomplexMalloc_dist(nnz_loc);
          temp_matrix_xa   =(int*) malloc( (m_loc+1)*sizeof(int) );
          temp_matrix_asub =(int*) malloc( (matrix_number)*sizeof(int) );

          /*create new Matrix */
          for (i=0;i<=matrix_number-1;i++)
               temp_matrix_asub[i]=matrix_asub[i];
          for (i=0;i<=m_loc;i++)
               temp_matrix_xa[i]=matrix_xa[i];

          /*real partion*/          
          for (i=0;i<=matrix_number-1;i++)
               nzval_loc[i].r=matrix_a_r[i];

          /*imag partion change according to omega*/
          for (i=0;i<=matrix_number-1;i++)
          { 
               nzval_loc[i].i=matrix_a_i[i]*omega;
          }

          /*set right hand b*/
          if (bound_mode==1)
          {    b_generate_Dirichlet();}
          else if (bound_mode==2)
          {    b_generate_Newman();}
          else
          {   
               printf( "Boundary condition mode set error !\n");
               exit(0);
          }


          zCreate_CompRowLoc_Matrix_dist(&A, n_loc, n_loc, nnz_loc, m_loc,\
          fst_row, nzval_loc, temp_matrix_asub, temp_matrix_xa , SLU_NR_loc, SLU_Z, SLU_GE);

          m = A.nrow;
          n = A.ncol;

          if (ifre==ifres)
          {
             /*Set the default input options:
               options.Fact = DOFACT;
               options.Equil = YES;
               options.ColPx_item_locerm = METIS_AT_PLUS_A;
               options.RowPerm = LargeDiag_MC64;
               options.ReplaceTinyPivot = NO;
               options.Trans = NOTRANS;
               options.IterRefine = DOUBLE;
               options.SolveInitialized = NO;
               options.RefineInitialized = NO;
               options.PrintStat = YES;
             */ 
               set_default_options_dist(&options);
               ScalePermstructInit(m, n, &ScalePermstruct);
               LUstructInit(n, &LUstruct);
          }
          else 
          {          
               options.Fact = SamePattern;
          }

          if (!iam)//iam==0; 
          {
	       print_sp_ienv_dist(&options);
	       print_options_dist(&options);
	       fflush(stdout);
          }  

          /* Initialize the statistics variables. */
          PStatInit(&stat);

          /* Call the linear equation solver Ex mode. */
          pzgssvx(&options, &A, &ScalePermstruct, b_x, ldb, nrhs, &grid,\
	  &LUstruct, &SOLVEstruct, berr, &stat, &info);

          //for (i=1;i<=m_loc;i++)
          //     printf("i=%d,br=%lf,bi=%lf\n",i,b_x[i-1].r,b_x[i-1].i);

          /* Print the statistics. */
          PStatPrint(&options, &stat, &grid); 
          PStatFree(&stat);

          /* Initialize the statistics variables. */
          options.Fact = FACTORED;
          PStatInit(&stat); 

          if (!iam)//iam==0; 
          {
	       print_sp_ienv_dist(&options);
	       print_options_dist(&options);
	       fflush(stdout);
          }  

          /* Call the linear equation solver Ey mode. */
          pzgssvx(&options, &A, &ScalePermstruct, b_y, ldb, nrhs, &grid,\
	  &LUstruct, &SOLVEstruct, berr, &stat, &info);

          //for (i=1;i<=m_loc;i++)
          //     printf("i=%d,br=%lf,bi=%lf\n",i,b_y[i-1].r,b_y[i-1].i);
         
          /* Print the statistics. */           
          PStatPrint(&options, &stat, &grid); 
          PStatFree(&stat);

          Destroy_CompRowLoc_Matrix_dist(&A);          
          Destroy_LU(n, &grid, &LUstruct);

          /* Gather b_x,b_y to main processes for calculating */           
          gather_b();

          //for (i=1;i<=edge_number_total;i++)
          //     printf("i=%d,br=%lf,bi=%lf\n",i,b_x_total[i-1].r,b_x_total[i-1].i);
          post_prepocess(ifre);

          MPI_Barrier(NewWorld_comp); 
          if (mynewrank==0)
          {
              time_free=MPI_Wtime();
              time_fre[ifre-ifres]=time_free-time_fres;
          }
     }
     
     

     MPI_Barrier(MPI_COMM_WORLD);
     if (myrank==0)
     {
          time2=MPI_Wtime();
          time_e=time2;
          Solver_time=time2-time1;
          Total_time=time_e-time_s;

          printf("Prepocess_time=%lf\n",Prepocess_time);
          printf("Fem_compute_time=%lf\n",Fem_compute_time);
          printf("Solver_time=%lf\n",Solver_time);
          printf("Total_time=%lf\n",Total_time);

     }
     
     if (mynewrank==0)
     {
         for (i=0;i<=ifree-ifres;i++)
              printf("Fre=%lf,Solver_time=%lf\n",frequency[i+ifres],time_fre[i]);
     } 

     ScalePermstructFree(&ScalePermstruct);  
     LUstructFree(&LUstruct);

     if ( options.SolveInitialized )
     {
          zSolveFinalize(&options, &SOLVEstruct);
     }

     free(matrix_asub);
     free(matrix_xa);
     free(matrix_a_r);
     free(matrix_a_i);

     SUPERLU_FREE(b_x);
     SUPERLU_FREE(b_y);
     //SUPERLU_FREE(berr);
     //SUPERLU_FREE(nzval_loc);

     superlu_gridexit(&grid);
     MPI_Finalize();
}
//-------------------------------------------------------------------------------------------
// right_rand:b_generate();
// 1.Dirichlet 
// 2.Newman
//-------------------------------------------------------------------------------------------
void b_generate_Dirichlet()
{
     int i;
     int bound_num;
     int b_num;
     complex double E0;
     complex double Ex;
     complex double Ey;
     complex double E_temp;

     /*right hand set 0 each frequence*/
     for (i=0;i<=m_loc-1;i++)
     {
          b_x[i].r=0.0;
          b_x[i].i=0.0;
          b_y[i].r=0.0;
          b_y[i].i=0.0;
     }
     
     if (num_out_column!=0)
     {           
          for (i=0;i<=num_out_column-1;i++)
          {               
               bound_num=bound_edge[matrix_index_rh[i]];
               Boundary_E0(bound_num,&E0);
               
               field_projection(bound_num,&E0,&Ex,&Ey);
              
               b_num=matrix_i_rh[i]-fst_row-1;

               E_temp=Ex* ( matrix_a_r_rh[i]  + matrix_a_i_rh[i]* omega* _Complex_I);

               b_x[b_num].r=b_x[b_num].r-creal(E_temp);
               b_x[b_num].i=b_x[b_num].i-cimag(E_temp);

               E_temp=Ey* ( matrix_a_r_rh[i]  + matrix_a_i_rh[i]* omega* _Complex_I);

               b_y[b_num].r=b_y[b_num].r-creal(E_temp);
               b_y[b_num].i=b_y[b_num].i-cimag(E_temp);

          } 
     }

     MPI_Barrier(MPI_COMM_WORLD);

     if (myrank==0)
     {
          printf("-------Boundary_Dirichlet_Done----------\n");
     }
          
}

//-------------------------------------------------------------------------------------------
void b_generate_Newman()
{
     int i;
     int bound_num;
     int b_num;

     /*right hand set 0 each frequence*/
     for (i=0;i<=m_loc-1;i++)
     {
          b_x[i].r=0.0;
          b_x[i].i=0.0;
          b_y[i].r=0.0;
          b_y[i].i=0.0;
     }

     generate_condition(face_up);
     generate_condition(face_down);
     generate_condition(face_left);
     generate_condition(face_right);
     generate_condition(face_forward);
     generate_condition(face_behind);

}
//------------------------------------------------------------------------------------------
void  generate_condition(struct boundary_face face)
{

     int ind_1,ind_2,ind_3,ind_4;
     int i,j,k;
     double x1,x2,x3,x4;
     double y1,y2,y3,y4;
     double z1,z2,z3,z4;
     double Se,Ve;
     double a1,a2,a3,a4;
     double b1,b2,b3,b4;
     double c1,c2,c3,c4;
     double d1,d2,d3,d4;
     double l1,l2,l3,l4,l5,l6;
     double depth1,depth2,depth3,depth4,depth5,depth6;
     int    mode;
     int    edge1,edge2,edge3,edge4,edge5,edge6;
     complex double H0_1,H0_2,H0_3;
     complex double H0_4,H0_5,H0_6;
     complex double H0_x,H0_y;
     complex double H0_b1,H0_b2,H0_b3;

     double N1,N2,N3,N4,N5,N6;

     for (i=1;i<=face.bound_face_number;i++)
     {
          j= face.bound_face_elem[i];
          mode = face.bound_face_mode;

          if (face.bound_face_edge1[i]<(fst_row+1) && face.bound_face_edge1[i]>(fst_row+m_loc+1)\
          &&  face.bound_face_edge2[i]<(fst_row+1) && face.bound_face_edge2[i]>(fst_row+m_loc+1)\
          &&  face.bound_face_edge3[i]<(fst_row+1) && face.bound_face_edge3[i]>(fst_row+m_loc+1))
          {
               continue;
          }

          ind_1=elem_node1[j];
          ind_2=elem_node2[j];
          ind_3=elem_node3[j];
          ind_4=elem_node4[j];
                    
          x1=node_x[ind_1];
          x2=node_x[ind_2];
          x3=node_x[ind_3];
          x4=node_x[ind_4];

          y1=node_y[ind_1];
          y2=node_y[ind_2];
          y3=node_y[ind_3];
          y4=node_y[ind_4];

          z1=node_z[ind_1];
          z2=node_z[ind_2];
          z3=node_z[ind_3];
          z4=node_z[ind_4];

          Ve= (x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2\
            - x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1\
            + x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1\
            - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1)/6.0;

          b1 =-y3*z4 + y4*z3 + y2*z4 - y4*z2 - y2*z3 + y3*z2;
          b2 = y3*z4 - y4*z3 - y1*z4 + y4*z1 + y1*z3 - y3*z1;
          b3 =-y2*z4 + y4*z2 + y1*z4 - y4*z1 - y1*z2 + y2*z1;
          b4 = y2*z3 - y3*z2 - y1*z3 + y3*z1 + y1*z2 - y2*z1;

          c1 =-x2*z4 + x2*z3 + x3*z4 - x3*z2 - x4*z3 + x4*z2;
          c2 = x1*z4 - x1*z3 - x3*z4 + x3*z1 + x4*z3 - x4*z1;
          c3 =-x1*z4 + x1*z2 + x2*z4 - x2*z1 - x4*z2 + x4*z1;
          c4 = x1*z3 - x1*z2 - x2*z3 + x2*z1 + x3*z2 - x3*z1;

          d1 =-x2*y3 + x2*y4 + x3*y2 - x3*y4 - x4*y2 + x4*y3;
          d2 = x1*y3 - x1*y4 - x3*y1 + x3*y4 + x4*y1 - x4*y3;
          d3 =-x1*y2 + x1*y4 + x2*y1 - x2*y4 - x4*y1 + x4*y2;
          d4 = x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2;

          l1 = sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
          l2 = sqrt(pow(x1-x3,2)+pow(y1-y3,2)+pow(z1-z3,2));
          l3 = sqrt(pow(x1-x4,2)+pow(y1-y4,2)+pow(z1-z4,2));
          l4 = sqrt(pow(x2-x3,2)+pow(y2-y3,2)+pow(z2-z3,2));
          l5 = sqrt(pow(x2-x4,2)+pow(y2-y4,2)+pow(z2-z4,2));
          l6 = sqrt(pow(x3-x4,2)+pow(y3-y4,2)+pow(z3-z4,2));

          edge1=face.bound_elem_edge1[i];
          edge2=face.bound_elem_edge2[i];
          edge3=face.bound_elem_edge3[i];
          edge4=face.bound_elem_edge4[i];
          edge5=face.bound_elem_edge5[i];
          edge6=face.bound_elem_edge6[i];

          //printf("face.bound_elem_edge1[i]=%d\n",face.bound_elem_edge1[i]);
          MPI_Barrier(MPI_COMM_WORLD);

          H0_1=0;H0_2=0;
          H0_3=0;H0_4=0;
          H0_5=0;H0_6=0;

          Boundary_H0(edge1, &H0_1);
          Boundary_H0(edge2, &H0_2);
          Boundary_H0(edge3, &H0_3);
          Boundary_H0(edge4, &H0_4);
          Boundary_H0(edge5, &H0_5);
          Boundary_H0(edge6, &H0_6);

          N1=l1/6/Ve;
          N2=l2/6/Ve;
          N3=l3/6/Ve;
          N4=l4/6/Ve;
          N5=l5/6/Ve;
          N6=l6/6/Ve;
          
          if (face.bound_face_point[i]==1)
          {
               Se=sqrt((pow((x2-x3),2)+pow((y2-y3),2)+pow((z2-z3),2))*\
               (pow((x2-x4),2)+pow((y2-y4),2)+pow((z2-z4),2))-\
               pow((x2-x3)*(x2-x4)+(y2-y3)*(y2-y4)+(z2-z3)*(z2-z4),2))/2;

             
               N4=N4*Se/3;
               N5=N5*Se/3;
               N6=N6*Se/3;
               H0_b1=N4*H0_4;
               H0_b2=N5*H0_5;
               H0_b3=N6*H0_6;

               if (ind_2>ind_3) {H0_b1=H0_b1*(-1.0);}
               if (ind_4>ind_2) {H0_b2=H0_b2*(-1.0);}
               if (ind_3>ind_4) {H0_b3=H0_b3*(-1.0);}

               p1_generate(mode,H0_b1,H0_b2,H0_b3,b2,b3,b4,c2,c3,c4,d2,d3,d4,i,face);
          }  
          else if (face.bound_face_point[i]==2)
          {
               Se=sqrt((pow((x1-x3),2)+pow((y1-y3),2)+pow((z1-z3),2))*\
                  (pow((x1-x4),2)+pow((y1-y4),2)+pow((z1-z4),2))-\
                  pow((x1-x3)*(x1-x4)+(y1-y3)*(y1-y4)+(z1-z3)*(z1-z4),2))/2;

               N2=N2*Se/3;
               N3=N3*Se/3;
               N6=N6*Se/3;

               H0_b1=N2*H0_2;
               H0_b2=N3*H0_3;
               H0_b3=N6*H0_6;

               if (ind_1>ind_3) {H0_b1=H0_b1*(-1.0);}
               if (ind_1>ind_4) {H0_b2=H0_b2*(-1.0);}
               if (ind_3>ind_4) {H0_b3=H0_b3*(-1.0);}

               p2_generate(mode,H0_b1,H0_b2,H0_b3,b1,b3,b4,c1,c3,c4,d1,d3,d4,i,face);
                       
          }
          else if (face.bound_face_point[i]==3)
          {
               Se=sqrt((pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2))*\
                  (pow((x1-x4),2)+pow((y1-y4),2)+pow((z1-z4),2))-\
                  pow((x1-x2)*(x1-x4)+(y1-y2)*(y1-y4)+(z1-z2)*(z1-z4),2))/2;

               N1=N1*Se/3;
               N3=N3*Se/3;
               N5=N5*Se/3;

               H0_b1=N1*H0_1;
               H0_b2=N3*H0_3;
               H0_b3=N5*H0_5;

               if (ind_1>ind_2) {H0_b1=H0_b1*(-1.0);}
               if (ind_1>ind_4) {H0_b2=H0_b2*(-1.0);}
               if (ind_4>ind_2) {H0_b3=H0_b3*(-1.0);}

               p3_generate(mode,H0_b1,H0_b2,H0_b3,b1,b2,b4,c1,c2,c4,d1,d2,d4,i,face);
   
          }
          else if (face.bound_face_point[i]==4)
          {

               Se=sqrt((pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2))*\
                  (pow((x1-x3),2)+pow((y1-y3),2)+pow((z1-z3),2))-\
                  pow((x1-x2)*(x1-x3)+(y1-y2)*(y1-y3)+(z1-z2)*(z1-z3),2))/2;

               N1=N1*Se/3;
               N2=N2*Se/3;
               N4=N4*Se/3;

               H0_b1=N1*H0_1;
               H0_b2=N2*H0_2;
               H0_b3=N4*H0_3;

               if (ind_1>ind_2) {H0_b1=H0_b1*(-1.0);}
               if (ind_1>ind_3) {H0_b2=H0_b2*(-1.0);}
               if (ind_2>ind_3) {H0_b3=H0_b3*(-1.0);}

               p4_generate(mode,H0_b1,H0_b2,H0_b3,b1,b2,b3,c1,c2,c3,d1,d2,d3,i,face);

          }
          else
          {
               printf("bound_face_point error !\n");
               exit(0);
  
          }
     }                         
}
//-------------------------------------------------------------------------------------------
// Background_fied_E0
//-------------------------------------------------------------------------------------------
void Boundary_E0(int edge_ind,complex double *E0)
{
     //only can calculate 1 layer now;

     int node1,node2;
     double depth;

     complex double k0;
     complex double k1;
     complex double K;
     complex double km,km2;
     complex double A0,B0,A1;

     node1=edge_node1[edge_ind];
     node2=edge_node2[edge_ind];

     depth=(node_z[node1]+node_z[node2])/2.0;
     
     //1D-analytical
     k0   = csqrt((-1.0)*mu0* layer_cond[1] * omega* _Complex_I);
     k1   = csqrt((-1.0)*mu0* layer_cond[2] * omega* _Complex_I);

     K =(k0-k1)/(k0+k1);
     km=K*cexp((-2.0)*k0*layer_depth[1]);

     A0=1.0/(1+km);
     B0=km*A0;
     A1=(A0*cexp((-1.0)*k0*layer_depth[1])+B0*cexp(k0*layer_depth[1]))/cexp((-1.0)*k1*layer_depth[1]);


     if (depth>0)
     {
          *E0  = A0*cexp((-1.0)*k0*(layer_depth[1]-depth))+\
                 B0*cexp(k0*(layer_depth[1]-depth));
     }
     else if(depth<=0)
     {
          *E0   =A1*cexp((-1.0)*k1*(layer_depth[1]+abs(depth)));
     }
}
//-------------------------------------------------------------------------------------------
// Background_fied_H0
//-------------------------------------------------------------------------------------------
void Boundary_H0(int edge_ind,complex double *H0)
{
     //only can calculate 1 layer now;
     int node1,node2;
     double depth;

     complex double k0;
     complex double k1;
     complex double K;
     complex double km,km2;
     complex double A0,B0,A1;

     node1=edge_node1[edge_ind];
     node2=edge_node2[edge_ind];

     depth=(node_z[node1]+node_z[node2])/2.0;

     //1D-analytical
     k0   = csqrt((-1.0)*mu0* layer_cond[1]* omega* _Complex_I);
     k1   = csqrt((-1.0)*mu0* layer_cond[2]* omega* _Complex_I);

     K =(k0-k1)/(k0+k1);
     km=K*cexp((-2.0)*k0*layer_depth[1]);

     A0=1.0/(1+km);
     B0=km*A0;
     A1=(A0*cexp((-1.0)*k0*layer_depth[1])+B0*cexp(k0*layer_depth[1]))/cexp((-1.0)*k1*layer_depth[1]);

     if (depth>0)//H0
     {
          *H0  = k0*(A0*cexp((-1.0)*k0*(layer_depth[1]-depth))-\
                 B0*cexp(k0*(layer_depth[1]-depth)));
     }
     else if(depth<=0)
     {
          *H0   =k1*A1*cexp((-1.0)*k1*(layer_depth[1]+abs(depth)));
     }

}
//------------------------------------------------------------------------------------------
// field_projection
//------------------------------------------------------------------------------------------
void field_projection(int edge_ind,complex double *E0,complex double *Ex,complex double *Ey)
{    
     int ind_1,ind_2;
     double x1,x2;
     double y1,y2;
     double z1,z2;
     double l;
     double ratio1,ratio2;
     
     ind_1=edge_node1[edge_ind];
     ind_2=edge_node2[edge_ind];
     
     x1=node_x[ind_1];
     x2=node_x[ind_2];
     y1=node_y[ind_1];
     y2=node_y[ind_2];
     z1=node_z[ind_1];
     z2=node_z[ind_2];
     
     l=sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
     
     ratio1=(x2-x1)/l;
     ratio2=(y2-y1)/l;
      
     *Ex = (*E0)*ratio1 ;
     
     *Ey = (*E0)*ratio2 ;

}

//-------------------------------------------------------------------------------------------
void p1_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b2,double b3,double b4,double c2,double c3 ,double c4,double d2,double d3,double d4,\
int i,struct boundary_face face)
{     

     if (mode==1)
     {
          
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
   
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c3-c2));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c3-c2));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b3-b2)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b3-b2)*(-1.0));

          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c2-c4));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c2-c4));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b2-b4)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b2-b4)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c4-c3));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c4-c3));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b4-b3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b4-b3)*(-1.0));
          }
     }
     else if (mode==2)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c3-c2)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c3-c2)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b3-b2));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b3-b2));
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c2-c4)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c2-c4)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b2-b4));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b2-b4));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c4-c3)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c4-c3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b4-b3));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b4-b3));
          }
     }
     else if (mode==3)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {   
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d2)*(-1.0)); 
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d2)*(-1.0));
          }   
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d2-d4)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d2-d4)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3)*(-1.0));
          }
     }
     else if (mode==4)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d2));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d2));
          
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d2-d4));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d2-d4));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3));
          }
     }
     else if (mode==5)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d2)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d2)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d2-d4)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d2-d4)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }
     else if (mode==6)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d2));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d2));         
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d2-d4));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d2-d4));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {  
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }

}

//-------------------------------------------------------------------------------------------
void p2_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b1,double b3,double b4,double c1,double c3 ,double c4,double d1,double d3,double d4,\
int i,struct boundary_face face)
{

     if (mode==1)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c3-c1));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c3-c1));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b3-b1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b3-b1)*(-1.0));

          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c4-c1));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c4-c1));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b4-b1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b4-b1)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c4-c3));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c4-c3));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b4-b3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b4-b3)*(-1.0));
          }
     }
     else if (mode==2)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c3-c1)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c3-c1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b3-b1));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b3-b1));
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c4-c1)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c4-c1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b4-b1));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b4-b1));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c4-c3)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c4-c3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b4-b3));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b4-b3));
          }
     }
     else if (mode==3)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {   
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d1)*(-1.0)); 
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d1)*(-1.0));
          }   
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3)*(-1.0));
          }
     }
     else if (mode==4)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d1));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d1));
          
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3));
          }
     }
     else if (mode==5)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d1)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }
     else if (mode==6)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d3-d1));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d3-d1));         
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {  
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d4-d3));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d4-d3));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }

}
//-------------------------------------------------------------------------------------------
void p3_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b1,double b2,double b4,double c1,double c2,double c4,double d1,double d2,double d4,\
int i,struct boundary_face face)
{

     if (mode==1)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c2-c1));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c2-c1));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b2-b1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b2-b1)*(-1.0));

          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c4-c1));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c4-c1));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b4-b1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b4-b1)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c2-c4));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c2-c4));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b2-b4)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b2-b4)*(-1.0));
          }
     }
     else if (mode==2)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c2-c1)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c2-c1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b2-b1));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b2-b1));
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c4-c1)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c4-c1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b4-b1));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b4-b1));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c2-c4)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c2-c4)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b2-b4));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b2-b4));
          }
     }
     else if (mode==3)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {   
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1)*(-1.0)); 
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1)*(-1.0));
          }   
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d2-d4)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d2-d4)*(-1.0));
          }
     }
     else if (mode==4)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1));
          
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d2-d4));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d2-d4));
          }
     }
     else if (mode==5)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d2-d4)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d2-d4)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }
     else if (mode==6)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1));         
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d4-d1));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d4-d1));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {  
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d2-d4));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d2-d4));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }

}
//-------------------------------------------------------------------------------------------
void p4_generate (int mode,complex double H0_b1,complex double H0_b2,complex double H0_b3,\
double b1,double b2,double b3,double c1,double c2,double c3,double d1,double d2,double d3,\
int i,struct boundary_face face)
{

     if (mode==1)
     {
         
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c2-c1));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c2-c1));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b2-b1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b2-b1)*(-1.0));

          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c3-c1));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c3-c1));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b3-b1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b3-b1)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c3-c2));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c3-c2));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b3-b2)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b3-b2)*(-1.0));
          }
     }
     else if (mode==2)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(c2-c1)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(c2-c1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(b2-b1));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(b2-b1));
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(c3-c1)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(c3-c1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(b3-b1));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(b3-b1));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(c3-c2)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(c3-c2)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(b3-b2));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(b3-b2));
          }
     }
     else if (mode==3)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {   
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1)*(-1.0)); 
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1)*(-1.0));
          }   
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d3-d1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d3-d1)*(-1.0));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d3-d2)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d3-d2)*(-1.0));
          }
     }
     else if (mode==4)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1));
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1));
          
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d3-d1));
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d3-d1));
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d3-d2));
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d3-d2));
          }
     }
     else if (mode==5)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1)*(-1.0));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1)*(-1.0));
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d3-d1)*(-1.0));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d3-d1)*(-1.0));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d3-d2)*(-1.0));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d3-d2)*(-1.0));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }
     else if (mode==6)
     {
          if (face.bound_face_edge1[i]>=fst_row+1 && face.bound_face_edge1[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge1[i]-fst_row-1].r=b_x[face.bound_face_edge1[i]-fst_row-1].r+creal(H0_b1*(d2-d1));
               b_x[face.bound_face_edge1[i]-fst_row-1].i=b_x[face.bound_face_edge1[i]-fst_row-1].i+cimag(H0_b1*(d2-d1));         
               b_y[face.bound_face_edge1[i]-fst_row-1].r=b_y[face.bound_face_edge1[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge1[i]-fst_row-1].i=b_y[face.bound_face_edge1[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge2[i]>=fst_row+1 && face.bound_face_edge2[i]<=(m_loc+fst_row+1))
          {
               b_x[face.bound_face_edge2[i]-fst_row-1].r=b_x[face.bound_face_edge2[i]-fst_row-1].r+creal(H0_b2*(d3-d1));
               b_x[face.bound_face_edge2[i]-fst_row-1].i=b_x[face.bound_face_edge2[i]-fst_row-1].i+cimag(H0_b2*(d3-d1));
               b_y[face.bound_face_edge2[i]-fst_row-1].r=b_y[face.bound_face_edge2[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge2[i]-fst_row-1].i=b_y[face.bound_face_edge2[i]-fst_row-1].i+0.0;
          }
          if (face.bound_face_edge3[i]>=fst_row+1 && face.bound_face_edge3[i]<=(m_loc+fst_row+1))
          {  
               b_x[face.bound_face_edge3[i]-fst_row-1].r=b_x[face.bound_face_edge3[i]-fst_row-1].r+creal(H0_b3*(d3-d2));
               b_x[face.bound_face_edge3[i]-fst_row-1].i=b_x[face.bound_face_edge3[i]-fst_row-1].i+cimag(H0_b3*(d3-d2));
               b_y[face.bound_face_edge3[i]-fst_row-1].r=b_y[face.bound_face_edge3[i]-fst_row-1].r+0.0;
               b_y[face.bound_face_edge3[i]-fst_row-1].i=b_y[face.bound_face_edge3[i]-fst_row-1].i+0.0;
          }
     }

}

