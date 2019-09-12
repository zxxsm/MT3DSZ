// Author:ZHUXIAOXIONG
// Institution:NUDT
// MTSZ Post_prepocss.c
#include "post_prepocess.h"

//-------------------------------------------------------------------------------------------
// 1.post_prepocess();
//-------------------------------------------------------------------------------------------
void post_prepocess(int ifre)
{
     int i,j,k;
     int ind_elem,ind_loc;
     int ind_edge1,ind_edge2,ind_edge3,ind_edge4,ind_edge5,ind_edge6;
     int node1,node2,node3,node4;
     double x0,x1,x2,x3,x4;
     double y0,y1,y2,y3,y4;
     double z0,z1,z2,z3,z4;
     double a1,a2,a3,a4;
     double b1,b2,b3,b4;
     double c1,c2,c3,c4;
     double d1,d2,d3,d4;
     double N1,N2,N3,N4,N5,N6;
     double l1,l2,l3,l4,l5,l6;
     double Ve;

     struct doublecomplex_z *ExA;
     struct doublecomplex_z *ExB;
     struct doublecomplex_z *EyA;
     struct doublecomplex_z *EyB;
     struct doublecomplex_z *HxA;
     struct doublecomplex_z *HxB;
     struct doublecomplex_z *HyA;
     struct doublecomplex_z *HyB;
         
     struct doublecomplex_z  Zxy,Zyx;
     double ap_ris_xy,ap_ris_yx;
     double alph_xy  ,alph_yx;

     char  output_file[50];
     FILE *out_file;
     /*for (i=1;i<=mesh_par.edge_number_total;i++)
          printf("i=%d,br=%lf,bi=%lf\n",i,b_x_total[i-1].r,b_x_total[i-1].i);

     for (i=1;i<=mesh_par.edge_number_total;i++)
          printf("i=%d,b1r=%lf,b1i=%lf\n",i,b_y_total[i-1].r,b_y_total[i-1].i);

     exit(0); 
     */
     if (mynewrank==0)
     {
          //open Fre[ifre]file
          if (bound_mode==1)
          {
              sprintf(output_file,"D_Output_FRE=%lf",frequency[ifre]);
          }
          else if  (bound_mode==2)
          {
              sprintf(output_file,"N_Output_FRE=%lf",frequency[ifre]); 
          }    
	  //printf("%s\n",output_file);
	  out_file=fopen(output_file,"w+");

	  fprintf( out_file," frequency = %12lf\n",frequency[ifre]);
		  
	  fprintf( out_file," num,    x0_loc,    y0_loc,    z0_loc,    xy_ris,    yx_ris,   alph_xy,   alph_xy\n");
		  
          ExA=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
          ExB=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
          EyA=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
          EyB=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));

          HxA=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
          HxB=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
          HyA=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
          HyB=(struct doublecomplex_z*)malloc(col_elem_num*sizeof(struct doublecomplex_z));
 
          printf ("col_elem_num=%d,myrank=%d,mynewrank=%d\n",col_elem_num,myrank,mynewrank);

          for (i=0;i<=col_elem_num-1;i++)
          {
               ExA[i].r=0.0;
               ExA[i].i=0.0;
               EyA[i].r=0.0;
               EyA[i].i=0.0;
               ExB[i].r=0.0;
               ExB[i].i=0.0;
               EyB[i].r=0.0;
               EyB[i].i=0.0;

               HxA[i].r=0.0;
               HxA[i].i=0.0;
               HyA[i].r=0.0;
               HyA[i].i=0.0;
               HxB[i].r=0.0;
               HxB[i].i=0.0;
               HyB[i].r=0.0;
               HyB[i].i=0.0;
          }
 
          for (i=0;i<=col_elem_num-1;i++)
          {
               ind_elem=col_elem[i];
               ind_loc =col_num [i];

               printf("ind_elem=%d,ind_loc=%d\n",ind_elem,ind_loc);

               ind_edge1=elem_edge1[ind_elem];
               ind_edge2=elem_edge2[ind_elem];
               ind_edge3=elem_edge3[ind_elem];
               ind_edge4=elem_edge4[ind_elem];
               ind_edge5=elem_edge5[ind_elem];
               ind_edge6=elem_edge6[ind_elem];
        
               x0=collect_x[ind_loc];
               y0=collect_y[ind_loc];
               z0=collect_z[ind_loc];

               printf("x0=%lf,y0=%lf,z0=%lf\n",x0,y0,z0);

               node1=elem_node1[ind_elem];
               node2=elem_node2[ind_elem];
               node3=elem_node3[ind_elem];
               node4=elem_node4[ind_elem];

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

               Ve= (x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2\
                 - x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1\
                 + x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1\
                 - x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1)/6.0;

               a1 = x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2;
               a2 =-x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1;
               a3 = x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1;
               a4 =-x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1;

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
           
               N1=l1*((a1+b1*x0+c1*y0+d1*z0)*b2-(a2+b2*x0+c2*y0+d2*z0)*b1)/36/Ve/Ve;
               N2=l2*((a1+b1*x0+c1*y0+d1*z0)*b3-(a3+b3*x0+c3*y0+d3*z0)*b1)/36/Ve/Ve;
               N3=l3*((a1+b1*x0+c1*y0+d1*z0)*b4-(a4+b4*x0+c4*y0+d4*z0)*b1)/36/Ve/Ve;

               N4=l4*((a2+b2*x0+c2*y0+d2*z0)*b3-(a3+b3*x0+c3*y0+d3*z0)*b2)/36/Ve/Ve;
               N5=l5*((a4+b4*x0+c4*y0+d4*z0)*b2-(a2+b2*x0+c2*y0+d2*z0)*b4)/36/Ve/Ve;
               N6=l6*((a3+b3*x0+c3*y0+d3*z0)*b4-(a4+b4*x0+c4*y0+d4*z0)*b3)/36/Ve/Ve;
      
               if (node1>node2) N1=N1*(-1.0);
               if (node1>node3) N2=N2*(-1.0);
               if (node1>node4) N3=N3*(-1.0);
               if (node2>node3) N4=N4*(-1.0);
               if (node4>node2) N5=N5*(-1.0);
               if (node3>node4) N6=N6*(-1.0);

               ExA[i].r=N1*b_x_total[ind_edge1-1].r+N2*b_x_total[ind_edge2-1].r+N3*b_x_total[ind_edge3-1].r\
                      + N4*b_x_total[ind_edge4-1].r+N5*b_x_total[ind_edge5-1].r+N6*b_x_total[ind_edge6-1].r;
 
               ExA[i].i=N1*b_x_total[ind_edge1-1].i+N2*b_x_total[ind_edge2-1].i+N3*b_x_total[ind_edge3-1].i\
                      + N4*b_x_total[ind_edge4-1].i+N5*b_x_total[ind_edge5-1].i+N6*b_x_total[ind_edge6-1].i;

               ExB[i].r=N1*b_y_total[ind_edge1-1].r+N2*b_y_total[ind_edge2-1].r+N3*b_y_total[ind_edge3-1].r\
                      + N4*b_y_total[ind_edge4-1].r+N5*b_y_total[ind_edge5-1].r+N6*b_y_total[ind_edge6-1].r;
 
               ExB[i].i=N1*b_y_total[ind_edge1-1].i+N2*b_y_total[ind_edge2-1].i+N3*b_y_total[ind_edge3-1].i\
                      + N4*b_y_total[ind_edge4-1].i+N5*b_y_total[ind_edge5-1].i+N6*b_y_total[ind_edge6-1].i;

               printf("ExA=%lf+%lfi,ExB=%lf+%lfi\n",ExA[i].r,ExA[i].i,ExB[i].r,ExB[i].i);
            
               N1=l1*((a1+b1*x0+c1*y0+d1*z0)*c2-(a2+b2*x0+c2*y0+d2*z0)*c1)/36/Ve/Ve;
               N2=l2*((a1+b1*x0+c1*y0+d1*z0)*c3-(a3+b3*x0+c3*y0+d3*z0)*c1)/36/Ve/Ve;
               N3=l3*((a1+b1*x0+c1*y0+d1*z0)*c4-(a4+b4*x0+c4*y0+d4*z0)*c1)/36/Ve/Ve;

               N4=l4*((a2+b2*x0+c2*y0+d2*z0)*c3-(a3+b3*x0+c3*y0+d3*z0)*c2)/36/Ve/Ve;
               N5=l5*((a4+b4*x0+c4*y0+d4*z0)*c2-(a2+b2*x0+c2*y0+d2*z0)*c4)/36/Ve/Ve;
               N6=l6*((a3+b3*x0+c3*y0+d3*z0)*c4-(a4+b4*x0+c4*y0+d4*z0)*c3)/36/Ve/Ve;
  
               if (node1>node2) N1=N1*(-1.0);
               if (node1>node3) N2=N2*(-1.0);
               if (node1>node4) N3=N3*(-1.0);
               if (node2>node3) N4=N4*(-1.0);
               if (node4>node2) N5=N5*(-1.0);
               if (node3>node4) N6=N6*(-1.0);

               EyA[i].r=N1*b_x_total[ind_edge1-1].r+N2*b_x_total[ind_edge2-1].r+N3*b_x_total[ind_edge3-1].r\
                      + N4*b_x_total[ind_edge4-1].r+N5*b_x_total[ind_edge5-1].r+N6*b_x_total[ind_edge6-1].r;
 
               EyA[i].i=N1*b_x_total[ind_edge1-1].i+N2*b_x_total[ind_edge2-1].i+N3*b_x_total[ind_edge3-1].i\
                      + N4*b_x_total[ind_edge4-1].i+N5*b_x_total[ind_edge5-1].i+N6*b_x_total[ind_edge6-1].i;

               EyB[i].r=N1*b_y_total[ind_edge1-1].r+N2*b_y_total[ind_edge2-1].r+N3*b_y_total[ind_edge3-1].r\
                      + N4*b_y_total[ind_edge4-1].r+N5*b_y_total[ind_edge5-1].r+N6*b_y_total[ind_edge6-1].r;
 
               EyB[i].i=N1*b_y_total[ind_edge1-1].i+N2*b_y_total[ind_edge2-1].i+N3*b_y_total[ind_edge3-1].i\
                      + N4*b_y_total[ind_edge4-1].i+N5*b_y_total[ind_edge5-1].i+N6*b_y_total[ind_edge6-1].i;

               printf("EyA=%lf+%lfi,EyB=%lf+%lfi\n",EyA[i].r,EyA[i].i,EyB[i].r,EyB[i].i);

               N1=2*l1*(c1*d2-d1*c2)/36/Ve/Ve;
               N2=2*l2*(c1*d3-d1*c3)/36/Ve/Ve;
               N3=2*l3*(c1*d4-d1*c4)/36/Ve/Ve;

               N4=2*l4*(c2*d3-d2*c3)/36/Ve/Ve;
               N5=2*l5*(c4*d2-d4*c2)/36/Ve/Ve;
               N6=2*l6*(c3*d4-d3*c4)/36/Ve/Ve;

               if (node1>node2) N1=N1*(-1.0);
               if (node1>node3) N2=N2*(-1.0);
               if (node1>node4) N3=N3*(-1.0);
               if (node2>node3) N4=N4*(-1.0);
               if (node4>node2) N5=N5*(-1.0);
               if (node3>node4) N6=N6*(-1.0);

               HxA[i].r=N1*b_x_total[ind_edge1-1].r+N2*b_x_total[ind_edge2-1].r+N3*b_x_total[ind_edge3-1].r\
                      + N4*b_x_total[ind_edge4-1].r+N5*b_x_total[ind_edge5-1].r+N6*b_x_total[ind_edge6-1].r;
 
               HxA[i].i=N1*b_x_total[ind_edge1-1].i+N2*b_x_total[ind_edge2-1].i+N3*b_x_total[ind_edge3-1].i\
                      + N4*b_x_total[ind_edge4-1].i+N5*b_x_total[ind_edge5-1].i+N6*b_x_total[ind_edge6-1].i;

               HxB[i].r=N1*b_y_total[ind_edge1-1].r+N2*b_y_total[ind_edge2-1].r+N3*b_y_total[ind_edge3-1].r\
                      + N4*b_y_total[ind_edge4-1].r+N5*b_y_total[ind_edge5-1].r+N6*b_y_total[ind_edge6-1].r;
 
               HxB[i].i=N1*b_y_total[ind_edge1-1].i+N2*b_y_total[ind_edge2-1].i+N3*b_y_total[ind_edge3-1].i\
                      + N4*b_y_total[ind_edge4-1].i+N5*b_y_total[ind_edge5-1].i+N6*b_y_total[ind_edge6-1].i;

               H_compute(&HxA[i].r,&HxA[i].i); 
               H_compute(&HxB[i].r,&HxB[i].i);  

               printf("HxA=%lf+%lfi,HxB=%lf+%lfi\n",HxA[i].r,HxA[i].i,HxB[i].r,HyB[i].i);

               N1=2*l1*(d1*b2-b1*d2)/36/Ve/Ve;
               N2=2*l2*(d1*b3-b1*d3)/36/Ve/Ve;
               N3=2*l3*(d1*b4-b1*d4)/36/Ve/Ve;

               N4=2*l4*(d2*b3-b2*d3)/36/Ve/Ve;
               N5=2*l5*(d4*b2-b4*d2)/36/Ve/Ve;
               N6=2*l6*(d3*b4-b3*d4)/36/Ve/Ve;

               if (node1>node2) N1=N1*(-1.0);
               if (node1>node3) N2=N2*(-1.0);
               if (node1>node4) N3=N3*(-1.0);
               if (node2>node3) N4=N4*(-1.0);
               if (node4>node2) N5=N5*(-1.0);
               if (node3>node4) N6=N6*(-1.0);

               HyA[i].r=N1*b_x_total[ind_edge1-1].r+N2*b_x_total[ind_edge2-1].r+N3*b_x_total[ind_edge3-1].r\
                      + N4*b_x_total[ind_edge4-1].r+N5*b_x_total[ind_edge5-1].r+N6*b_x_total[ind_edge6-1].r;
 
               HyA[i].i=N1*b_x_total[ind_edge1-1].i+N2*b_x_total[ind_edge2-1].i+N3*b_x_total[ind_edge3-1].i\
                      + N4*b_x_total[ind_edge4-1].i+N5*b_x_total[ind_edge5-1].i+N6*b_x_total[ind_edge6-1].i;

               HyB[i].r=N1*b_y_total[ind_edge1-1].r+N2*b_y_total[ind_edge2-1].r+N3*b_y_total[ind_edge3-1].r\
                      + N4*b_y_total[ind_edge4-1].r+N5*b_y_total[ind_edge5-1].r+N6*b_y_total[ind_edge6-1].r;
 
               HyB[i].i=N1*b_y_total[ind_edge1-1].i+N2*b_y_total[ind_edge2-1].i+N3*b_y_total[ind_edge3-1].i\
                      + N4*b_y_total[ind_edge4-1].i+N5*b_y_total[ind_edge5-1].i+N6*b_y_total[ind_edge6-1].i;

               H_compute(&HyA[i].r,&HyA[i].i);  
               H_compute(&HyB[i].r,&HyB[i].i);  
               printf("HyA=%lf+%lfi,HyB=%lf+%lfi\n",HyA[i].r,HyA[i].i,HyB[i].r,HyB[i].i);

               Zxy=Zxy_compute(ExA[i],EyA[i],HxA[i],HyA[i],ExB[i],EyB[i],HxB[i],HyB[i]);
               Zyx=Zyx_compute(ExA[i],EyA[i],HxA[i],HyA[i],ExB[i],EyB[i],HxB[i],HyB[i]);
            
               ap_ris_xy=(Zxy.i*Zxy.i+Zxy.r*Zxy.r)/omega/mu0;
               ap_ris_yx=(Zyx.i*Zyx.i+Zyx.r*Zyx.r)/omega/mu0;
               alph_xy=-atan(Zxy.i/Zxy.r)*180/pi;
               alph_yx=-atan(Zyx.i/Zyx.r)*180/pi;
               
               printf("xy_ris=%lf,xy_aph=%lf\n",ap_ris_xy,alph_xy);
               printf("yx_ris=%lf,yx_aph=%lf\n",ap_ris_yx,alph_yx);
			   
			   
               fprintf(out_file,"%4d,  %8lf,  %8lf,  %8lf,  %8lf,  %8lf,  %8lf,  %8lf\n",\
               i+1,x0,y0,z0,ap_ris_xy,ap_ris_yx,alph_xy,alph_yx);
			   
	  } 
	  
	 fclose(out_file);
     }

}

//-------------------------------------------------------------------------------------------
// H_compute
//-------------------------------------------------------------------------------------------
void H_compute(double *real, double *imag)
{

     complex double temp;
     
     temp=*real + *imag* _Complex_I;
     temp=temp/((-1.0)*omega *mu0 * _Complex_I);
 
     *real = creal(temp);
     *imag = cimag(temp);

}
//-------------------------------------------------------------------------------------------
// Zxy_compute
//-------------------------------------------------------------------------------------------
struct doublecomplex_z Zxy_compute(struct doublecomplex_z ExA, struct doublecomplex_z EyA,\
                           struct doublecomplex_z HxA, struct doublecomplex_z HyA,\
                           struct doublecomplex_z ExB, struct doublecomplex_z EyB,\
                           struct doublecomplex_z HxB, struct doublecomplex_z HyB)
{
     complex double temp;               
     struct doublecomplex_z ctemp;
  
     temp=((ExB.r+ExB.i* _Complex_I)*(HxA.r+HxA.i* _Complex_I)-  \
           (ExA.r+ExA.i* _Complex_I)*(HxB.r+HxB.i* _Complex_I))/ \
          ((HxA.r+HxA.i* _Complex_I)*(HyB.r+HyB.i* _Complex_I)-  \
           (HxB.r+HxB.i* _Complex_I)*(HyA.r+HyA.i* _Complex_I));    

     ctemp.r=creal(temp);
     ctemp.i=cimag(temp);
     return ctemp;
}


struct doublecomplex_z Zyx_compute(struct doublecomplex_z ExA, struct doublecomplex_z EyA,\
                           struct doublecomplex_z HxA, struct doublecomplex_z HyA,\
                           struct doublecomplex_z ExB, struct doublecomplex_z EyB,\
                           struct doublecomplex_z HxB, struct doublecomplex_z HyB)
{
     complex double temp;               
     struct doublecomplex_z ctemp;
 
     temp=((EyA.r+EyA.i* _Complex_I)*(HyB.r+HyB.i* _Complex_I)-  \
           (EyB.r+EyB.i* _Complex_I)*(HyA.r+HyA.i* _Complex_I))/ \
          ((HxA.r+HxA.i* _Complex_I)*(HyB.r+HyB.i* _Complex_I)-  \
           (HxB.r+HxB.i* _Complex_I)*(HyA.r+HyA.i* _Complex_I));    

     ctemp.r=creal(temp);
     ctemp.i=cimag(temp);

     return ctemp;
}

