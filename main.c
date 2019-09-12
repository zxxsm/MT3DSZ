// Author:ZHUXIAOXIONG
// Institution:NUDT 
// MTSZ main.c 
#include "main.h" 

//-------------------------------------------------------------------------------------------
// MTSZ MAIN
//-------------------------------------------------------------------------------------------

int main ( int argc, char **argv )
{
 
        // ---------------Task_init------------------------//
        // MPI initialize;
        MTSZ_task_init(&argc,&argv);

	// ---------------Prepoces-------------------------//
        // Prepocess before compute;
   	Prepocess();

        // ---------------Elem_compute---------------------//
        // Compute elem;
        // Elem assemble;
        Elem_compute();

        // ---------------Main_Solver----------------------//
        // Matrix solver;
        // Post prepocess;
        Main_Solver();

        // ---------------Out_put_result-------------------//
        // Write result to file;

}
//-------------------------------------------------------------------------------------------
