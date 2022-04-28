#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "math.h"

int main()
{

  int m=0;
  printf("    Input m: ");
  scanf("%i",&m);

  matrix Coords = new_matrix(m*m, 2); // x,y
  matrix NodesNum = new_matrix(2 * (m - 1) * (m - 1), 3);

  vector boundary = new_vector(4*m-4);

  /// Coords from 1 ~ m*m
  int k=0;
  for (int i=1;i<=m;i++){
    for (int j=1;j<=m;j++){
      k++;
      mget(Coords, k, 1) = (double)(j-1)/(m-1);
      mget(Coords, k, 2) = (double)(i-1)/(m-1);
    }
  }

  //print_matrix(&Coords);

  /// Nodes number for every triangle

/*
  for (int i = 1; i <= m * m - (m + 2)+1; i++)
  {
    //printf("%d\n",i);
    if (i % m != 0)
    {
      true_num++;
      mget(NodesNum, true_num, 1) = 1 + i - 1;
      mget(NodesNum, true_num, 2) = 2 + i - 1;
      mget(NodesNum, true_num, 3) = m + 2 + i - 1;
    }
  }


  for (int i = 1; i <= m * m - (m + 2)+1; i++)
  {
    if (i % m != 0)
    {
      true_num++;
      mget(NodesNum, true_num, 1) = 1 + i - 1;
      mget(NodesNum, true_num, 2) = m + 2 + i - 1;
      mget(NodesNum, true_num, 3) = m + 1 + i - 1;
    }
  }
*/

  // merge into one
  int true_num = 1;
  for (int i = 1; i <= m * m - (m + 2)+1; i=i+1)
  {
    //printf("%d\n",i);
    if (i % m != 0)
    {
      mget(NodesNum, true_num, 1) = 1 + i - 1;
      mget(NodesNum, true_num, 2) = 2 + i - 1;
      mget(NodesNum, true_num, 3) = m + 2 + i - 1;
      mget(NodesNum, true_num+1, 1) = 1 + i - 1;
      mget(NodesNum, true_num+1, 2) = m + 2 + i - 1;
      mget(NodesNum, true_num+1, 3) = m + 1 + i - 1;
      true_num=true_num+2;
    }
  }


  //print_matrix(&NodesNum);


  /// Boundary
  k=0;

  for (int i=1;i<=m;i++){ //y-
    k++;
    vget(boundary, k) = i;
  }

  for (int i=m+1;i<=m*m;i=i+m){ //x-
    //printf("%d\n",i);
    k++;
    vget(boundary, k) = i;
  }

  for (int i=2*m;i<=m*m;i=i+m){ // x+
    k++;
    vget(boundary, k) = i;
  }

  for (int i=m*m-m+2;i<=m*m-1;i=i+1){ //y+
    k++;
    vget(boundary, k) = i;
  }
  //print_vector(&boundary);

  /// Stiffness matrix and force vector
  int Nnodes=m*m;
  int Nele=2 * (m - 1) * (m - 1);
  matrix K = new_matrix(Nnodes, Nnodes);
  vector F = new_vector(Nnodes);


  /// Loop over element to do the assembly
  vector tri_nodes_number = new_vector(3);
  matrix Me = new_matrix(3, 3);
  matrix Be = new_matrix(2, 3);
  matrix Ke = new_matrix(3, 3);
  vector Fe = new_vector(3);


  for (int e=1;e<=Nele;e++){  // Loop over element
    for (int i=1;i<=3;i++){
      vget(tri_nodes_number,i) = mget(NodesNum,e,i);
    }
    //print_vector(&tri_nodes_number);
    for (int i=1;i<=3;i++){
      int NodesOnTri = vget(tri_nodes_number,i);
      mget(Me, i, 1) = 1;
      mget(Me, i, 2) = mget(Coords,NodesOnTri,1); //x
      mget(Me, i, 3) = mget(Coords,NodesOnTri,2); //y
    }
    //print_matrix(&Me);

    double det_33(const matrix* A);
    double Area = det_33(&Me)/2;
    //printf("Area=%10.9e\n",Area);

    // Be = Element shape function derivative matrices
    mget(Be,1,1)=(mget(Coords,(int)vget(tri_nodes_number,2),2)-mget(Coords,(int)vget(tri_nodes_number,3),2))/2/Area;
    mget(Be,1,2)=(mget(Coords,(int)vget(tri_nodes_number,3),2)-mget(Coords,(int)vget(tri_nodes_number,1),2))/2/Area;
    mget(Be,1,3)=(mget(Coords,(int)vget(tri_nodes_number,1),2)-mget(Coords,(int)vget(tri_nodes_number,2),2))/2/Area;
    mget(Be,2,1)=(mget(Coords,(int)vget(tri_nodes_number,3),1)-mget(Coords,(int)vget(tri_nodes_number,2),1))/2/Area;
    mget(Be,2,2)=(mget(Coords,(int)vget(tri_nodes_number,1),1)-mget(Coords,(int)vget(tri_nodes_number,3),1))/2/Area;
    mget(Be,2,3)=(mget(Coords,(int)vget(tri_nodes_number,2),1)-mget(Coords,(int)vget(tri_nodes_number,1),1))/2/Area;
    //print_matrix(&Be);

    // build local stiffness matrix
    matrix matrix_mult_transpose(const matrix* A, const matrix* B);
    matrix temp_Ke = matrix_mult_transpose(&Be,&Be);
    matrix scale_matrix(const matrix* A, double scaling);
    Ke= scale_matrix(&temp_Ke,Area);
    //print_matrix(&Ke);

    // build local force vector
    // center point
    double xg =(mget(Coords,(int)vget(tri_nodes_number,1),1)+mget(Coords,(int)vget(tri_nodes_number,2),1)
    +mget(Coords,(int)vget(tri_nodes_number,3),1))/3;
    double yg =(mget(Coords,(int)vget(tri_nodes_number,1),2)+mget(Coords,(int)vget(tri_nodes_number,2),2)
    +mget(Coords,(int)vget(tri_nodes_number,3),2))/3;
    double eval_f = 2*M_PI*M_PI*sin(M_PI*xg)*sin(M_PI*yg);

    for (int i=1;i<=3;i++){
      vget(Fe,i)=Area*eval_f/3;
    }

    // assemby for stiffness
    for (int i=1;i<=3;i++){
      for (int j=1;j<=3;j++){
      //printf("i=%d \n",(int)vget(tri_nodes_number,i));
      //printf("j=%d \n",(int)vget(tri_nodes_number,j));
      mget(K,(int)vget(tri_nodes_number,i),(int)vget(tri_nodes_number,j))
      =mget(K,(int)vget(tri_nodes_number,i),(int)vget(tri_nodes_number,j))+mget(Ke,i,j);
      }
    }

    // assemby for force
    for (int i=1;i<=3;i++){
      vget(F,(int)vget(tri_nodes_number,i))
      =vget(F,(int)vget(tri_nodes_number,i))+vget(Fe,i);
    }
  }

  // strong enforce BC => set LHS = 1 and set RHS = 0
  for (int i=1;i<=4*m-4;i++){
    vget(F,(int)vget(boundary,i))=0;
    for (int j=1;j<m*m;j++){
      mget(K,(int)vget(boundary,i),j)=0;
      mget(K,j,(int)vget(boundary,i))=0;
    }
    mget(K,(int)vget(boundary,i),(int)vget(boundary,i))=1;
  }


  //print_matrix(&K);
  //print_vector(&F);

  /// solve Linear System
  //vector u = solve(&K,&F);

  //TO DO: try to save stiffness matrix as sparse matrix!!!!!!!!!!!!!!!
  vector solveCG(const matrix* A, const vector* b);
  vector u = solveCG(&K,&F);
  //print_vector(&u);

  /// postprosessing
/*
  matrix u_OnGrid = new_matrix(m, m);

  int k=0;
  for (int i=1;i<=m;i++){
    for (int j=1;j<=m;j++){
      k++;
      mget(u_OnGrid,i,j)=vget(u,k);
    }
  }
*/
    /// Print solution to file
    char filename[] = "output.data";

    // open file
    FILE* outfile = fopen(filename,"w");

    // output data
    for (int i=1;i<=m*m;i++){
      fprintf(outfile,"%25.20e  ",mget(Coords,i,1));
      fprintf(outfile,"%25.20e  ",mget(Coords,i,2));
      fprintf(outfile,"%25.20e\n",vget(u,i));
    }

    // close file
    fclose(outfile);

    // Call python script to plot
    system("python3.8 phase_plot.py");


    /// Output Tecplot
    char filenametec[] = "output.tec";
    FILE* outfiletec = fopen(filenametec,"w");
    fprintf(outfiletec,"%s \n","Titile = 'Poisson' ");
    fprintf(outfiletec,"%s \n","VARIABLES = X, Y, U, Error");
    fprintf(outfiletec,"zone N= %d, E = %d",Nnodes,Nele);
    fprintf(outfiletec,"\n");
    fprintf(outfiletec,"DATAPACKING=POINT, ZONETYPE=FETRIANGLE");
    fprintf(outfiletec,"\n");
    for (int i=1;i<=Nnodes;i++){
      double Exact =sin(M_PI*mget(Coords,i,1))*sin(M_PI*mget(Coords,i,2));
      fprintf(outfiletec,"%25.20e  %25.20e  %25.20e  %25.20e\n",mget(Coords,i,1),mget(Coords,i,2),vget(u,i),fabs(Exact-vget(u,i)));
    }
    for (int i=1;i<=Nele;i++){
      fprintf(outfiletec,"%d  %d  %d",(int)mget(NodesNum,i,1),(int)mget(NodesNum,i,2),(int)mget(NodesNum,i,3));
      if (i != Nele){
        fprintf(outfiletec,"\n");
      }
    }


}
