#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <time.h>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <bits/stdc++.h>

const short int N = 200; //Tamaño del arreglo
const short int T = 200; // Tiempo máximo
const short int r = 40; //Radio en celdas de la región caliente
const double E_A = pow(10,6); // Energía de activación, como multiplo de la constante de Boltzmann
const double R = 8.31446261; // Constante del Gas ideal

/*------------------Clase Paper1-------------------*/

class Material1
{
private:
  int h_old[N][N];
  int h_new[N][N];
  int areas[N*N+1000];

public:
  void fill (void);
  void fill_circle(void);
  void evolution(int t);
  void evolution_aux1 (int a, int b, int c, int d, int i, int j);
  void evolution_aux2 (int a, int b, int c, int d, int i, int j, int prob);
  void array_change(void);
  void print_array(const char * Arreglo);
  int h_old_get(int i, int j);
};
  
void Material1::fill (void){
  int count = 1;
  std::mt19937 gen(N);
  std::bernoulli_distribution rand(0.75);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      h_old[i][j] = 1;
    }
  }
  for(int u = 0; u<N; u++){
    h_old[0][u] = count;
    count++;
  }
  
  for(int i=1; i<N; i++){
    for(int j=0; j<N; j++){
      if (i == (N-1)){
	h_old[i][j] = h_old[0][j];
      }
      else if (j == (N-1)){
	h_old[i][j] = h_old[i][0];
      }
      else {
	if(rand(gen) == true){
	h_old[i][j]= count+1000;
	}
	else {
	  h_old[i][j]= count;
	}
	count++;
      }
    }
  }
}

void Material1::evolution_aux1 (int a, int b, int c, int d, int i, int j){
 
  // a=(i,j-1) ; b = (i-1,j) ; c = (i,j+1) ; d = (i+1,j)
  if (a==b && b==c && c==d){
    h_old[i][j] = a;
    h_new[i][j]=h_old[i][j];
  }

  else{
    if((a == b) && (b == c)){//Prueba si a,b,c son iguales
      if(i==0){
	h_new[i][j] = a;
	h_new[N-1][j] = a;
    }
      else if(j==0){
	h_new[i][j] = a;
	h_new[i][N-1] = a;
      }
      else{
	h_new[i][j] = a;
      }
      
    }
    
    else  if((b == c) && (c == d)){//Prueba si b,c,d son iguales
      if(i==0){
	h_new[i][j] = b;
	h_new[N-1][j] = b;
      }
      else if(j==0){
	h_new[i][j] = b;
	h_new[i][N-1] = b;
    }
      else{
	h_new[i][j] = b;
      }
      
    }
    
    else if((c == d) && (d == a)){//Prueba si c,d,a son iguales
      if(i==0){
	h_new[i][j] = c;
	h_new[N-1][j] = c;
      }
      else if(j==0){
	h_new[i][j] = c;
	h_new[i][N-1] = c;
      }
      else{
	h_new[i][j] = c;
      }
      
    }
    else if((d == a) && (a == b)){//Prueba si d,a,b son iguales
      if(i==0){
	h_new[i][j] = d;
	h_new[N-1][j] = d;
      }
      else if(j==0){
	h_new[i][j] = d;
	h_new[i][N-1] = d;
      }
      else{
	h_new[i][j] = d;
      }
      
    }
  }
}

void Material1::evolution_aux2 (int a, int b, int c, int d, int i, int j, int prob){
 if(1<=prob && prob<=25){
	h_new[i][j] = a;
      }
      else if (25<prob && prob<=50){
	h_new[i][j] = b;
      }
      else if (50<prob && prob<=75){
	h_new[i][j] = c;
      }
      else if (75<prob && prob<=100) {
	h_new[i][j] = d;
      }
      else{
	h_new[i][j] = h_old[i][j];
      }
}
void Material1::evolution(int t){
  int a,b,c,d; //variables auxiliares: a=(i,j-1) ; b = (i-1,j) ; c = (i,j+1) ; d = (i+1,j)
 std::mt19937 gen(t);
  std::uniform_int_distribution<> rand(1, 100);
  for(int i=0; i<(N-1); i++){
    for(int j=0; j<(N-1); j++){
      if (i == 0){
	if (j== 0){
	  a = h_old[i][N-2];
	}
	else{
	  a = h_old[i][j-1];
	}
	b = h_old[N-2][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
	evolution_aux2 (a, b, c,  d,  i, j, rand(gen));
      }
      else if (j == 0){
	a = h_old[i][N-2];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
	evolution_aux2 (a, b, c,  d,  i, j, rand(gen));
      }
      else {
	a = h_old[i][j-1];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
	evolution_aux2 (a, b, c,  d,  i, j, rand(gen));
      }
    }
  }

  for(int i=0; i<(N-1); i++){
    for(int j=0; j<(N-1); j++){
      if (i == 0){
	if (j== 0){
	  a = h_old[i][N-2];
	}
	else{
	  a = h_old[i][j-1];
	}
	b = h_old[N-2][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
      }
      else if (j == 0){
	a = h_old[i][N-2];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
      }
      else {
	a = h_old[i][j-1];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
      }
    }
  }
  

  
}
void Material1::array_change(void){
 for(int i=0; i<(N-1); i++){
    for(int j=0; j<(N-1); j++){
      if (i==0){
	if(j==0){
	  h_old[i][j] = h_new[i][j];
	h_old[N-1][j] = h_new[i][j];
	h_old[i][N-1] = h_new[i][j];
	}
	else{
	  h_old[i][j] = h_new[i][j];
	  h_old[N-1][j] = h_new[i][j];
	}
      }
      else if(j==0 && i!=0){
	h_old[i][j] = h_new[i][j];
	h_old[i][N-1] = h_new[i][j];
      }
      else{
	h_old[i][j] = h_new[i][j];
      }
    }
 }
 h_old[N-1][N-1] = h_old[0][0];
}

void Material1::print_array(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<h_old[i][j]<<"  ";
    }
    MiArchivo<<std::endl;
  }
  MiArchivo.close();
}

int Material1::h_old_get(int i, int j){
  return h_old[i][j];
}
  

/*------------------Clase Paper2-------------------*/


class Material
{
private:
  short int h_old[N][N];
  short int h_new[N][N];
  double Delta_E[N][N];
  short int neighborhood [25]; //Numero de segundos vecinos de Moore
  int areas [N*N];
  double tempt [N][N];
  bool frontier [N][N];
public:
  void fill_tempt (double T_min, double T_max);
  void fill_circle(void);
  double probability (int i, int j, double M_max, double D_Emin);
  void neighbor_def (int i, int j);
  int neighbor_get (int ii);
  bool is_boundary (int i, int j);
  long int Boundary_energy (void);
  void Delta_E_min (int i, int j, int t);
  void evolve (int t);
  void h_old_set(int i, int j, int v);
  void size_count (void);
  void print_array(const char * Arreglo);
  void print_gradient(const char * Arreglo);
  void metrics(double &mean_area, double &mean_size, bool dist);
  int mean_area (void);
};

void Material::fill_circle(void){
    
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      h_old[i][j] = 1;
    }
  }
  for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
	if(((i-(N/2))*(i-(N/2)))+((j-(N/2))*(j-(N/2)))<r*r){
	  h_old[i][j] = 2;
	}
      }
  }
}


void Material::fill_tempt(double T_min, double T_max){
  double a = r*(T_max-T_min);
  double b = T_min;
  double x = 0.0;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      tempt [i][j] = 500;
    }
  }
  
  for(int i=1; i<(N-1); i++){
    for(int j=1; j<(N-1); j++){
      if((((i-(N/2))*(i-(N/2)))+((j-(N/2))*(j-(N/2))))<(r*r)){
	tempt [i][j] = T_max ;
      }
       else{
	 x = sqrt(((1.0*i)-(1.0*N)/2)*((1.0*i)-(1.0*N)/2)+((1.0*j)-(1.0*N)/2)*((1.0*j)-(1.0*N)/2));
	 tempt [i][j] = a/(x)+b;
       }
    }
  }
}

double Material::probability (int i, int j, double M_max, double D_Emin){
  double M_ij = exp ((-E_A)/tempt [i][j]);
  if (Delta_E[i][j] < 0){
  return ((M_ij/M_max)*(Delta_E[i][j]/D_Emin));
  }
  else {
    return 0;
  }
}
void Material::neighbor_def (int i, int j){
  int n_index=0; //indice que reccore los vecinos
  for (n_index = 0; n_index<25; n_index++){
    neighborhood[n_index] = 0;
  }
  // Frontera interior
  if((i == 1 || i == (N-2)) && (j>=1 && j<=(N-2))){
    int i_aux = 0;
    int j_aux = 0;
    for (int ii = i-1; ii <= i+1; ii++){
      for(int jj = j-1; jj <= j+1; jj++){
	neighborhood[6+(3*i_aux)+j_aux] = h_old[ii][jj];
	j_aux++;
      }
      j_aux--;
      i_aux++;
    }
    if (i == 1){
      if(j==1){
	neighborhood[0] = h_old[N-2][N-2];
	for (int u = 1; u <= 4; u++){
	  neighborhood [u*5] = h_old[u-1][N-2];
	  neighborhood [(u*5)+4] = h_old [u-1][3];
	  neighborhood [20+u] = h_old[3][u-1];
	  neighborhood [u] = h_old [N-2] [u-1];
	}
      }
      else if (j==(N-2)){
	neighborhood[4] = h_old[1][1];
	for (int u = 1; u <= 4; u++){
	  neighborhood[u*5] = h_old[u-1][N-4];
	  neighborhood [(u*5)+4] = h_old [u-1][1];
	  neighborhood [19+u] = h_old[3][(N-5)+u];
	  neighborhood [u-1] = h_old[N-2][(N-5)+u];
	}
      }
      else{
	for (int u = 0; u < 5; u++){
	  neighborhood [u] = h_old[N-2][(j-2)+u];
	  neighborhood [5+u] = h_old[i-1][(j-2)+u];
	  neighborhood [10+u] = h_old [i][(j-2)+u];
	  neighborhood [15+u] = h_old [i+1][(j-2)+u];
	  neighborhood [20+u] = h_old [i+2][(j-2)+u];
	}
      }
    }
    else if(i == (N-2)){
      if (j==1){
	neighborhood[20] = h_old [1][N-2];
	for (int u = 1; u <= 4; u++){
	  neighborhood [((u-1)*5)] = h_old [(N-5)+u][N-2];
	  neighborhood [((u-1)*5)+4] = h_old [(N-5)+u][3];
	  neighborhood [u] = h_old [N-4][u-1];
	  neighborhood [20+u] = h_old [1][u-1];
	}
      }
      else if (j==(N-2)){
	neighborhood[24]=h_old[1][1];
	for(int u = 1; u<=4; u++){
	  neighborhood [u-1] = h_old[N-4][(N-5)+u];
	  neighborhood [(u-1)*5] = h_old [(N-5)+u][N-4];
	  neighborhood [((u-1)*5)+4] = h_old [(N-5)+u][1];
	  neighborhood [19+u] = h_old [1][(N-5)+u];
	}
      }
    } 
  }
  else if ((j==1||j==(N-2)) && (i > 1 && i < (N-2))){
    if (j==1){
      for (int u = 0; u<5 ; u++){
	neighborhood [(u*5)] = h_old[(i-2)+u][N-2];
	neighborhood [(u*5)+1] = h_old [(i-2)+u][j-1];
	neighborhood [(u*5)+2] = h_old [(i-2)+u][j];
	neighborhood [(u*5)+3] = h_old [(i-2)+u][j+1];
	neighborhood [(u*5)+4] = h_old [(i-2)+u][j+2];
      }
    }
    else if (j==(N-2)){
      for (int u = 0; u<5; u++){
	neighborhood [(u*5)] = h_old [(i-2)+u][j-2];
	neighborhood [(u*5)+1] = h_old [(i-2)+u][j-1];
	neighborhood [(u*5)+2] = h_old [(i-2)+u][j];
	neighborhood [(u*5)+3] = h_old [(i-2)+u][j+1];
	neighborhood [(u*5)+4] = h_old[(i-2)+u][1];
      }
    }
  }
  //Frontera exterior
  else if (i == 0 || i == (N-1)){ 
    if (i==0){//terminado
      if(j==0){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j; jj <= j+2; jj++){
	    neighborhood[12+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u < 3; u++){
	  neighborhood [7+u] = h_old [N-2][u];
	  neighborhood [2+u] = h_old [N-3][u];
	  neighborhood [(5*u)+11] = h_old [u][N-2];
	  neighborhood [(5*u)+10] = h_old [u][N-3];
	}
	neighborhood [0] = h_old [N-3][N-3];
	neighborhood [1] = h_old [N-3][N-2];
	neighborhood [5] = h_old [N-2][N-3];
	neighborhood [6] = h_old [N-2][N-2];
      }
      else if (j==1){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-1; jj <= j+2; jj++){
	    neighborhood[11+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for(int u = 1; u < 5; u++){
	  neighborhood [u] = h_old [N-3][u-1];
	  neighborhood [5+u] = h_old [N-2][u-1];
	  neighborhood [(5*u)+5] = h_old [u-1][N-2];
	}
	neighborhood [0] = h_old [N-3][N-2];
	neighborhood [5] = h_old [N-2][N-2];
      }
      
      else if (j==(N-1)){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-2; jj <= j; jj++){
	    neighborhood[10+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [5+u] = h_old [N-2][(N-3)+u];
	  neighborhood [u] = h_old [N-3][(N-3)+u];
	  neighborhood [(5*u)+13] = h_old [u][1];
	  neighborhood [(5*u)+14] = h_old [u][2];
	}	
	neighborhood [3] = h_old [N-3][1];
	neighborhood [4] = h_old [N-3][2];
	neighborhood [8] = h_old [N-2][1];
	neighborhood [9] = h_old [N-2][2];
      }
      else if (j == (N-2)){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-2; jj <= j+1; jj++){
	    neighborhood[10+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for (int u = 0 ; u<4 ; u++){
	  neighborhood [u] = h_old [N-3][(N-4)+u];
	  neighborhood [5+u] = h_old [N-2][(N-4)+u];
	  neighborhood [(5*u)+14] = h_old [u][1];
	}
	neighborhood [4] = h_old [N-3][1];
	neighborhood [9] = h_old [N-2][1];
      }
      
      else {// terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[10+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
        i_aux = 0;
	j_aux = 0;
	for (int ii = 0; ii < 2; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[N-3+ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
      }
    }
    
    else if (i == (N-1)){
      if (j==0){ // terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j; jj <= j+2; jj++){
	    neighborhood[2+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [5*u] = h_old [(N-3)+u][N-3];
	  neighborhood [(5*u)+1] = h_old [(N-3)+u][N-2];
	  neighborhood [17+u] = h_old [1][u];
	  neighborhood [22+u] = h_old [2][u];
	}
	neighborhood [15] = h_old [1][N-3];
	neighborhood [16] = h_old [1][N-2];
	neighborhood [20] = h_old [2][N-3];
	neighborhood [21] = h_old [2][N-2];
      }
      else if (j==1){//terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-1; jj <= j+2; jj++){
	    neighborhood[1+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for (int u =0; u < 4; u++){
	  neighborhood [16+u] = h_old [1][u];
	  neighborhood [21+u] = h_old [2][u];
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [5*u] = h_old [(N-3)+u][N-2];
	}
	neighborhood [15] = h_old [1][N-2];
	neighborhood [20] = h_old [2][N-2];
      }
      else if (j==N-1){//terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-2; jj <= j; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u < 3; u++){
	  neighborhood [15+u] = h_old [1][(N-3)+u];
	  neighborhood [20+u] = h_old [2][(N-3)+u];
	  neighborhood [(5*u)+3] = h_old [(N-3)+u][1];
	  neighborhood [(5*u)+4] = h_old [(N-3)+u][2];
	}
	neighborhood [18] = h_old [1][1];
	neighborhood [19] = h_old [1][2];
	neighborhood [23] = h_old [2][1];
	neighborhood [24] = h_old [2][2];
      }
      else if (j == N-2){//terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-2; jj <= j+1; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for (int u =0; u < 4; u++){
	  neighborhood [15+u] = h_old [1][(N-4)+u];
	  neighborhood [20+u] = h_old [2][(N-4)+u];
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [(5*u)+4] = h_old [(N-3)+u][1];
	}
	neighborhood [19] = h_old [1][1];
	neighborhood [24] = h_old [2][1];
      }
      else { // terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
        i_aux = 0;
	j_aux = 0;
	for (int ii = 0; ii < 2; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[15+(3*i_aux)+j_aux] = h_old[1+ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
      }
    }
  }
   else  if ((i == 1 || i == N-2) && (j == 0 || j == N-1)){
     if ( i == 1 && j==0){ //terminado
      int i_aux = 0;
      int j_aux = 0;
      for (int ii = i-1; ii <= i+2; ii++){
	for(int jj = j; jj <= j+2; jj++){
	  neighborhood[7+(4*i_aux)+j_aux] = h_old[ii][jj];
	  j_aux++;
	}
	j_aux-=2;
	i_aux++;
      }
      for (int u = 0; u<4 ; u++){
	neighborhood [(5*u)+5] = h_old [u][N-3];
	neighborhood [(5*u)+6] = h_old [u][N-2];
      }
      for(int u = 0; u<3 ; u++){
	neighborhood [2+u] = h_old [N-2][u];
      }
      neighborhood [0] = h_old [N-2][N-3];
      neighborhood [1] = h_old [N-2][N-2];
     }
     else if (i==1 && j==N-1){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-1; ii <= i+2; ii++){
	 for(int jj = j-2; jj <= j; jj++){
	   neighborhood[5+(4*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=2;
	i_aux++;
       }
       for (int u = 0; u<4 ; u++){
	 neighborhood [(5*u)+8] = h_old [u][1];
	 neighborhood [(5*u)+9] = h_old [u][2];
       }
       for(int u = 0; u<3 ; u++){
	 neighborhood [u] = h_old [N-2][(N-3)+u];
       }
       neighborhood [3] = h_old [N-2][1];
       neighborhood [4] = h_old [N-2][2];
     }
     else if (i==(N-2) && j == 0){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-2; ii <= i+1; ii++){
	 for(int jj = j; jj <= j+2; jj++){
	   neighborhood[2+(4*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
       for (int u = 0; u<4 ; u++){
	 neighborhood [(5*u)] = h_old [(N-4)+u][N-3];
	 neighborhood [(5*u)+1] = h_old [(N-4)+u][N-2];
       }
        for(int u = 0; u<3 ; u++){
	 neighborhood [22+u] = h_old [1][u];
       }
       neighborhood [20] = h_old [1][N-3];
       neighborhood [21] = h_old [1][N-2];
     }
     else if (i == (N-2) && j== (N-1)){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-2; ii <= i+1; ii++){
	 for(int jj = j-2; jj <= j; jj++){
	   neighborhood[(4*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
       for (int u = 0; u<4 ; u++){
	 neighborhood [(5*u)+3] = h_old [(N-4)+u][1];
	 neighborhood [(5*u)+4] = h_old [(N-4)+u][2];
       }
       for(int u = 0; u<3 ; u++){
	 neighborhood [20+u] = h_old [1][(N-3)+u];
       }
       neighborhood [23] = h_old [1][1];
       neighborhood [24] = h_old [1][2];
     }
   }
   else if ((j==0 || j==(N-1)) && (i>1 && i <N-2)){
     if (j == 0) {//terminado
       int i_aux = 0;
       int j_aux = 0;
        for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = j; jj <= j+2; jj++){
	   neighborhood[2+(5*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=3;
	 i_aux++;
       }
       i_aux = 0;
       j_aux = 0;
       for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = 0; jj < 2; jj++){
	   neighborhood[(5*i_aux)+j_aux] = h_old[ii][N-3+jj];
	    j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
     }
     else if (j == (N-1)){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = j-2; jj <= j; jj++){
	   neighborhood[(5*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=3;
	 i_aux++;
       }
       i_aux = 0;
       j_aux = 0;
       for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = 0; jj < 2; jj++){
	   neighborhood[3+(5*i_aux)+j_aux] = h_old[ii][1+jj];
	    j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
     }
   }
   
   else{
     n_index=0;
     for(int ii = i-2; ii <= i+2; ii++){
       for(int jj = j-2; jj <= j+2; jj++){
	 neighborhood[n_index] = h_old[ii][jj];
	 n_index++;
       }
     }
   }
}

bool Material::is_boundary (int i, int j){
  int actual = h_old [i][j];
  int test = 0;
  int i_aux =  0;
  int j_aux = 0;
  neighbor_def (i,j);
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj<3; jj++){
      test = neighborhood [6+(3*i_aux)+j_aux];
      if (test != actual){
	frontier [i][j] = true;
	return true;
      }
      j_aux++;
    }
    j_aux--;
    i_aux++; 
  }
  frontier [i][j] = false;
  return false;
}

long int Material::Boundary_energy (void){
   long int Energy = 0;
#pragma omp parallel for collapse (2) reduction (+:Energy)
 for (int ii= 0 ;ii<N;ii++){
   for(int jj=0;jj<N;jj++){
     for(int u=0;u<25;u++){
       if (h_old[ii][jj]==neighborhood[u]){
	 continue;
       }
       else {
	 Energy++;
       } 
     }
   }
 }
  return Energy/2;
  }

void Material::Delta_E_min(int i, int j,int t){
  long int EB_0 = 0;
  long int EB_f = 0;
  int EB = 0;
  int E_aux = 0;
  int value_aux = h_old[i][j];
  int count = 0;
  int rand_part = 0;
  double Delta_ET = 0.0;
  std::vector <int> S(neighborhood,neighborhood+25);
  

  std::mt19937 gen(t);
  std::uniform_real_distribution<> dis(0, 1);
  
  Delta_ET = -R*tempt [i][j]*log (dis(gen));
  EB_0 = Boundary_energy ();
  
  sort(S.begin(),S.end());
  S.erase(std::remove(S.begin(),S.end(), 0),S.end());
  S.shrink_to_fit();
  std::vector <int> Delta_EB(S.size());
  for (int u = 0; u<S.size(); u++){
    //std::cout<<u<<std::endl;
    if(u == 0){
      neighborhood [12] = S[0];
      h_old[i][j] = S[0];
      EB_f = Boundary_energy ();
      EB = EB_f-EB_0;
      if (EB<E_aux){
	E_aux=EB;
	h_new[i][j] = h_old[i][j];
      }
      Delta_EB[0] = EB;
    }
    else{
      if (S[u] == value_aux){
	Delta_EB[u] = 0;
	}
      else if(S[u-1]==S[u]){
	Delta_EB[u] = Delta_EB[u-1];
      }
      else{
	neighborhood [12] = S[u];
	h_old [i][j] = S[u];
	EB_f = Boundary_energy ();
	EB = EB_f-EB_0;
	if (EB<E_aux){
	  E_aux=EB;
	  h_new[i][j] = h_old[i][j];
	}
	Delta_EB[u]=EB;
      }
    }
  }
  
  for(int u = 0; u<S.size(); u++){
    if(Delta_EB[u] == E_aux){
      count++;
      }
  }
  if (count > 1){
    double part = 100.0/count;
    double rand_num = dis(gen) *100.0;
    rand_part = round(rand_num/part);
    
    if (rand_part==0){
      rand_part++;
    }
    count = 0;
    for(int u = 0; u<S.size(); u++){
      if(Delta_EB[u] == E_aux){
	count++;
	if(count == rand_part){
	  
	  h_new [i][j] = S [u];
	  break;
	}
      }  
    }
  }
  Delta_E [i][j] = EB-Delta_ET;
  h_old[i][j] = value_aux;
  }

void Material::evolve (int t){
  double D_Emin = 0.0;
  double M_max = exp ((-E_A)/tempt [N/2][N/2]);

  std::mt19937 gen(t);
  
 
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (is_boundary(i,j)==true){
	Delta_E_min(i,j, t);
      }
      else {
	Delta_E [i][j] = 0;
      }
    }
  }
  

  
  D_Emin = Delta_E [0][0];
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){  
      if(Delta_E [i][j] < D_Emin){
	D_Emin = Delta_E[i][j];
      }
    }
  }

 
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (Delta_E [i][j] != 0){
	std::bernoulli_distribution prob(probability (i, j, M_max, D_Emin));
	if (prob(gen) == false){
	  h_new[i][j] = h_old [i][j];
	}
      }
      else {
	h_new[i][j] = h_old[i][j];
      }
    }
  }
  
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (j == N-1){
	h_old [i][j] = h_new[i][0];
      }
      else if (i== N-1){
	h_old [i][j] = h_new [0][j];
      }
      else {
	h_old [i][j] = h_new [i][j];
      }
    }
  }
  
}

void Material::size_count (void){
  for(int u =0; u<(N*N+1000); u++){
    areas[u] = 0;
  }
  
  for (int i=0; i<N-1; i++){
    for(int j=0; j<N-1; j++){
      areas[(h_old [i][j])-1]++;
    }
  }
}

void Material::metrics(double &mean_area, double &mean_size, bool dist){
  
  double max_freq=0;;
  int size = N*N+1000;
  std::vector <int> areas_calc(areas,areas+(size));
  
  sort(areas_calc.begin(),areas_calc.end());
  areas_calc.erase(std::remove(areas_calc.begin(),areas_calc.end(), 0),areas_calc.end());
  areas_calc.shrink_to_fit();
  mean_area = 1.0* std::accumulate(areas_calc.begin(),areas_calc.end(),0);
  mean_area = 1.0*mean_area/areas_calc.size();
  
  std::map<int, double> counts;
  for (auto v : areas_calc)
    ++counts[v];
  
  
  for (auto v : areas_calc){
    if(max_freq<counts[v]){
      max_freq = counts[v];
    }
  }
  std::vector<double> radius(counts.size());
  std::vector<double> frecuency(counts.size());
  int u = 0;
  for (auto const &p : counts ){
    radius [u] = sqrt(p.first);
    frecuency[u] = p.second/max_freq;
    u++;
  }
  mean_size = 1.0* std::accumulate(radius.begin(),radius.end(),0);
  mean_size = mean_size/radius.size()*1.0;
  
  for (int v = 0; v<u; v++){
    radius[v] = radius[v]/mean_size;
  }
  if(dist == true){
    std::ofstream MiArchivo("distribucion_2.dat");
    for( int v= 0; v<u; v++){
      MiArchivo<<radius[v]<<'\t'<<frecuency[v]<<std::endl;
    }
     MiArchivo.close();
  }
}
int Material::neighbor_get (int ii){
  return neighborhood[ii];
}
void Material::h_old_set(int i, int j, int v){
  h_old[i][j] = v;
}

void Material::print_array(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<h_old [i][j]<<"  ";
    }
    MiArchivo<<std::endl;
  }
  MiArchivo.close();
}
void Material::print_gradient(const char * Arreglo)
{
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<=N/2;i++){
    MiArchivo<<i<<'\t'<<tempt [i][N/2]<<std::endl;;
  }
  MiArchivo.close();
}

/*------------------Main-------------------*/
int main (void){
  double mean_area = 0;
  double mean_size = 0;

  std::ofstream MiArchivo("size_evolv2.dat");
  
   Material1 initial_sys;
  initial_sys.fill();
  // initial_sys.print_array("Inicial1.dat");
  for(int t =0 ; t<103; t++){
    initial_sys.evolution(t);
    initial_sys.array_change();
  }
  // initial_sys.print_array("Inicial2.dat");
  Material granos;
  for (int i = 0; i<N ; i++){
    for (int j = 0; j<N; j++){
      granos.h_old_set(i,j,initial_sys.h_old_get(i,j));
    }
  }
  
  granos.fill_tempt(1400.0,1400.0);
  granos.size_count();
  granos.metrics(mean_area, mean_size, false);
  std::cout<<0<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
  MiArchivo<<0<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
  granos.print_array("C2_A1.dat");
  for (int t = 1; t <=40; t++){
    granos.evolve(t);

    if(t%10==0){
      granos.size_count();
      granos.metrics(mean_area, mean_size, false);
      std::cout<<t<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
      MiArchivo<<t<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
      if(t == 100){
	granos.print_array("C2_A2.dat");
      }
      else if (t == 200){
	granos.metrics(mean_area, mean_size, true);
	granos.print_array("C2_A3.dat");
      }
    }
  }
  granos.print_array("C2_A4.dat");
  
  return 0;
  
}
