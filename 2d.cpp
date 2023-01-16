#include "util.hpp"

/*****************************************************************************
 * 
 *                  Edited by Jiwon Choi Jan 6 2023
 *
 *      File description:
 *      Stochastic Series Expansion(SSE) for the one-dimensional random
 *      transverse-field Ising model.
 *      
 *      Reference: Anders W. Sandvik, PRE 68, 056701 (2003)
 * 
*****************************************************************************/

int main(int argc, char **argv)
{
  int i,n;
  int L = atoi(argv[1]);
  int M = atoi(argv[2]);
  double beta = atof(argv[3]);
  double Gmax = atof(argv[4]);
  int L2,nbond;
  double sumofJG,beta_sumofJG;
	L2 = L*L;
	nbond = 2*L2;
  n=0;


  int *spin = new int[L2];
  double *J = new double[nbond];
  double *G = new double[L2];
  double *CDtable = new double[L2+nbond];
  int *opstring = new int[M];
  int **bsites = new int*[2];
  bsites[0] = new int[nbond];
  bsites[1] = new int[nbond];

	int *first = new int[L2];	
	int *last = new int[L2];	
  int *vertex = new int[M];
	int *link = new int[4*M];
	int *stack = new int[4*M];

  //std::vector<int> *constterm = new std::vector<int>[L2];
  //std::vector<int> *flipterm = new std::vector<int>[L2];

  // Load interaction and transverse-field from file
  LoadInteractionFromFile("./Jconfig",2*L2,J);
  LoadInteractionFromFile("./Gconfig",L2,G);
  for (i=0;i<L2;++i) G[i] *= Gmax;
  // Initialize spin configuration
  for (i=0;i<L2;++i) spin[i] = (real(gen)>0.5 ? 1 : -1);
  for (i=0;i<M;++i) opstring[i] = -1;
  

  // Construct lattice structure
	for (int y=0;y<L;++y){
	  for (int x=0;x<L;++x){
			i = x+y*L;
			bsites[0][i] = i;
			bsites[1][i] = (i%L == L-1 ? i-L+1 : i+1);
			bsites[0][i+L2] = i;
			bsites[1][i+L2] = (i/L == L-1 ? i+L-L2 : i+L);
	}
	}
  
  // make the cumulant distribution table
  // 0 to nbond-1 : Ising term
  // nbond to nbond+L2-1 : gamma term
  sumofJG = 0.0;
  for (i=0;i<nbond;++i){
    CDtable[i] = J[i];
    sumofJG += J[i];
  }
  for (i=0;i<L2;++i){
    CDtable[i+nbond] = G[i];
    sumofJG += G[i];
  }
  for (i=0;i<nbond+L2-1;++i) CDtable[i+1] += CDtable[i];
  for (i=0;i<nbond+L2;++i) CDtable[i] /= sumofJG;

  beta_sumofJG = beta*sumofJG;
  
  // nbond to nbond+L2-1 : gamma term
  // nbond+L2 to nbond+2*L2-1 : Spin-flipping term
  for (int eq=0;eq<1000;++eq){
    DiagonalUpdate(L2,M,nbond,&n,spin,bsites,opstring,CDtable,beta_sumofJG);
    AdjustM(&M,n,opstring,vertex,link);
  }
  for (int eq=0;eq<1000;++eq){
		ConstructVertexAndLink(L2,M,n,nbond,bsites,opstring,vertex,link,first,last,stack);
  }
	//for (i=0;i<4*M;i+=4){
		//if (link[i] != -1)
		//std::cout << i << " " << link[link[i]] << "\n";
		//printf("%d %d %d %d\n",link[i],link[i+1],link[i+2],link[i+3]);
	//}
	//for (i=0;i<n;i++){
	//	std::cout << first[i] << " " << last[i] << "\n";
	//}
	
	
	
	
  return 0;
}

