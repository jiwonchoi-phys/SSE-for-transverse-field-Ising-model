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
	if (argc<9){
		std::cout << "Usage: ./2d L M beta Gmax Equilibration binsize repeat sampleidx\n";
		std::exit(1);
	}
  int i,n;
  int L = atoi(argv[1]);
  int M = atoi(argv[2]);
  double beta = atof(argv[3]);
  double Gmax = atof(argv[4]);
	int Equilibration = atoi(argv[5]);
	int binsize = atoi(argv[6]);
	int repeat = atoi(argv[7]);
	int sidx = atoi(argv[8]);
  int L2,nbond;
  double sumofJG,beta_sumofJG;
	double mag,mag2,mag4,magt0;
	double Mag,Mag2,Mag4,Magt0;
	char Jfile[80],Gfile[80],ObservableFile[80];
	L2 = L*L;
	nbond = 2*L2;
  n=0;

  char *spin = new char[L2];
  double *J = new double[nbond];
  double *G = new double[L2];
  double *CDtable = new double[L2+nbond];
  int *opstring = new int[M];
  int **bsites = new int*[2];
  bsites[0] = new int[nbond];
  bsites[1] = new int[nbond];

	int *first = new int[L2];	
	int *last = new int[L2];	
  char *vertex = new char[M];
	int *link = new int[4*M];
	int *stack = new int[4*M];
	char *visitedleg = new char[4*M];

  //std::vector<int> *constterm = new std::vector<int>[L2];
  //std::vector<int> *flipterm = new std::vector<int>[L2];

  // Load interaction and transverse-field from file
	sprintf(Jfile,"sample%d/Jconfig",sidx);
	sprintf(Gfile,"sample%d/Gconfig",sidx);
  LoadInteractionFromFile(Jfile,2*L2,J);
  LoadInteractionFromFile(Gfile,L2,G);
  for (i=0;i<L2;++i) G[i] *= Gmax;
  // Initialize spin configuration
  for (i=0;i<L2;++i) spin[i] = (real(gen)>0.5 ? 1 : -1);
  for (i=0;i<M;++i) opstring[i] = -1;
	
	//for (i=0;i<nbond;++i) std::cout << J[i] << "\n";

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

	// Equilibration
  for (int eq=0;eq<Equilibration;++eq){
    DiagonalUpdate(L2,M,nbond,&n,spin,bsites,opstring,CDtable,beta_sumofJG);
		ConstructVertexAndLink(L2,M,nbond,spin,bsites,opstring,vertex,link,first,last,stack);
		SwendsenWangUpdate(L2,M,nbond,opstring,link,visitedleg,stack,vertex);
		UpdateSpinAndOpstring(L2,M,nbond,spin,opstring,vertex,first,last);
    AdjustM(&M,n,opstring,vertex,link,stack,visitedleg);
		//std::cout << M << " " << n << "\n";
  }
	// Measurement
	sprintf(ObservableFile,"sample%d/L%dGmax%.4fdS",sidx,L,Gmax);
	std::ofstream output;
	output.open(ObservableFile,std::ios::app);
	for (int count=0;count<repeat;++count){
		Mag = Mag2 = Mag4 = Magt0 = 0.0;
		for (int bin=0;bin<binsize;++bin){
  	  DiagonalUpdate(L2,M,nbond,&n,spin,bsites,opstring,CDtable,beta_sumofJG);
			ConstructVertexAndLink(L2,M,nbond,spin,bsites,opstring,vertex,link,first,last,stack);
			SwendsenWangUpdate(L2,M,nbond,opstring,link,visitedleg,stack,vertex);
			//ClusterUpdate(L2,M,n,nbond,opstring,link,visitedleg,stack,vertex);
			UpdateSpinAndOpstring(L2,M,nbond,spin,opstring,vertex,first,last);
			//std::cout << M << " " << n << "\n";
			
			// Measurement (m,mt0,U,E,...)
			mag = mag2 = mag4 = magt0 = 0.0;
			for (i=0;i<L2;++i) magt0 += spin[i];
			Magt0 += fabs(magt0)/L2;
  	}
		Magt0 /= binsize;	
		output << beta << " " << Magt0 << "\n";
	}
	output.close();
	
  return 0;
}

