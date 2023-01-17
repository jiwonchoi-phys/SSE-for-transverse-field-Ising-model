#include "util.hpp"

// Use for load disorder configuration from file
void LoadInteractionFromFile(std::string filename, int size, double *coupling)
{
  std::filesystem::path path_of_file(filename);
  if (std::filesystem::exists(path_of_file)){
    // If file exists, open the file 
    std::ifstream filestream(path_of_file);
    char temp[30];
    filestream.getline(temp,sizeof(temp));
    if (atoi(temp) != size){
      std::cout << "parameter 'size' and file size are not matching!\n";
      std::exit(1);
    } else{
      for (int i=0;i<size;++i){
        filestream.getline(temp,sizeof(temp));
        coupling[i] = atof(temp);
      }
    }
    filestream.close();
  } else{
    // If file does not exist, write new coupling on the file
    std::ofstream filestream(path_of_file);
    filestream << size << "\n";
    double temp;
    for (int i=0;i<size;++i){
      temp = real(gen);
      coupling[i] = temp;
      filestream << temp << "\n";
    }
    filestream.close();
  }
}

// Use for finding the first index from cumulative distribution table
int LowerBound(double *list, int length, double target){
  int start = 0;
  int end = length;
  int mid = 0;
  while (start < end){
    mid = (start + end)/2;
    if (list[mid] > target) end = mid;
    else start = mid + 1;
  }
  return start; 
}

void DiagonalUpdate(int L2, int M, int nbond, int *n, char *spin, int **bsites, int *opstring, double *CDtable, double aprob)
{    
  int nop, op, b;
  double dprob = 1/aprob;
  nop = L2+nbond;
  for (int m=0;m<M;++m){
    op = opstring[m];
    if (op == -1){ // If empty, add new operator
      if (real(gen)*(M-(*n))<aprob){
        b = LowerBound(CDtable,nop,real(gen));
        if (b<nbond){
          if (spin[bsites[0][b]] == spin[bsites[1][b]]){
            opstring[m] = b;
            *n += 1;
          }
        } else{
          opstring[m] = b;
          *n += 1;
        }
      }
    } else if(op<nbond+L2){
      if (real(gen) < (M-(*n)+1)*dprob){
        opstring[m] = -1;
        *n -= 1;
      }
    } else spin[op-nbond-L2] *= -1;
  }
}


void AdjustM(int *M, int n, int* &opstring, char* &vertex, int* &link, int* &stack, char* &visitedleg)
{
  int newM;
  // adjust opstring
  if (n < (int)(0.7*(*M))) return;

  newM = 1.2*(*M);
	// replace opstring
  int *newopstring = new int[newM];
  for (int i=0;i<(*M);++i) newopstring[i] = opstring[i];
  for (int i=*M;i<newM;++i) newopstring[i] = -1;

  delete[] opstring;
  opstring = newopstring;
  newopstring = NULL;
	
	// replace vertex
	char *newvertex = new char[newM];
	// For now, do not copy value of vertex
	vertex = newvertex;
	newvertex = NULL;

	// replace link
	int *newlink = new int[4*newM];
	link = newlink;
	newlink = NULL;
	
	// replace link
	int *newstack = new int[4*newM];
	stack = newstack;
	newstack = NULL;

	// replace link
	char *newvisitedleg = new char[4*newM];
	visitedleg = newvisitedleg;
	newvisitedleg = NULL;

  *M = newM;
}
/*
void Partition(int L2, int M, int nbond, int *opstring, std::vector<gammaterm> *gammaseq)
{
	int op;
	// clear subsequence
	for (int i=0;i<L2;++i){
		if (gammaseq[i].size()>0) std::vector<gammaterm>().swap(gammaseq[i]);
	}

	// divide opstring into subsequences
  for (int m=0;m<M;++m){
		op = opstring[m];	
    if (op >= nbond) gammaseq[(op-nbond)%L2].push_back(gammaterm(m,(op-nbond)/L2));
  }
}
*/

void ConstructVertexAndLink(int L2, int M, int n, int nbond, char* spin, int** bsites, int* opstring, char* vertex, int* link, int* first, int* last, int* stack)
{
	//int i,m,n,op,s1,s2;
	int i,s1,s2,v0,v1,v2,p;
	int op,sp;
	// Initialize first and last
	for (i=0;i<L2;++i) first[i] = last[i] = -1;
	for (i=0;i<4*M;++i) link[i] = -1;
	for (i=0;i<M;++i) vertex[i] = -1;

	// Construct link and vertex
	for (v0=0;v0<4*M;v0+=4){
		p = v0/4;
		op = opstring[p];
		if (op==-1){
			// ACtion for empty operator
		} else if (op<nbond){
			// Ising term
			s1 = bsites[0][op];
			s2 = bsites[1][op];
			v1 = last[s1];
			v2 = last[s2];
			if (v1 != -1){
				link[v1] = v0;
				link[v0] = v1;
			} else first[s1] = v0;
			if (v2 != -1){
				link[v2] = v0+1;
				link[v0+1] = v2;
			} else first[s2] = v0+1;
			last[s1] = v0+2;
			last[s2] = v0+3;
			//Assign vertex index
			if (spin[s1]>0) vertex[p] = 0;
			else vertex[p] = 1;

		} else if (op<nbond+L2){
			// constant term
			op -= nbond;
			v1 = last[op];
			if (v1 != -1){
				link[v1] = v0;
				link[v0] = v1;
			} else first[op] = v0+1;
			last[op] = v0+2;
			if (spin[op]>0) vertex[p] = 2;
			else vertex[p] = 3;

		} else{
			// flip term
			op -= L2+nbond;
			v1 = last[op];
			if (v1 != -1){
				link[v1] = v0;
				link[v0] = v1;
			} else first[op] = v0+1;
			last[op] = v0+2;
			if (spin[op]>0) vertex[p] = 4;
			else vertex[p] = 5;
			spin[op] *= -1;
		}
	}
	for (i=0;i<L2;++i){
		v1 = first[i];
		if (v1!=-1){
			v2 = last[i];
			link[v2] = v1;
			link[v1] = v2;
		}
	}
}

void LoopUpdate(int L2,int M,int n,int nbond,int *opstring,int *link,char *visitedleg,int *stack,char *vertex)
{
	int i,v,leg,l,pos;
	int sp=1;
	int spmax=1;
	for (i=0;i<4*M;++i) visitedleg[i] = 0;
	// pick one of the legs
	v = M*real(gen);
	while(opstring[v] == -1) v = M*real(gen);
	if (v<nbond) leg = 4*v+4*real(gen);
	else	leg = 4*v+(real(gen)<0.5 ? 2 : 0);
	stack[0] = leg;

	while(sp){
		leg = stack[--sp];
		pos = leg/4;
		if (visitedleg[leg] == 1) continue;
		visitedleg[leg] = 1;
		if (sp>spmax) spmax = sp;
		
		// If vertex is Ising type	
		if (vertex[pos]/2 == 0){
			// add all leg to stack
			for (l=0;l<4;++l){
				if (visitedleg[link[4*pos+l]] == 0) stack[sp++] = link[4*pos+l];
			}
			// change vertex type
			vertex[pos] = (vertex[pos]+1)%2;
		// If vertex is constant type
		} else if (vertex[pos]/2 == 1){
			if (visitedleg[link[leg]] == 0) stack[sp++] = link[leg];
			if (vertex[pos] == 2) vertex[pos] = 5;
			else vertex[pos] = 4;
		// If vertex is flip type
		} else{
			if (visitedleg[link[leg]] == 0) stack[sp++] = link[leg];
			if (vertex[pos] == 4) vertex[pos] = 3;
			else vertex[pos] = 2;
		}

		visitedleg[leg] = 1;
	}
	std::cout << spmax << "\n";
}

void UpdateSpinAndOpstring



