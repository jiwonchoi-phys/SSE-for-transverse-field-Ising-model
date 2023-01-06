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

void DiagonalUpdate(int L2, int M, int nbond, int *n, int *spin, int **bsites, int *opstring, double *CDtable, double aprob)
{    
  int nop, op, b;
  double dprob = 1/aprob;
  nop = L2+nbond;
  for (int m=0;m<M;++m){
    op = opstring[m];
    if (op == -1){
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


void AdjustM(int *M, int n, int* &opstring, int* &vertex, int* &link)
{
  int newM;
  // adjust opstring
  if (n < (int)(0.7*(*M))) return;

  newM = 1.2*(*M);
	// replace opstring
  int *newopstring = new int[newM];
  //std::cout << "new opsting size: " << sizeof(newopstring) << "\n";
  for (int i=0;i<(*M);++i) newopstring[i] = opstring[i];
  for (int i=*M;i<newM;++i) newopstring[i] = -1;

  delete[] opstring;
  opstring = newopstring;
  newopstring = NULL;
	
	// replace vertex
	int *newvertex = new int[newM];
	// For now, do not copy value of vertex
	vertex = newvertex;
	newvertex = NULL;

	// replace link
	int *newlink = new int[4*newM];
	link = newlink;
	newlink = NULL;
	
  *M = newM;
}

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

void ConstructVertexAndLink(int L2, int M, int nbond, int** bsites, int* opstring, int* vertex, int* link, int* first, int* last, int* flag)
{
	//int i,m,n,op,s1,s2;
	int m,s1,s2,v0,v1,v2;
	int op;
	n=0;
	// Initialize first and last
	for (i=0;i<L2;++i) first[i] = last[i] = -1;
	for (m=0;m<4*M;m+=4){
		op = opstring[m/4];
		if (op==-1){
			v1 = first[s1];
			if (v1 != -1){
				v2 = last[s1];
				vertex[v2] = v1;
				vertex[v1] = v2;
			}
		} else if (op<nbond){
			s1 = bsites[0][op];
			s2 = bsites[1][op];
			v1 = last[s1];
			v2 = last[s2];
			if (v1 != -1){
				vertex[v1] = v0;
				vertex[v0] = v1;
			} else first[s1] = v0;
			if (v2 != -1){
				vertex[v2] = v0+1;
				vertex[v0+1] = v2;
			} else first[s2] = v0+1;
			last[s1] = v0+2;
			last[s2] = v0+3;
		} else if (op<nbond+L2){
			// To do:
			//	 1. construct vertex and link for the constant and flip term
			
		}
}


