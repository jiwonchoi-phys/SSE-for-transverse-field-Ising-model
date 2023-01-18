#include "util.hpp"

// Use for load disorder configuration from file
void LoadInteractionFromFile(char* filename, int size, double *coupling)
{
  std::filesystem::path path_of_file(filename);
  if (!std::filesystem::exists(path_of_file)){
		std::cout << "Configuration file does not exist.\n";
		std::exit(1);
  } else{
    std::ifstream inputfile(path_of_file);
    //char temp[60];
		std::string temp;
    for (int i=0;i<size;++i){
			//inputfile.getline(temp);
			std::getline(inputfile,temp);
			coupling[i] = std::stod(temp);
    }
    inputfile.close();
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

void ConstructVertexAndLink(int L2, int M, int nbond, char* spin, int** bsites, int* opstring, char* vertex, int* link, int* first, int* last, int* stack)
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

void SwendsenWangUpdate(int L2,int M,int nbond,int *opstring,int *link,char *visitedleg,int *stack,char *vertex)
{
	int i,v,leg,nseed;
	std::vector<int> seedlist;
	// pick a seed leg
	v = M*real(gen);
	while(opstring[v] == -1) v = M*real(gen);
	if (v<nbond) leg = 4*v+4*real(gen);
	else	leg = 4*v+(real(gen)<0.5 ? 2 : 0);
	for (i=0;i<4*M;++i) visitedleg[i] = 0;
	seedlist.push_back(leg);

	while(seedlist.size()>0){	
		leg = seedlist.back();
		seedlist.pop_back();
		if (real(gen)<0.5) SingleClusterUpdate(L2,M,nbond,opstring,link,visitedleg,stack,vertex,seedlist,leg,true);
		else SingleClusterUpdate(L2,M,nbond,opstring,link,visitedleg,stack,vertex,seedlist,leg,false);
	}
	std::vector<int>().swap(seedlist);
}

void SingleClusterUpdate(int L2,int M,int nbond,int *opstring,int *link,char *visitedleg,int *stack,char *vertex,std::vector<int> &seedlist,int leg,bool AllowFlip)
{
	int i,l,pos,legidx;
	int sp=0;
	int spmax=1;
	// pick one of the legs
	
	// start cluster spanning process
	visitedleg[leg] = 1;
	pos = leg/4;
	legidx = leg%4;
	stack[sp] = leg;
	sp++;
	if (vertex[pos]/2 == 0){
		for (l=0;l<4;++l){
			if (visitedleg[4*pos+l]==0){
				stack[sp] = 4*pos+l;
				sp++;
				//if (sp>spmax) spmax=sp;
			}
		}
		// update vertex type
		vertex[pos] = (vertex[pos]+1)%2;
	} else if (vertex[pos]/2 == 1){
		if (legidx==0){
//////////////////////////// Should add new seed to SW list
			if (visitedleg[4*pos+2] == 0) seedlist.push_back(4*pos+2);
			if (AllowFlip){
				if (vertex[pos] == 2) vertex[pos] = 4;
				else vertex[pos] = 5;
			}
		} else{
			if (visitedleg[4*pos] == 0) seedlist.push_back(4*pos);
			if (AllowFlip){
				if (vertex[pos] == 2) vertex[pos] = 5;
				else vertex[pos] = 4;
			}
		}
	} else{
		if (legidx==0){
			if (visitedleg[4*pos+2] == 0) seedlist.push_back(4*pos+2);
			if (AllowFlip){
				if (vertex[pos] == 4) vertex[pos] = 2;
				else vertex[pos] = 3;
			}
		} else{
			if (visitedleg[4*pos] == 0) seedlist.push_back(4*pos);
			if (AllowFlip){
				if (vertex[pos] == 4) vertex[pos] = 3;
				else vertex[pos] = 2;
			}
		}
	}

	while(sp){
		leg = link[stack[--sp]];
		if (visitedleg[leg] == 1) continue;
		legidx = leg%4;
		pos = leg/4;
		visitedleg[leg] = 1;
		
		// If vertex is Ising type	
		if (vertex[pos]/2 == 0){
			// add all leg to stack
			for (l=0;l<4;++l){
				if (visitedleg[4*pos+l]==0){
					stack[sp] = 4*pos+l;
					sp++;
					//if (sp>spmax) spmax=sp;
				}
			}
			// update vertex type
			vertex[pos] = (vertex[pos]+1)%2;
		} else if (vertex[pos]/2 == 1){ // If vertex is constant type
			//if (visitedleg[leg] == 0) stack[sp++] = link[leg];
			if (legidx==0){
				if (visitedleg[4*pos+2] == 0) seedlist.push_back(4*pos+2);
				if (AllowFlip){
					if (vertex[pos] == 2) vertex[pos] = 4;
					else vertex[pos] = 5;
				}
			} else{
				if (visitedleg[4*pos] == 0) seedlist.push_back(4*pos);
				if (AllowFlip){
					if (vertex[pos] == 2) vertex[pos] = 5;
					else vertex[pos] = 4;
				}
			}
		// If vertex is flip type
		} else{
			//if (visitedleg[leg] == 0) stack[sp++] = link[leg];
			if (legidx == 0){
				if (visitedleg[4*pos+2] == 0) seedlist.push_back(4*pos+2);
				if (AllowFlip){
					if (vertex[pos] == 4) vertex[pos] = 2;
					else vertex[pos] = 3;
				}
			} else{
				if (visitedleg[4*pos] == 0) seedlist.push_back(4*pos);
				if (AllowFlip){
					if (vertex[pos] == 4) vertex[pos] = 3;
					else vertex[pos] = 2;
				}
			}
		}
	}
}

void UpdateSpinAndOpstring(int L2, int M, int nbond, char* spin, int *opstring, char* vertex, int* first, int* last)
{
	int i,m;
	// update spin
	for (i=0;i<L2;++i){
		//std::cout << first[i] << " " << (int)vertex[first[i]/4] << " " << (int)spin[i] << " ";
		switch (vertex[first[i]/4])
		{
			case 0: spin[i] = 1; break;
			case 1: spin[i] = -1; break;
			case 2: spin[i] = 1; break;
			case 3: spin[i] = -1; break;
			case 4: spin[i] = -1; break;
			case 5: spin[i] = 1; break;
		}
	}

	// update opstring
	for (m=0;m<M;++m){
		if (vertex[m]<2) continue;
		// constant type operator
		if (vertex[m]/2 == 1) opstring[m] = nbond+(opstring[m]-nbond)%L2;
		// flip type operator
		else opstring[m] = nbond+L2+(opstring[m]-nbond)%L2;
	}
}

