#include "util.hpp"

// Use for load disorder configuration from file
void LoadInteractionFromFile(char* filename, int size, double *coupling)
{
  std::filesystem::path path_of_file(filename);
  if (!std::filesystem::exists(path_of_file)){
		////std::cout << "Configuration file does not exist.\n";
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
    } else if(op<nbond+L2){ // If constant op
      if (real(gen) < (M-(*n)+1)*dprob){
        opstring[m] = -1;
        *n -= 1;
      }
    // If operator is spin flip op
    } else spin[op-nbond-L2] *= -1;
  }
}


void AdjustM(int *M, int n, int* &opstring, char* &vertex, int* &link, int* &stack, char* &visitedleg)
{
  int newM;
  // adjust opstring
  if (n < (int)(0.85*(*M))) return;

  newM = 1.2*(*M);
	// replace opstring
  int *newopstring = new int[newM];
  for (int i=0;i<newM;++i) newopstring[i] = -1;
  for (int i=0;i<(*M);++i) newopstring[i] = opstring[i];

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
	int *newstack = new int[newM];
	stack = newstack;
	newstack = NULL;

	// replace link
	char *newvisitedleg = new char[4*newM];
	visitedleg = newvisitedleg;
	newvisitedleg = NULL;

  *M = newM;
}

void ConstructVertexAndLink(int L2, int M, int nbond, char* spin, int** bsites, int* opstring, char* vertex, int* link, int* first, int* last, int* stack)
{
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
    // If op is empty, skip it.
		if (op==-1) continue;
    // If op is Ising operator 
		if (op<nbond){
      // spins corresponding to bond operator 'op'.
			s1 = bsites[0][op];
			s2 = bsites[1][op];
      // previous leg acting on spin 's1' and 's2'.
			v1 = last[s1];
			v2 = last[s2];
      // construct link
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

    // If 'op' is constant operator
		} else if (op<nbond+L2){
      // now 'op' is spin site
			op -= nbond;
      // In cases of the contant and spin-flip operator,
      // there are only two leg.
			v1 = last[op];
			if (v1 != -1){
				link[v1] = v0;
				link[v0] = v1;
			} else first[op] = v0;
			last[op] = v0+2;
			if (spin[op]>0) vertex[p] = 2;
			else vertex[p] = 3;

    // If 'op' is spin-flip operator
		} else{
      // now 'op' is acting site
			op -= L2+nbond;
			v1 = last[op];
			if (v1 != -1){
				link[v1] = v0;
				link[v0] = v1;
			} else first[op] = v0;
			last[op] = v0+2;
			if (spin[op]>0) vertex[p] = 5;
			else vertex[p] = 4;
			spin[op] *= -1;
		}
	}
  // connect link across the time boundary
	for (i=0;i<L2;++i){
		v1 = first[i];
		if (v1!=-1){
			v2 = last[i];
			link[v2] = v1;
			link[v1] = v2;
		}
	}
}

void SwendsenWangUpdate(int L2,int M,int nbond,char *spin,int *opstring,int *link,char* &visitedleg,int *stack,char *vertex,int *first)
{
	int i,l,v,seedleg;
  // Initialize all legs to unvisited state.
  for (i=0;i<4*M;++i) visitedleg[i] = 0;
	// pick a seed leg
  // visit all spins
  for (l=0;l<L2;++l){
    // If spin[l] is not isolated
    if (first[l] > -1){
      // If spin[l] is not visited
      if (visitedleg[first[l]] == 0){
        // perform SW update
        //printf("spin%2d not visited\n",l);
	      std::vector<int> seedlist;
        seedleg = first[l];
        seedlist.push_back(seedleg);
        //std::cout << "start SW algorithm with seed leg: " << seedleg << "\n";
        // repeat single cluster algorithm until 'seedlist' become emtpy
	      while(seedlist.size()>0){	
	      	seedleg = seedlist.back();
	      	seedlist.pop_back();
          //std::cout << "start single cluster algorithm with seed leg: " << seedleg << "\n";
          //std::cout << "length of seedlist: " << seedlist.size() << "\n";

          // flip cluster with probability 1/2
	      	if (real(gen)<0.5) SingleClusterUpdate(L2,M,nbond,opstring,link,visitedleg,stack,vertex,seedlist,seedleg,true);
	      	else SingleClusterUpdate(L2,M,nbond,opstring,link,visitedleg,stack,vertex,seedlist,seedleg,false);
	      }
	      std::vector<int>().swap(seedlist);
      }
      // If spin[k] is isolated
    } else {if (real(gen)<0.5) spin[l]*=-1;}
  }
  //for (i=0;i<L2;++i) std::cout << (int)visitedleg[first[i]] << " ";
  //std::cout << "\n";
}

void SingleClusterUpdate(int L2,int M,int nbond,int *opstring,int *link,char* &visitedleg,int *stack,char *vertex,std::vector<int> &seedlist,int leg,bool AllowFlip)
{
	int i,l,pos,legidx;
	int sp=0;
	int spmax=1;

  // start with given seed leg
	
  // add seed leg to stack
	stack[sp] = leg;
  sp++;
	while(sp){
    //std::cout << "selected leg " << leg << "\n";
    sp--;
		leg = stack[sp];
    // If leg have been visited, ignore it.
		if (visitedleg[leg] == 1) continue;
    //std::cout << sp << " " << seedlist.size() << "\n";
		legidx = leg%4;
		pos = leg/4;
		//printf("pos: %d, legidx: %d\n",pos,legidx);
    // mark this leg to be visited.
		visitedleg[leg] = 1;
		
		// If vertex is Ising operator	
		if (vertex[pos]/2 == 0){
      //std::cout << "Ising\n";
      // all legs corresponding to this vertex will not be visited anymore.
		  for (l=0;l<4;++l){
        visitedleg[4*pos+l] = 1;
        // if linked leg was not visited yet, add it to stack
        if (visitedleg[link[4*pos+l]] == 0){
          stack[sp++] = link[4*pos+l];
          //std::cout << "current link: " << 4*pos+l << ", added to stack: " << link[4*pos+l] << "\n";
        }
		  }
      // if spin-flip is allowed, then change type of vertex
			if (AllowFlip) vertex[pos] = (vertex[pos]+1)%2;

    // If vertex is constant operator 
		} else if (vertex[pos]/2 == 1){ 
      //std::cout << "Constant\n";
      // if linked leg is not visited, add it to stack
			if (visitedleg[link[leg]] == 0) stack[sp++] = leg;
      // if leg from lower side
			if (legidx==0){
				if (visitedleg[4*pos+2] == 0) seedlist.push_back(4*pos+2);
				if (AllowFlip){
					if (vertex[pos] == 2) vertex[pos] = 4;
					else vertex[pos] = 5;
				}
      // if leg from upper side
			} else{
				if (visitedleg[4*pos] == 0) seedlist.push_back(4*pos);
				if (AllowFlip){
					if (vertex[pos] == 2) vertex[pos] = 5;
					else vertex[pos] = 4;
				}
			}

		// If vertex is spin-flip operator
		} else{
      //std::cout << "Flip\n";
      // if linked leg is not visited, add it to stack
			if (visitedleg[link[leg]] == 0) stack[sp++] = leg;
      // if leg from lower side
			if (legidx == 0){
				if (visitedleg[4*pos+2] == 0) seedlist.push_back(4*pos+2);
				if (AllowFlip){
					if (vertex[pos] == 4) vertex[pos] = 2;
					else vertex[pos] = 3;
				}
      // if leg from upper side
			} else{
				if (visitedleg[4*pos] == 0) seedlist.push_back(4*pos);
				if (AllowFlip){
					if (vertex[pos] == 4) vertex[pos] = 3;
					else vertex[pos] = 2;
				}
			}
		}
		//std::cout << "next leg(left: " << sp <<")\n";
	}
}

void UpdateSpinAndOpstring(int L2, int M, int nbond, char* spin, int *opstring, char* vertex, int* first, int* last)
{
	int i,m;
	// update spin
	for (i=0;i<L2;++i){
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

