#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <iostream>
#include <cstdlib>
#include <random>
#include <vector>
#include <filesystem>
#include <fstream>
#include <algorithm>

namespace{
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> real(0,1);
}


void LoadInteractionFromFile(std::string filepath, int size, double *coupling);
int LowerBound(double *list, int length, double target);
void DiagonalUpdate(int L2, int M, int nbond, int *n, char *spin, int **bsites, int *opstring, double *CDtable, double aprob);
void AdjustM(int *M, int n, int* &opstring, char* &vertex, int* &link, int* &stack, char* &visitedleg);
void ConstructVertexAndLink(int L2, int M, int n, int nbond, char* spin, int** bsites, int* opstring, char* vertex, int* link, int* first, int* last, int* stack);
void LoopUpdate(int L2,int M,int n,int nbond,int *opstring,int *link,char *visitedleg,int *stack,char *vertex);
//void Partition(int L2, int M, int nbond, int *opstring, std::vector<gammaterm> *gammaseq);
void UpdateSpinAndOpstring(int L2, char* spin, int *opstring, char* vertex, int* first, int* last);
bool Visited(int leg, char* visitedleg, int* link);
#endif
