/*  ALLEN DANIEL
 *  YESA
 *  ALLENDAN
 */

#ifndef A1_HPP
#define A1_HPP

#include <vector>
#include <mpi.h>
#include <algorithm>


int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
std::vector<signed char> B,C,M,Q,Z,D,E,F,N,G,H,I;
int b = n/q;
B.resize(b*b);
C.resize(b*b);//Equivalent to P
M.resize(b*b);//Helper for pointer jumping
Q.resize(b*b);
Z.resize(b*b);
D.resize(b*b,0);
E.resize(b*b);
F.resize(b*b);//Equivalent to P_Prime
N.resize(b*b,0);//Helper for tree hanging
G.resize(b*b);
H.resize(b*b);//Q matrix as per lecture
I.resize(b*b);
//getting the adjacent nodes
for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(A[i * b + j] != 0)
                                {B[i * b + j]=i;}
                        else
                                B[i * b + j]=A[i * b + j];
                }
        }

//inside every block finding max in column sequentially in parallel
for(int j=0;j<(b*b);j++)
{
	for(int k=j+b;k<(b*b);k=k+b)
	{	
		if(B[j]>B[k])
			B[k]=B[j];
		else
		B[j]=B[k];
	}
}


//Communication column wise and getting the max
int rank, size;


MPI_Comm_size(comm, &size);
MPI_Comm_rank(comm, &rank);


int col = rank % q;
int row = rank / q;


MPI_Comm col_comm;
MPI_Comm_split(comm, col, rank, &col_comm);


int nrank;
MPI_Comm_rank(col_comm, &nrank);



MPI_Comm row_comm;
MPI_Comm_split(comm, row, rank, &row_comm);


int mrank;
MPI_Comm_rank(row_comm, &mrank);

MPI_Allreduce(B.data(), C.data(), (b*b) , MPI_CHAR , MPI_MAX, col_comm);

do
{
//Creating helper matrix for pointer jumping
for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(A[i * b + j] == 1)
                                M[i * b +j] = C[i * b +j];
                        else
                                M[i * b +j] = 0;
                }
        }
    
//Find max in each row in each block sequentially
for (int k=0;k<(b*b);k=k+b)
{int max=0;
	for (int j=k;j<k+b;j++)	
	{
			if(max>M[j])
                		max=max;   
                	else
                		max=M[j];
	}
	for (int a=k;a<k+b;a++)
	Q[a]=max;
}
        

MPI_Allreduce(Q.data(), Z.data(), (b*b) , MPI_CHAR , MPI_MAX, row_comm);



for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(Z[i * b + j] == j)
                                {D[i * b +j]=C[i*b+j];}
			else
				D[i * b +j]=0;
                }
        }

//Each block sequentially finding max
for (int k=0;k<(b*b);k=k+b)
{int max=0;
        for (int j=k;j<k+b;j++)
        {
                        if(max>D[j])
                                max=max;
                        else
                                max=D[j];
        }
        for (int a=k;a<k+b;a++)
        E[a]=max;
}


MPI_Allreduce(E.data(), F.data(), (b*b) , MPI_CHAR , MPI_MAX, row_comm);


//Helper matrix for tree hanging
for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(C[i * b + j] == i)
                                N[i * b +j] = F[i * b +j];
                        else
                                N[i * b +j] = 0;
                }
        }

//Each block sequentially finding max
for (int k=0;k<(b*b);k=k+b)
{int max=0;
        for (int j=k;j<k+b;j++)
        {
                        if(max>N[j])
                                max=max;
                        else
                                max=N[j];
        }
        for (int a=k;a<k+b;a++)
        G[a]=max;
}



MPI_Allreduce(G.data(), H.data(), (b*b) , MPI_CHAR , MPI_MAX, row_comm);


//Final Comparison

for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(F[i * b + j] > H[i * b + j])
                                I[i * b + j]=F[i * b + j]; 
                	else
				I[i * b + j]=H[i * b + j];
		}
			
        }
//Comparing and changing row and colum communicators
if(I!=F)
{

//Creating helper matrix for pointer jumping
for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(A[i * b + j] == 1)
                                M[i * b +j] = C[i * b +j];
                        else
                                M[i * b +j] = 0;
                }
        }
    
//Find max in each column in each block sequentially


for(int j=0;j<b;j++)
{int max=0;
	for(int k=j;k<(b*b);k=k+b)
	{	
		if(M[k]>max)
			max=M[k];
		else
			max=max;
	}
	for(int a=j;a<(b*b);a=a+b)
	Q[a]=max;
}


        
MPI_Allreduce(Q.data(), Z.data(), (b*b) , MPI_CHAR , MPI_MAX, col_comm);



for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(Z[i * b + j] == j)
                                {D[i * b +j]=C[i*b+j];}
			else
				D[i * b +j]=0;
                }
        }



//Each block sequentially finding max
for(int j=0;j<b;j++)
{int max=0;
	for(int k=j;k<(b*b);k=k+b)
	{	
		if(D[k]>max)
			max=D[k];
		else
			max=max;
	}
	for(int a=j;a<(b*b);a=a+b)
	E[a]=max;
}



MPI_Allreduce(E.data(), F.data(), (b*b) , MPI_CHAR , MPI_MAX, col_comm);


//Helper matrix for tree hanging
for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(C[i * b + j] == i)
                                N[i * b +j] = F[i * b +j];
                        else
                                N[i * b +j] = 0;
                }
        }


//Each block sequentially finding max
for(int j=0;j<b;j++)
{int max=0;
	for(int k=j;k<(b*b);k=k+b)
	{	
		if(N[k]>max)
			max=N[k];
		else
			max=max;
	}
	for(int a=j;a<(b*b);a=a+b)
	G[a]=max;
}



MPI_Allreduce(G.data(), H.data(), (b*b) , MPI_CHAR , MPI_MAX, col_comm);


//Final Comparison

for (int i = 0; i < b; ++i)
        {
        for (int j = 0; j < b; ++j)
                {       if(F[i * b + j] > H[i * b + j])
                                I[i * b + j]=F[i * b + j]; 
                	else
				I[i * b + j]=H[i * b + j];
		}
			
        }

}
}while(I!=F);


std::sort(F.begin(), F.end());
int count = std::unique(F.begin(), F.end()) - F.begin();

MPI_Comm_free(&col_comm);
MPI_Comm_free(&row_comm);
return count;
} // connected_components

#endif // A1_HPP
