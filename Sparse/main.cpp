#include "os.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <math.h>
#include "linalg.h"

#ifdef WINDOWS_OS
#include <windows.h>
#endif

using namespace std;
using namespace alglib;

int main ()
{
	FILE *readFile; //the file to be read
	char line[100]; //used for input from the file
	char *tokenizer; //tokenizer used for input 
	int numVertex, numFaces, nPoints; //the numbers of vertices and faces
	int *laplacian[2], *diagonal; //used for reading the laplacian
	int *faceP; //used to read points in a face
	bool bRelation;
	int numRelations=0;

	//lancosz algorithm variables
	real_1d_array lapDiagonal, lapSubDiagonal; //the tridiagonilized laplacian matrix
	real_1d_array v, vPrev, w; //helpers
	real_2d_array some;
	float sumW, beta;

	//HeatKernelSignature variables
	real_1d_array eigenVector, eigenValues, x, hks;
	float eSum=0;

	readFile = fopen ("a.txt","r"); //open the file for reading
	if(readFile==NULL)
			printf("Eroare la deschiderea fisierului");
	
	fgets(line, 100, readFile); //jump the first line (.OFF line)
	
	fgets(line, 100, readFile); //read the second line
	tokenizer = strtok (line, " "); //read the number of vertices
	numVertex=atoi(tokenizer); 

	tokenizer = strtok (NULL," "); //read the number of faces
	numFaces=atoi(tokenizer);

	for(int i=0; i<numVertex; i++) //we jump over the next numVertex lines that contaign the coordinates of the points
	{
		fgets(line, 100, readFile);
		//in case we want to compute the matrix using the distances between the vertices
		/*
		tokenizer = strtok (line, " "); //read the x coordinate
		distanceV[i][0]=atoi(tokenizer); 
		tokenizer = strtok (NULL," "); //read the y coordinate
		distanceV[i][1]=atoi(tokenizer); 
		tokenizer = strtok (NULL," "); //read the z coordinate
		distanceV[i][2]=atoi(tokenizer); */
	}

	diagonal=new int[numVertex];

	for(int i=0; i<numVertex; i++)
		diagonal[i]=0;

	laplacian[0]=new int[numFaces*2];
	laplacian[1]=new int[numFaces*2];

	for(int i=0; i<numFaces*2; i++)
	{
		laplacian[0][i]=0;
		laplacian[1][i]=0;
	}

	//we construct the laplacian
	//*********************************************//
	for(int i=0; i<numFaces; i++)
	{
		fgets(line, 100, readFile);
		tokenizer=strtok(line," ");
		nPoints=atoi(tokenizer);  //we read the number of vertices that define the face
		
		faceP=new int[nPoints];   //initialize a vector that will keep the nPoints vertices
		
		//read the points in the faceP vector
		for(int i=0; i<nPoints; i++)
		{
			tokenizer=strtok(NULL," ");
			faceP[i]=atoi(tokenizer);
		}

		//we compare the points to the points stored in laplacian
		for(int i=0; i<nPoints-1; i++)
			for(int j=i+1; j<nPoints; j++)
			{
				bRelation=true;
				for(int k=0; k<numRelations; k++)
				{
					if(faceP[i]==laplacian[0][k])
						if(faceP[j]==laplacian[1][k])
							bRelation=false;
					if(faceP[i]==laplacian[1][k])
						if(faceP[j]==laplacian[0][k])
							bRelation=false;
				}
				//if the relation is not yet stored in the laplacian
				if(bRelation)
				{
					diagonal[faceP[i]]++;
					diagonal[faceP[j]]++; //add to the diagonal elements +1, for the extra relation the vertex has
					laplacian[0][numRelations]=faceP[i];
					laplacian[1][numRelations]=faceP[j]; //we add the two vertexes on the laplacian list for -1 elements
					numRelations++;
				}
			}
		delete[] faceP;
	}
	
	//tridiagonize the laplacian
	//*********************************************//
	//initializations
	lapDiagonal.setlength(numVertex);
	lapSubDiagonal.setlength(numVertex-1);
	v.setlength(numVertex);
	vPrev.setlength(numVertex);
	w.setlength(numVertex);
	some.setlength(2,2);

	for(int i=0; i<numVertex; i++)
		lapSubDiagonal[i]=0;
	for(int i=0; i<numVertex; i++)
	{
		vPrev[i]=0;
		v[i]=0;
		lapDiagonal[i]=0;
	}
	v[0]=1;
	beta=0;
	
	for(int k=0; k<numVertex; k++)
	{
		for(int i=0; i<numVertex; i++)
		{
			sumW=0;
			for(int j=0; j<numVertex; j++)
				for(int z=0; z<numRelations; z++)
					if(((laplacian[0][z]==i)&&(laplacian[1][z]==j))||((laplacian[0][z]==j)&&(laplacian[1][z]==i)))
						sumW-=v[j];
			sumW+=diagonal[i]*v[i];
			w[i]=sumW-beta*vPrev[i];
		}
		for(int i=0; i<numVertex; i++)
			lapDiagonal[k]+=w[i]*v[i];
		for(int i=0; i<numVertex; i++)
			w[i]=w[i]-lapDiagonal[k]*v[i];
		for(int i=0; i<numVertex; i++)
			lapSubDiagonal[k]+=w[i]*w[i];
		lapSubDiagonal[k]=sqrt(lapSubDiagonal[k]);
		beta=lapSubDiagonal[k];
		for(int i=0; i<numVertex; i++)
		{
			vPrev[i]=v[i];
			v[i]=w[i]/beta;
		}
	}

	
	//use the function to compute the eigenvalues
	//*********************************************//
	eigenValues.setlength(numVertex);
	for(int i=0; i<numVertex; i++)
		eigenValues[i]=lapDiagonal[i]; ///actually they are swapped
	smatrixtdevd(lapDiagonal,lapSubDiagonal,numVertex,0,some);
	
	//compute the eigenvectors and calculate the coresponding heat kernel
	//*********************************************//
	eigenVector.setlength(numVertex);
	hks.setlength(100);
	for(int i=0; i<100; i++)
		hks[i]=0;


	//x.setlength(numVertex);
	//for(int i=0; i<numVertex; i++)
		//x[i]=1;

	/*for(int j=0; j<numVertex; j++)
	{
		for(int k=1;k<5;k++)
		{
			bRelation=true;
			for (int i=0; i<numVertex; i++)
			{
				w[i]=eigenValues[i]-lapDiagonal[j];
				v[i]=x[i];
				eigenVector[i]=v[i];
			}

			for (int i=1; i<numVertex; i++)
			{
				double m = lapSubDiagonal[i-1]/w[i-1];
				w[i]=w[i]-m*lapSubDiagonal[i-1];
				v[i]=v[i]-m*v[i-1];
				cout<<w[i]<<" ";
			}

			cout<<"\n";

			x[numVertex-1]=v[numVertex-1]/w[numVertex-1];
			for(int i=numVertex-2; i>=0; i--)
				x[i]=(v[i]-lapSubDiagonal[i]*x[i+1])/w[i];
		}
	}*/

	float tmin=fabs(1/lapDiagonal[numVertex-1]);
	float tmax=fabs(1/lapDiagonal[0]);
	float step=(log(tmax)-log(tmin))/100;
	float t=tmin;

	for(int k=0; k<100; k++)
	{
		t=tmin+exp(k*step);
		eSum=0;

		for(int j=0; j<numVertex; j++)
		{
			eigenVector[0]=1;
			eigenVector[1]=((lapDiagonal[j]-eigenValues[0])*eigenVector[0])/lapSubDiagonal[0];
			for(int i=2; i<numVertex; i++)
				eigenVector[i]=(-lapSubDiagonal[i-2]*eigenVector[i-2]+eigenVector[i-1]*(lapDiagonal[j]-eigenValues[i-1]))/lapSubDiagonal[i-1];	
			
			w[0]+=eigenVector[0]*eigenVector[0]*pow(2.71828,-lapDiagonal[j]*t);
			for(int i=1; i<numVertex; i++)
				w[i]+=eigenVector[i]*eigenVector[i]*pow(2.71828,-lapDiagonal[j]*t);
		}

		for(int i=0; i<numVertex; i++)
			eSum+=pow(2.71828,-lapDiagonal[i]*t);

		for(int i=0; i<numVertex; i++)
			hks[k]+=w[i]/eSum;
	}

	

	for(int k=0; k<100; k++)
		cout<<hks[k]<<" ";


	int a;
	cin>>a;
	fcloseall();

}