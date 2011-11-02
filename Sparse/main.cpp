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

#ifdef WINDOWS_OS
int main ()
{
#else
int main (int argc, char **argv)
{
#if 0
}  // fool most auto-indenters
#endif
#endif
	FILE *readFile; //the file to be read
	char line[100]; //used for input from the file
	char *tokenizer; //tokenizer used for input
	int numVertex, numFaces, nPoints; //the numbers of vertices and faces
	int **laplacian, *diagonal; //used for reading the laplacian
	int *faceP; //used to read points in a face
	bool bRelation;

	char fileName[100];
	char fullPath[150];
#ifdef WINDOWS_OS
	cin>>fileName;
	strcpy (fullPath,"c:\\test\\");
	strcat (fullPath, fileName);
	strcat (fullPath, ".off");
#else
        if (argc < 2) {
          fprintf(stderr, "oops, arguments!");
          return 1;
        }
        else {
          strcpy(fullPath, argv[1]);
          strcpy(fileName, argv[1]);
        }
#endif

	readFile = fopen (fullPath,"r"); //open the file for reading
	if(readFile==NULL)
			printf("Error opening the file");

	fgets(line, 100, readFile); //jump the first line (.OFF line)

	fgets(line, 100, readFile); //read the second line
	tokenizer = strtok (line, " "); //read the number of vertices
	numVertex=atoi(tokenizer);

	tokenizer = strtok (NULL," "); //read the number of faces
	numFaces=atoi(tokenizer);

	for(int i=0; i<numVertex; i++) //we jump over the next numVertex lines that contaign the coordinates of the points
		fgets(line, 100, readFile);

	//initialization
	diagonal=new int[numVertex];

	for(int i=0; i<numVertex; i++)
		diagonal[i]=0;

	laplacian=new int* [numVertex];

	for(int i=0; i<numVertex; i++)
	{
		laplacian[i]=new int[20];
		laplacian[i][0]=0;
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
				for(int k=0; k<diagonal[faceP[i]]; k++)
					if(faceP[j]==laplacian[faceP[i]][k])
						bRelation=false;

				//if the relation is not yet stored in the laplacian
				if(bRelation)
				{
					laplacian[faceP[i]][diagonal[faceP[i]]]=faceP[j];
					laplacian[faceP[j]][diagonal[faceP[j]]]=faceP[i];

					diagonal[faceP[i]]++;
					diagonal[faceP[j]]++; //add to the diagonal elements +1, for the extra relation the vertex has
				}
			}
		delete[] faceP;
	}
	//close the file
	fclose(readFile);

	real_2d_array  matrix;
	matrix.setlength(numVertex, numVertex);

	for(int i=0; i<numVertex; i++)
		for(int j=0; j<numVertex; j++)
			if(i==j)
				matrix[i][j]=diagonal[i];
			else
			{
				bRelation=false;
				for(int k=0; k<diagonal[i]; k++)
					if(laplacian[i][k]==j)
						bRelation=true;
				if(bRelation)
					matrix[i][j]=-1;
				else
					matrix[i][j]=0;
			}

	//for(int i=0; i<numVertex; i++)
			//cout<<matrix[i][i]<<" ";


	//use the function to compute the eigenvalues
	//*********************************************//

	//HeatKernelSignature variables
	real_1d_array eigenValues, w; //eigenvalues vector
	real_2d_array eigenVector; //eigenvector matrix
	double eSum=0; //will hold e^(-t*eigenValue)
	double tmin, tmax, step, t; //time variables
	double HKS[100];

	for(int i=0; i<100; i++)
		HKS[i]=0;

	//initializations
	eigenValues.setlength(numVertex);
	w.setlength(numVertex);
	for(int i=0; i<numVertex; i++)
		eigenValues[i]=0;

	eigenVector.setlength(numVertex,numVertex);


	cout<<"\n";

	smatrixevd(matrix, numVertex, 1, false, eigenValues, eigenVector);

	//for(int i=0; i<numVertex; i++)
		//cout<<eigenValues[i]<<" ";


	//compute the eigenvectors and calculate the coresponding heat kernel


	//initialize time variables
	tmin=4*log(10.0)/fabs(eigenValues[numVertex-1]);
	tmax=4*log(10.0)/fabs(eigenValues[numVertex-300]);
	step=(log(tmax)-log(tmin))/100;

	//compute the HKS for 100 t's
	for(int k=0; k<100; k++)
	{
		//compute the t
		t=tmin+exp(k*step);
		eSum=0;
		for(int i=0; i<numVertex; i++)
			w[i]=0;

		//e^(-t*eigenValue)
		for(int i=numVertex-300; i<numVertex; i++)
			eSum+=fabs(pow(2.71828,-eigenValues[i]*t));

		//in w we store each heatKernel for each vertex
		for(int i=0; i<numVertex; i++)
			for(int j=numVertex-300; j<numVertex; j++)
				w[i]+=eigenVector[i][j]*eigenVector[i][j]*pow(2.71828,-eigenValues[i]*t);

		//compute the heat kernel signature for each time t
		for(int i=0; i<numVertex; i++)
			HKS[k]+=w[i];//eSum;
	}

	FILE *writeFile; //the file to be read

#ifdef WINDOWS_OS
	strcpy (fullPath,"c:\\test\\signatures\\");
	strcat (fullPath, fileName);
	strcat (fullPath, ".txt");
#else
        strcpy (fullPath, fileName);
        strcat (fullPath, ".signature");
#endif
	writeFile = fopen (fullPath, "w"); //open the file for writing
	if(writeFile==NULL)
			printf("Error opening the file");

	for(int i=0; i<100; i++)
		fprintf(writeFile,"%f\n", HKS[i]);

	fclose(writeFile);

#ifdef WINDOWS_OS
        system("pause");
#endif
	fcloseall();
	return 0;
}
