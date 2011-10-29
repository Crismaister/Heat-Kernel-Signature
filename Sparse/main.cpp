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
	float *diagonal, **sValues;
	int **laplacian, *nElements; //used for reading the laplacian
	int faceP[3]; //used to read points in a face
	real_2d_array pCoords;
	int bRelation;
	float sA, sB, sC, A, B, C, *area, semi; //sides of the triangles and angles

	char fileName[100];
	char fullPath[150];
#ifdef WINDOWS_OS
	cin>>fileName;
	strcpy (fullPath,"c:\\test\\");
	strcat (fullPath, fileName);
	strcat (fullPath, ".off");
#else
        if (argc < 2) {
          fprintf(stderr, "oops, arguments!\n");
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

	pCoords.setlength(numVertex, 3);
	for(int i=0; i<numVertex; i++) //we read the coords of the points in pCoords
	{
		fgets(line, 100, readFile);
		tokenizer = strtok (line, " ");
		pCoords[i][0]=atof(tokenizer); //read the x coordinate

		tokenizer = strtok (NULL," ");
		pCoords[i][1]=atof(tokenizer); //read the y coordinate

		tokenizer = strtok (NULL," ");
		pCoords[i][2]=atof(tokenizer); //read the z coordinate
	}

	//initialization
	diagonal=new float[numVertex];
	nElements=new int[numVertex];

	for(int i=0; i<numVertex; i++)
	{
		diagonal[i]=0;
		nElements[i]=0;
	}

	laplacian=new int* [numVertex];
	sValues=new float*[numVertex];

	for(int i=0; i<numVertex; i++)
	{
		laplacian[i]=new int[50];

		sValues[i]=new float[50];
	}

	for(int i=0; i<numVertex; i++)
		for(int j=0; j<50; j++)
			laplacian[i][j]=-1;

	area=new float[numVertex];
	for(int i=0; i<numVertex; i++)
		area[i]=0;

	//we construct the laplacian
	//*********************************************//
	for(int i=0; i<numFaces; i++)
	{
		fgets(line, 100, readFile);
		tokenizer=strtok(line," ");
		nPoints=atoi(tokenizer);  //we read the number of vertices that define the face

		//read the points in the faceP vector
		for(int i=0; i<3; i++)
		{
			tokenizer=strtok(NULL," ");
			faceP[i]=atoi(tokenizer);
		}

		//compute the sides of the triangle
		sA=sqrt((pCoords[faceP[1]][0]-pCoords[faceP[2]][0])*(pCoords[faceP[1]][0]-pCoords[faceP[2]][0])+
				(pCoords[faceP[1]][1]-pCoords[faceP[2]][1])*(pCoords[faceP[1]][1]-pCoords[faceP[2]][1])+
				(pCoords[faceP[1]][2]-pCoords[faceP[2]][2])*(pCoords[faceP[1]][2]-pCoords[faceP[2]][2]));
		sB=sqrt((pCoords[faceP[0]][0]-pCoords[faceP[2]][0])*(pCoords[faceP[0]][0]-pCoords[faceP[2]][0])+
				(pCoords[faceP[0]][1]-pCoords[faceP[2]][1])*(pCoords[faceP[0]][1]-pCoords[faceP[2]][1])+
				(pCoords[faceP[0]][2]-pCoords[faceP[2]][2])*(pCoords[faceP[0]][2]-pCoords[faceP[2]][2]));
		sC=sqrt((pCoords[faceP[0]][0]-pCoords[faceP[1]][0])*(pCoords[faceP[0]][0]-pCoords[faceP[1]][0])+
				(pCoords[faceP[0]][1]-pCoords[faceP[1]][1])*(pCoords[faceP[0]][1]-pCoords[faceP[1]][1])+
				(pCoords[faceP[0]][2]-pCoords[faceP[1]][2])*(pCoords[faceP[0]][2]-pCoords[faceP[1]][2]));

		//compute the areas of the triangle needed for the diagonal
		semi=(sA+sB+sC)/2;
		area[faceP[0]]+=sqrt(semi*(semi-sA)*(semi-sB)*(semi-sC));
		area[faceP[1]]+=sqrt(semi*(semi-sA)*(semi-sB)*(semi-sC));
		area[faceP[2]]+=sqrt(semi*(semi-sA)*(semi-sB)*(semi-sC));

		//compute the cos of the angles between the sides
		A=(sB*sB+sC*sC-sA*sA)/(2*sB*sC);
		B=(sA*sA+sC*sC-sB*sB)/(2*sA*sC);
		C=(sA*sA+sB*sB-sC*sC)/(2*sA*sB);
		if(A<0.0001)
			A=0;
		if(B<0.0001)
			B=0;
		if(C<0.0001)
			C=0;

		//first side
		bRelation=-1;
		for(int k=0; k<nElements[0]+1; k++)
			if(faceP[1]==laplacian[faceP[0]][k])
				bRelation=k;

		if(bRelation==-1)
		{
			laplacian[faceP[0]][nElements[faceP[0]]]=faceP[1];
			laplacian[faceP[1]][nElements[faceP[1]]]=faceP[0]; //add the elements that become non 0 on the sparse matrix

			sValues[faceP[0]][nElements[faceP[0]]]=C/sqrt(1-C*C);
			sValues[faceP[1]][nElements[faceP[1]]]=C/sqrt(1-C*C); //add the value of the element on the respective position of the sparse matrix

			nElements[faceP[0]]++;
			nElements[faceP[1]]++;
		}else{
			sValues[faceP[0]][bRelation]+=C/sqrt(1-C*C);

			for(int k=0; k<nElements[0]+1; k++)
				if(faceP[0]==laplacian[faceP[1]][k])
					bRelation=k;

			sValues[faceP[1]][bRelation]+=C/sqrt(1-C*C); //add the other angle to the respective poasition of the sparse matrix
		}

		//second side
		bRelation=-1;
		for(int k=0; k<nElements[0]+1; k++)
			if(faceP[2]==laplacian[faceP[0]][k])
				bRelation=k;

		if(bRelation==-1)
		{
			laplacian[faceP[0]][nElements[faceP[0]]]=faceP[2];
			laplacian[faceP[2]][nElements[faceP[2]]]=faceP[0]; //add the elements that become non 0 on the sparse matrix

			sValues[faceP[0]][nElements[faceP[0]]]=B/sqrt(1-B*B);
			sValues[faceP[2]][nElements[faceP[2]]]=B/sqrt(1-B*B); //add the value of the element on the respective position of the sparse matrix

			nElements[faceP[0]]++;
			nElements[faceP[2]]++;
		}else{
			sValues[faceP[0]][bRelation]+=B/sqrt(1-B*B);

			for(int k=0; k<nElements[0]+1; k++)
				if(faceP[0]==laplacian[faceP[2]][k])
					bRelation=k;

			sValues[faceP[2]][bRelation]+=B/sqrt(1-B*B); //add the other angle to the respective poasition of the sparse matrix
		}

		//third side
		bRelation=-1;
		for(int k=0; k<nElements[1]+1; k++)
			if(faceP[2]==laplacian[faceP[1]][k])
				bRelation=k;

		if(bRelation==-1)
		{
			laplacian[faceP[1]][nElements[faceP[1]]]=faceP[2];
			laplacian[faceP[2]][nElements[faceP[2]]]=faceP[1]; //add the elements that become non 0 on the sparse matrix

			sValues[faceP[1]][nElements[faceP[1]]]=A/sqrt(1-A*A);
			sValues[faceP[2]][nElements[faceP[2]]]=A/sqrt(1-A*A); //add the value of the element on the respective position of the sparse matrix

			nElements[faceP[1]]++;
			nElements[faceP[2]]++;
		}else{
			sValues[faceP[1]][bRelation]+=A/sqrt(1-A*A);

			for(int k=0; k<nElements[1]+1; k++)
				if(faceP[1]==laplacian[faceP[2]][k])
					bRelation=k;

			sValues[faceP[2]][bRelation]+=A/sqrt(1-A*A); //add the other angle to the respective poasition of the sparse matrix
		}
	}

	//calculate the diagonal
	for(int i=0; i<numVertex; i++)
		for(int j=0; j<nElements[i]; j++)
			diagonal[i]+=(2*sValues[i][j])/area[i];


	//close the file
	fclose(readFile);

	//tridiagonize the laplacian using the lancosz algorithm
	//*********************************************//
	real_1d_array lapDiagonal, lapSubDiagonal; //the tridiagonilized laplacian matrix
	real_1d_array v, vPrev, w; //helpers
	float sumW, beta;

	//initializations
	lapDiagonal.setlength(numVertex);
	lapSubDiagonal.setlength(numVertex-1);
	v.setlength(numVertex);
	vPrev.setlength(numVertex);
	w.setlength(numVertex);

	for(int i=0; i<numVertex; i++)
		lapSubDiagonal[i]=0;

	//v=v[1]-random vector with norm 1, vprev=v[0]-0 vector
	for(int i=0; i<numVertex; i++)
	{
		vPrev[i]=0;
		v[i]=0;
		lapDiagonal[i]=0;
	}
	v[0]=1;
	//beta-first sub/superdiagonal element does not exist (there are only n-1 sub/superdiagonal elements compared to n diagonal elements)
	beta=0;

	for(int k=0; k<numVertex; k++)
	{
		//compute the w vector
		for(int i=0; i<numVertex; i++)
		{
			sumW=0;
			for(int j=0; j<nElements[i]; j++)
					sumW+=sValues[i][j]*v[laplacian[i][j]];
			sumW+=diagonal[i]*v[i];
			w[i]=sumW-beta*vPrev[i];
		}
		//compute the diagonal element from the w*v
		for(int i=0; i<numVertex; i++)
			lapDiagonal[k]+=w[i]*v[i];
		//finish computing the w vector
		for(int i=0; i<numVertex; i++)
			w[i]=w[i]-lapDiagonal[k]*v[i];
		//compute the sub/superdiagonal elements
		for(int i=0; i<numVertex; i++)
			lapSubDiagonal[k]+=w[i]*w[i];
		lapSubDiagonal[k]=sqrt(lapSubDiagonal[k]);
		//put the newly sub/superdiagonal element in beta for the next iteration
		beta=lapSubDiagonal[k];
		//vprev=v and compute the new v
		for(int i=0; i<numVertex; i++)
		{
			vPrev[i]=v[i];
			v[i]=w[i]/beta;
		}
	}



	//release the memory of elements that we don't need
	laplacian=NULL;
	diagonal=NULL;
	cout<<"here";
	//use the function to compute the eigenvalues
	//*********************************************//

	//HeatKernelSignature variables
	real_1d_array eigenValues; //eigenvalues vector
	real_1d_array eigenVector; //eigenvector matrix
	double eSum=0; //will hold e^(-t*eigenValue)
	double tmin, tmax, step, t; //time variables
	double HKS[100];

	real_2d_array ceva;

	for(int i=0; i<100; i++)
			HKS[i]=0;

	//initializations
	eigenValues.setlength(numVertex);
	for(int i=0; i<numVertex; i++)
		eigenValues[i]=lapDiagonal[i];

	smatrixtdevd(eigenValues, lapSubDiagonal, numVertex, 0, ceva);

	//compute the eigenvectors and calculate the coresponding heat kernel
	for(int i=0; i<numVertex; i++)
		w[i]=0;

	//initialize time variables
	tmin=fabs(4*log(10.0)/eigenValues[numVertex-1]);
	tmax=fabs(4*log(10.0)/eigenValues[0]);
	step=(log(tmax)-log(tmin))/100;
	/*
	//compute the HKS for 100 t's
	eigenVector.setlength(numVertex,numVertex);

	for(int k=0; k<100; k++)
	{
		//compute the t
		t=tmin+exp(k*step);
		eSum=0;

		//e^(-t*eigenValue)
		for(int i=0; i<numVertex; i++)
			eSum+=fabs(pow(2.71828,-eigenValues[i]*t));

		//in w we store each heatKernel for each vertex
		for(int i=0; i<numVertex; i++)
			for(int j=0; j<numVertex; j++)
				w[i]+=fabs(eigenVector[i][j])*fabs(eigenVector[i][j])*fabs(pow(2.71828,-eigenValues[i]*t));

		//compute the heat kernel signature for each time t
		for(int i=0; i<numVertex; i++)
			HKS[k]+=w[i]/eSum;
	}
	*/
	//eigenvector not computed using the function
	eigenVector.setlength(numVertex);

	for(int k=0; k<10; k++)
	{
		//compute the t
		t=tmin+exp(k*step);
		eSum=0;

		//e^(-t*eigenValue)
		for(int i=0; i<numVertex; i++)
			eSum+=pow(2.71828,-eigenValues[i]*t);

		//compute the eigenvector coresponding to j'th-eigenvalue
		for(int j=1; j<numVertex; j++)
		{
			eigenVector[0]=pow(10.0,-1000);
			eigenVector[1]=((eigenValues[j]-lapDiagonal[0])*eigenVector[0])/lapSubDiagonal[0];
			for(int i=2; i<numVertex; i++)
				eigenVector[i]=(-lapSubDiagonal[i-2]*eigenVector[i-2]+eigenVector[i-1]*(eigenValues[j]-lapDiagonal[i-1]))/lapSubDiagonal[i-1];

			//in w we store each heatKernel for each vertex
			for(int i=0; i<numVertex; i++)
				w[i]+=eigenVector[i]*eigenVector[i]*pow(2.71828,-lapDiagonal[j]*t);
		}

		//compute the heat kernel signature
		for(int i=0; i<numVertex; i++)
			HKS[k]+=w[i]/eSum;

	}
	for(int i=0; i<numVertex; i++)
			cout<<lapSubDiagonal[i]<<" ";

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
