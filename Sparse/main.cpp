#include "os.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include "linalg.h"

#include <vector>

#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
extern "C" {
#include "trlan.h"
#include "trl_map.h"
#include "trl_comm_i.h"
}

#ifdef WINDOWS_OS
#include <windows.h>
#endif

using namespace std;
using namespace alglib;

#define TARGET_NUMEIGEN 300

vector<vector<pair<int, double> > > matrixRows;
int numVertex, numFaces, nPoints; //the numbers of vertices and faces

#ifdef TRL_FOTRAN_COMPATIBLE
#warning Need to rewrite for different fortran abi!
#endif

/* Matrix-vector multiplication
 *
 * The matrix is stored in a global
 * Multiple vectors are multiplied in one call -- ncol × nrow is a block of transposed vectors
 */
void diag_op(const int nrow, const int ncol, const double *xin, const int ldx,
	     double *yout, const int ldy, void* mvparam) {
    int i, j;
    /* ncol: number of vectors */
    /* nrow: number of rows in vector block */
    /* xin: vectors */
    /* ldx: vector stride */
    /* yout: multiplication result */
    /* ldy: multiplication result stride */

    for( j=0; j<ncol; j++ ) {
	for( i=0; i<nrow; i++ ) {
          double r = 0.0;
          vector<pair<int, double> >::iterator rowEnd = matrixRows[i].end();
          for (vector<pair<int, double> >::iterator mi = matrixRows[i].begin();
               mi != rowEnd;
               ++mi) {
            r += xin[j*ldx+mi->first] * mi->second;
          }
          yout[j*ldy+i] = r;
	}
    }
}

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
        //
	int **laplacian, *diagonal; //used for reading the laplacian
	int *faceP; //used to read points in a face
	bool bRelation;

	char fileName[100];
	char fullPath[150];
        const char *outPath = NULL;
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
          if (argc >= 3)
            outPath = argv[2];
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

        // Analyze sparseness of the matrix
        int zeroes = 0;
        int values = 0;

        matrixRows.resize(numVertex);
        // numVertex rows initialized to empty => all-zeroes

	for(int i=0; i<numVertex; i++)
          for(int j=0; j<numVertex; j++)
            if(i==j) {
              // matrix[i][j] = diagonal[i];
              matrixRows[i].push_back(pair<int,double>(j,diagonal[i]));
              values ++;
            }
            else {
              bRelation=false;
              for(int k=0; k<diagonal[i]; k++)
                if(laplacian[i][k]==j)
                  bRelation=true;
              if(bRelation) {
                // matrix[i][j] = -1;
                matrixRows[i].push_back(pair<int,double>(j, -1));
                values++;
              }
              else {
                // matrix[i][j] = 0;
                zeroes++;
              }
            }

        cerr << "Sparseness: " << endl;
        cerr << "Zeroes:     " << zeroes << endl;
        cerr << "Non-zeroes: " << values << endl;

        // forget zeroes and values

	//use the function to compute the eigenvalues
	//*********************************************//

        /* maxlan:
           maximum lanczos base size
           rule of thumb: >= ned + min(6, ned)
           large numbers increase performance
           too large numbers make orthogonalization expensive
           ned:
           desired number of eigenvalues
        */
        double *eval, *evec;
        int numEigen;

        { /* TODO: refactor into subroutine */
          static const int nrow=numVertex, lohi=1, ned=TARGET_NUMEIGEN, maxlan=400, mev=ned;
          int lwrk;
          // local variable declaration
          //double eval[mev], evec[mev*nrow], exact[mev];
          eval = new double[mev];
          evec = new double[mev * nrow];
          double *exact = new double[mev];

          double *res, *wrk;
          trl_info info;
          int i, j, k, fp, check;
          char name2[133], name[150];
          int tmp1, tmp2, nlen;
          // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
          // of the diag_op
          lwrk=maxlan*(maxlan+10);
          if( lwrk > 0 ) {
            res = (double*)malloc(lwrk*sizeof(double));
            wrk = (double*)malloc(lwrk*sizeof(double));
          }
          trl_init_info( &info, nrow, maxlan, lohi, ned, 1.4901e-8, 1, 2000000, -1 );
          trl_set_iguess( &info, 0, 1, 0, NULL );
          // the Lanczos recurrence is set to start with [1,1,...,1]^T
          memset(eval, 0, mev*sizeof(double) );
          for( i=0; i<nrow; i++ ) {
            evec[i] = 1.0;
          }
          info.verbose = -1; // Be somewhat quiet. Still outputs timing and a bit of other stuff
          //info.verbose =  8; // no need, really
          // call TRLAN to compute the eigenvalues
          trlan(diag_op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
          trl_print_info(&info, 3*nrow);
          for( i=0; i<mev; i++ ) { // need to match with the definition in diag_op
            exact[i] = (i+1)*(i+1);
          }
          if( info.nec > 0 ) {
            i = info.nec;
          } else {
            i = mev - 1;
          }
          if( info.verbose >= 0) {
            trl_check_ritz( diag_op, &info, nrow, i, evec, nrow, eval,
                            &check, res, exact, wrk, lwrk);
          } else {
            trl_terse_info( &info, 0 );
          }
          if( info.verbose > 1 ) {
            trl_rayleigh_quotients( diag_op, &info, i, evec, nrow, res, NULL );
            trl_check_ritz( diag_op, &info, nrow, i, evec, nrow, res, &check,
                            &res[i], exact, wrk, lwrk );
          }
          if( info.nec == 0 ) info.nec = min(info.ned, mev-1);
          if( lwrk > 0 ) {
            free(res);
            free(wrk);
          }
          numEigen = info.nec;
        }

        // @post numEigen <= TARGET_NUMEIGEN

        // Eigen computation done

        // @pre i in [0 .. numEigen)
        // @pre j in [0 .. numVertex)
#define EVector(i,j) evec[i * numVertex + j]

        // @pre x in [0 .. numEigen)
#define EValue(x) eval[x]

        assert(numEigen > 0);
        cerr << "Number of eigenValues determined: " << numEigen << endl;

	//HeatKernelSignature variables
        real_1d_array w; //eigenvalues vector
	double eSum=0; //will hold e^(-t*eigenValue)
	double tmin, tmax, step, t; //time variables
	double HKS[100];

	for(int i=0; i<100; i++)
		HKS[i]=0;

	w.setlength(numVertex);
	//for(int i=0; i<numVertex; i++)
	//	eigenValues[i]=0;

	//eigenVector.setlength(numVertex,numVertex);


	cout<<"\n";

	/*smatrixevd(matrix, numVertex, 1, false, eigenValues, eigenVector);

	for(int i=0; i<numVertex; i++)
        cout<<eigenValues[i]<<"\n"; */


	//compute the eigenvectors and calculate the coresponding heat kernel


	//initialize time variables
	tmin=4*log(10.0)/fabs(EValue(numEigen-1));
	tmax=4*log(10.0)/fabs(EValue(0));
	step=(log(tmax)-log(tmin))/100;

	//compute the HKS for 100 t's
	for(int k=0; k<100; k++)
	{
		//compute the t
		t=tmin+exp(k*step);

		for(int i=0; i<numEigen; i++)
			w[i]=0;
#ifdef USE_ESUM
		eSum=0;
		//e^(-t*eigenValue)
                for (int i = 0; i < numEigen; i++)
                  eSum+=fabs(pow(2.71828,-EValue(i)*t));
#endif

		//in w we store each heatKernel for each vertex
                /*for (int j = 0; j < 100; j++)
                  for (int i = 0 ; i < numVertex; i++)
                  w[i]+=EVector(j, i)*EVector(j, i)*pow(2.71828,-EValue(j)*t);*/
                for (int i = 0; i < numEigen; i++)
                  for (int j = 0; j < numVertex; j++)
                    w[i] += exp(-EValue(i)*t) * EVector(i,j)*EVector(i,j);

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
        if (outPath) {
          strcpy (fullPath, outPath);
          strcpy (fullPath, "/");
        }
        strcpy (fullPath, fileName);
        strcat (fullPath, ".signature");
#ifdef GIT_COMMIT
        strcat (fullPath, ".");
        strcat (fullPath, GIT_COMMIT);
#endif
#endif
	writeFile = fopen (fullPath, "w"); //open the file for writing
	if(writeFile==NULL)
			printf("Error opening the file");

	for(int i=0; i<100; i++)
		fprintf(writeFile,"%g\n", HKS[i]);

	fclose(writeFile);

#ifdef WINDOWS_OS
        system("pause");
#endif
	fcloseall();
	return 0;
}
