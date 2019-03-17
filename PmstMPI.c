/*
* Prim Minimum Spanning Tree implementation Parallel implementation using MPI.
*
* Based on the @elahehrashedi version.
*/

//Libraries
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//Constants
int VERTICES = 0;
int MASTER = 0; 	//PROCESS ID OF MASTER PROCESS

//+++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++
/*
*Generate a random matrix given the number of vertices. The matrix
*is rerpresented in form an array
*
*The matrix generated has the same values of the upper triangle
*for the botom part put in inverse form.
*/
//-----------------------------------------------------------
int* generateRandomGraph(int numVertices)
{
	int *randomGraph;
	int randomValue;
	
	VERTICES = numVertices;
	
	//Asignar espacio de memoria que ocupa toda la matriz en una sola fila. 
	//.:. Espacio a reservar para todas las fila = NumeroFilas*NumeroClumnas(SizeOfPointerInt*)
	randomGraph = (int *)malloc( (numVertices*numVertices)*sizeof(int*));
	
	//Generar numeros aleatorios del 0 al 20 y llenar con ellos todos los espacios en matriz
	for(int i = 0; i < numVertices; i++)
	{
		for(int j = 0; j < numVertices; j++)
		{
			//srand(time(NULL));	//Seed to generate random numbers
			randomValue = rand() % 20;
			
			//Los valores que den 0 se guardan como INFINITO
			if(randomValue == 0)
			{
				randomValue = 9999;
			}//Fin if
			
			//Hacer que los valores de los caminos sean los mismos
			//en ambas direcciones
			randomGraph[(i*numVertices)+j] = randomValue;
			//----
			//printf("%i) %i \n", (i*numVertices+j),randomValue);
			//----
			randomGraph[(j*numVertices)+i] = randomValue; //Transpuesta de matriz
			//----
			//printf("%i) %i \n", (j*numVertices+i),randomValue);
			//----
			
		}//Fin for 3
		
	}//Fin for 2
	
	//Asegurarse que la diagonal principal de la matriz es de puros 0s(Los vertices no pueden tener camino asi mismos)
	for(int l = 0; l < numVertices; l++)
	{
		randomGraph[l*numVertices+l] = 0;
		
	}//FIn for 4
	
	return randomGraph;
	
}//FIn metodo generateRandomGraph
//--------------------------------------------------------------
/*
* Print in form of matrix the grahp that is recieved as an array
*/
void printMatrix(int *matrix, int vertices)
{
	for(int i = 0; i < (vertices*vertices); i++)
	{
		
		if(i%vertices == 0)
		{
			printf("\n");

		}//End if	

		printf("%4i ", matrix[i]);

	}//End for	

	//End of matrix
	printf("\n");

}//End function printMatrix
//++++++++++++++++++++++ END FUNCTIONS ++++++++++++++++++++++++++++

//Strart of main function/program
int main(int argc, char *argv[])
{
	//Declaration of variables
	int *graph;
	int resourceProcessId; 		//ID current processor
	int totalProcessorsInComm; 	//Total number processors in communicator
	int *numVerticesToSend;		//# of vertices send to each Procesor
	int *displacement;			//Specifies the displacement of vertices for each Procesor
	int remainingVertices;		//#of vertices after distribute in equal parts the vertices according to the # of processors

	//++++ INITIALIZATION OF VARIABLES ++++
	graph = generateRandomGraph(5);
	printMatrix(graph,5);

	//+++++++++ INPUT DATA +++++++++++++++

	//+++++++++++++ PROCESS ++++++++++++

	/*
	*1) Initialize the structure of communication between Processors(Indispenable for communication)
	*/
	MPI_Init(&argc, &argv);

	/*
	*2) Save the rank/Id actual process in variable resourceProcessId
	*
	*Input: desired communicator that we want to know it is ID->MPI_COMM_WORLD
	*Output: int ID communicator -> save in resourceProcessId
	*/
	MPI_Comm_rank(MPI_COMM_WORLD, &resourceProcessId);

	/*
	*3) Determine the current number of processors associated with specific communicator
	*
	*Input: desired communicator that we want to know its size ->MPI_COMM_WORLD
	*Output: int #Processors in current communicator -> save in totalProcessorsInComm
	*/
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcessorsInComm);

	//-------
	printf("current process ID is: %i and the total # of processors is: %i \n", resourceProcessId, totalProcessorsInComm);
	//-------

	/*
	* Master/Root process broadcast the number of vertices in the matrix to
	* all the communicators
	* INPUT
	* 2)counter = #of elements to be sent-> 1 {INT}
	* 3)datatype = type of data sent -> MPI_INT {MPI DATA TYPE}
	* 4)root = id of root processor -> MASTER {int}
	* 5)comm = Communicator throught wich the communication is performed -> MPI_COMM_WORLD {MPI CONSTANT}
	*
	* OUTPUT/INPUT
	* 1)buffer = starting address of buffer -> VERTICES {*memory pointer}
	*/
	MPI_Bcast(&VERTICES, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	//------------ WORK FOR ALL THE PROCESSORS ------------------

	//Each procesor has a number of vertices assigned
	numVerticesToSend = (int*)malloc(sizeof(int)*VERTICES);

	//Each procesor has a different number of displacement
	displacement = (int*)malloc(sizeof(int)*VERTICES);

	//Determine the # that are assigned to each processor and the offset (first vertice)

	//Master process starts always from 0
	displacement[0] = 0;
	//Number of processors assigned to master depending of number of processors
	numVerticesToSend[0] = VERTICES/totalProcessorsInComm;
	//If number of vertices can not be devided proportionally divide the remainining vertices
	//to the processors. 
	remainingVertices = totalProcessorsInComm - (VERTICES % totalProcessorsInComm);

	//**Distribute Vertices that were divided proportionally succesfully.
	//Start from 1 becauset master has already set.
	//i < totalProcessorsInComm ?
	for(int i = 1; i < remainingVertices; i++)
	{
		//Actual processor has the same number of asigned vertices as the master
		numVerticesToSend[i] = numVerticesToSend[0];
		//Calculate offset of actual processor = offsetPreviousProcessor+#VerticesAssignedToPrviousProcessor
		displacement[i] = displacement[i-1] + numVerticesToSend[i-1];

	}//End for 1

	//**Distribute remaining vertices to processors
	for(int i = remainingVertices; i < totalProcessorsInComm; i++)
	{
		//Add remaining vertice to stadnar # of vertices of each processor
		numVerticesToSend[i] = numVerticesToSend[0] + 1;
		//Rcalculate the offset of this processor reconsidering the offset of prvious
		//procesor and the #of vertices that has assigned to it
		displacement[i] = displacement[i - 1] + numVerticesToSend[i-1];

	}//End for 2	

	//-------
	printf("PROCESID: %i # of asined vertices: %i with offset: %i \n", resourceProcessId, numVerticesToSend[resourceProcessId], displacement[resourceProcessId]);
	//-------



	//------------ END OF WORK FOR ALL THE PROCESSORS ----------

	//------------ WORK EXCLUSIVE OF THE MASTER -----------------

	//------------ END OF WORK EXCLUSIVE OF THE MASTER ----------

	//++++++++++ OUTPUT DATA ++++++++++++


	//Ends the prallel communication between processors(After this processors
	//can not continue sending messages)
	MPI_Finalize();

	//++++++++++ END OF PROGRAM ++++++++
	printf("End of program V1\n");

	//Indicate to the OS that all ends OK
	return 0;

}//End of main
