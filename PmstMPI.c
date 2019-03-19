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

//Structs
/*
*Tuple structure that represents the minimum weight found by a processor
*with the id specidfied with the attribute processorId.
*/
typedef struct 
{
	int minimumWeight;
	int processorId;
	
}Tuple;
/*
*Edge structure that represents origin vertex and destination vertex
*/
typedef struct 
{
	int origin;
	int destination;
	
} Edge;

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
				randomValue = 99999;
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
	int *matrixSegment;			//Segment of the matrix that is assigned to each Processor
	int *minumumSpanningTree;	//Result
	int actualMinWeight;		

	double startTime;
	double endTime;
	double totalTime;

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
	*4) Master/Root process broadcast the number of vertices in the matrix to
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

	//Assign the corresponding segment of data to each processor: the number of rows of the matrix assigned to each processor
	//NumberVerticesAssigned*NumVerticesToWhichIsConnectedAvertex*sizeInt = TotalIntegers For TheRowsOfNumberOfVertexAssigned to each processor
	matrixSegment = (int*)malloc(numVerticesToSend[resourceProcessId]*VERTICES*sizeof(int));

	//Define an MPI_Datatype
	MPI_Datatype continuousMatrix;

	/*
	*5) Create a contiguous datatype that represents a Row of integers
	*
	*INPUT:
	*1)count = replication count -> VERTICES {int} 
	*2)oldtype = old data type -> MPI_INT {*MPI_Datatype}
	*
	*OUTPUT:
	*3)newtype = new data type -> continuousMatrix {*MPI_Datatype}
	*/
	MPI_Type_contiguous(VERTICES,MPI_INT,&continuousMatrix);

	/*
	*6) Commit the data type
	*INPUT: 
	*1)datatype -> continuousMatrix {MPU_Datatype}
	*/
	MPI_Type_commit(&continuousMatrix);

	/*
	*7) Send to each processor the chunk (Rows) that they require from all the 
	*matrix content.
	*
	*INPUT
	*1)sendbuf = address of send buffer -> graph {memory address}
	*2)sendcounts = sepecify the number of elements to send to each processor ->numVerticesToSend {int array}
	*3)displs = Sepecifies the displacemnt/offset from wich start to send data to process i -> displacement {int array}
	*4)sendtype = Data type sent in buffer -> continuousMatrix {MPI_Datatype}
	*6)recvcount = number of elements recived in recvbuf -> numVerticesToSend[resourceProcessId] {int}
	*7)recvtype = Data type recived in buffer -> continuousMatrix {MPI_Datatype}
	*8)root = ID of sending process -> MASTER {int}
	*9)comm = Communicator throught wich the communication is performed -> MPI_COMM_WORLD {MPI CONSTANT}
	*
	*OUTPUT
	*5)recvbuf = Address to recieve the sended buffer -> matrixSegment {memory address}
	*/
	MPI_Scatterv(graph,numVerticesToSend,displacement,continuousMatrix,matrixSegment,numVerticesToSend[resourceProcessId],continuousMatrix,MASTER,MPI_COMM_WORLD);

	//----------
	//Test that chunks were successfully send to a specific processor
	/*if(resourceProcessId == MASTER)
	{
		for(int i = 0; i < VERTICES; i++)
		{
			for(int j = 0; j < VERTICES; j++)
			{
				printf("%i ", matrixSegment[VERTICES*i+j] );

			}//End for 2
			printf("\n");	

		}//End for 1	

	}//End if*/	
	//----------

	//Initialize array to save results
	minumumSpanningTree = (int*)malloc(sizeof(int)*VERTICES);

	for(int i = 0; i < VERTICES; i++)
	{
		minumumSpanningTree[i] = -1;

	}//End for initialize mst	

	//------- Start of MST Algorithm ------
	startTime = MPI_Wtime();

	//Set root of spanning tree
	minumumSpanningTree[0] = 0;
	int nextClosestVertext = 0;
	int actualVertex = 0;		//Vertex # from wich nextClosetVertext was reached

	Tuple tupleSend;
	Tuple tupleRecieved;
	Edge edge;
	
	//1)For each Vertex in the graph that has not been visited yet
	//(TotalVertices-RootVertex)
	for(int k = 0; k < VERTICES -1; k++)
	{
		//For each new vertex start by seting minimum wight equal to the largest
		//possible weigh in order to not consider not possible paths between vertices
		actualMinWeight = 99999;	

		//2)For each vertex assigned to the actual processor
		for(int i = 0; i < numVerticesToSend[resourceProcessId]; i++)
		{
			//---------
			printf("i: %d < numVerticesSend:%d\n", i,numVerticesToSend[resourceProcessId]);

			for(int i = 0; i < VERTICES; i++)
			{
				printf("%i ", minumumSpanningTree[i] );
					
			}//End for 1

			printf("\n");

			printf("¿ mst:%i != -1 ?\n",minumumSpanningTree[ i+displacement[resourceProcessId] ] );
			//---------

			//3)If the next vertex for the actual vertex i assigned to this processor
			//has been visited
			if( minumumSpanningTree[ i+displacement[resourceProcessId] ] != -1)
			{
				//4)For each vertex that has not been visited yet and that is next to
				//the vertex of 3)
				for(int j = 0; j < VERTICES; j++)
				{
					//Verify that has not been visited the vertex j
					if(minumumSpanningTree[j] == -1)
					{
						//-------
						printf("(i:%i,j:%i) %i < %i and  %i != 0 \n", i,j,matrixSegment[VERTICES*i+j],actualMinWeight,matrixSegment[VERTICES*i+j] );
						//------
						//5)Verify if the distance to vertex j is less than the
						//actualMinWeightt and that is a valid edge (Not 0)
						if( (matrixSegment[VERTICES*i+j] < actualMinWeight) && (matrixSegment[VERTICES*i+j] != 0) )
						{
							//New minimum/closer vertex found
							//1)Save new minimumWeight
							actualMinWeight = matrixSegment[VERTICES*i+j];
							//2)Save new closest vertex number
							nextClosestVertext = j;
							//3)Save the vertex from which 2) was reached
							actualVertex = i;

						}//End if 3	

					}//End if 2

				}//End for 5	

			}//End if 1

		}//End for 4	

		//Each processor has to send it is local solution to other processors
		tupleSend.minimumWeight = actualMinWeight;
		tupleSend.processorId = resourceProcessId;

		/*
		*8)Combines values from all processes and distributes the minimum value 
		*of all the elements back to all processers
		*
		*INPUT
		*1)sendbuf = address of send buffer -> &tupleSend {memory address}
		*3)count = # of elements in send buffer -> 1 {int}
		*4)datatype = -> MPI_2INT {MPI Dataype}
		*5)op = operation of reduction -> MPI_MINLOC {MPI constant}
		*6)comm = Communicator throught wich the communication is performed -> MPI_COMM_WORLD {MPI CONSTANT}
		*
		*
		*OUTPUT
		*2)recvbuf = address of receive buffer -> &tupleRecieved {memory address}
		*/
		MPI_Allreduce(&tupleSend, &tupleRecieved,1,MPI_2INT,MPI_MINLOC,MPI_COMM_WORLD);

		//Each processor has to send the min row to others processors
		edge.origin = actualVertex + displacement[resourceProcessId];
		edge.destination = nextClosestVertext;

		/*
		*9)Sends a message from root process to all the processes in the same group.
		*
		*INPUT/OUTPUT
		*1)buf = address of buffer -> &edge {memory address}
		*
		*INPUT
		*2)count = # of elements in send buffer -> 1 {int}
		*3)datatype = -> MPI_2INT {MPI Dataype}
		*4)root = ID of sending process -> tupleRecieved.processorId {int}
		*5)comm = Communicator throught wich the communication is performed -> MPI_COMM_WORLD {MPI CONSTANT}
		*/
		MPI_Bcast(&edge,1,MPI_2INT,tupleRecieved.processorId,MPI_COMM_WORLD);

		//Update the minimum spanning tree with the minimum vertex decided
		//by all the processors with the MINLOC operation
		minumumSpanningTree[edge.destination] = edge.origin; 


	}//End for 3

	//---------- END MSPT ALGORITHM ---------

	endTime = MPI_Wtime();
	totalTime = endTime - startTime;

	//------------ END OF WORK FOR ALL THE PROCESSORS ----------

	//------------ WORK EXCLUSIVE OF THE MASTER -----------------

	//Test result
	//if(resourceProcessId == MASTER)
	//{
		/*for(int i = 0; i < VERTICES; i++)
		{
			printf("%i ", minumumSpanningTree[i] );
				
		}//End for 1

		printf("\n");*/	

	//}//End if
	//----------

	//------------ END OF WORK EXCLUSIVE OF THE MASTER ----------

	//++++++++++ OUTPUT DATA ++++++++++++


	//10)Ends the prallel communication between processors(After this processors
	//can not continue sending messages)
	MPI_Finalize();

	//++++++++++ END OF PROGRAM ++++++++
	printf("End of program V1\n");

	//Indicate to the OS that all ends OK
	return 0;

}//End of main
