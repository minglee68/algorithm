#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

#define INFINITY INT_MAX

struct adj_list
{
	int v_index;
	int w;
	struct adj_list* next;
};

const int max_num_lists = 512; //  We expect no more than 512 lists
struct adj_list* all_lists;
int num_all_lists = 0; // initialization

enum v_color {
	WHITE, // 0
	GRAY, // 1
	BLACK // 2
};

struct DFS_vertex
{
	int dis; // discovery time
	int fin; // finish time
	int color;
	int pi; // index of parent vertex
};

int getVeticesNames(char** v_names, char* filename);
void getAdjMat(int** Adj_mat, int num_v, char* filename);
void printMat(int** Adj_mat, int num_v);
void getAdjArray(struct adj_list** Adj_array, int num_v, int** Adj_mat);
void printAdjArray(struct adj_list** Adj_array, int num_v, char* v_names);
int DFS_VISIT(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int u, int time);
void DFS(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int source_index, int* indexOrder);
void printDFS(struct DFS_vertex* DFS_vertices, int num_v, char* v_names);
void getTransAdjMat(int** Adj_mat, int** Adj_mat_trans, int num_v);
void swapItems(int* Array, int i, int j);
int partition(int* Array, int p, int r, int* originalIndex);
void quickSort(int* Array, int p, int r, int* originalIndex);
void DownAdjust(int * Array, int * originalIndex, int mid, int arraySize);
void heapSort(int* Array, int * originalIndex, int arraySize);
void countingSort(int * Array, int * originalIndex, int arraySize);
int findDecendants(struct DFS_vertex* DFS_vertices, int num_v, int rootIndex, int*decendants);
int findSSCs(struct DFS_vertex* DFS_vertices, int num_v, int**SSCs);
void printSSCs(int** SSCs, int num_SSCs, char* v_names);
void init_SPs(int**SPs, int**SP_pi, int arraySize, int s);
void Relax(int**SPs, int**SP_pi, int u, int v, int w, int s);
void BellmanFord(int**SPs, int**SP_pi, struct adj_list** Adj_array, int arraySize);
int isNotEmpty(int Q[], int arraySize);
int ExtractMin(int**SPs, int Q[], int arraySize, int s);
void Dijkstra(int**SPs, int**SP_pi, struct adj_list** Adj_array, int arraySize);
void printSPs(int**SPs, int arraySize, char* v_names);

int main()
{
	// input data
	char filename[] = "HW4.dat";

	// getting the number and names of vertices
	char* v_names = NULL; // assume each vertex has a single-charactor name
	int num_v = getVeticesNames(&v_names, filename); // num_v: number of all vertices in the given graph

	// getting adjacency matrix
	int** Adj_mat = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) Adj_mat[i] = malloc(sizeof(int)*num_v);
	getAdjMat(Adj_mat, num_v, filename);
	printf("Adjacency Matrix\n");
	printMat(Adj_mat, num_v);

	printf("\n");

	// getting the array of adjacency list from the adjacency matrix
	struct adj_list** Adj_array;
	Adj_array = malloc(sizeof(struct adj_list*)*num_v);
	all_lists = malloc(sizeof(struct adj_list)*max_num_lists);
	getAdjArray(Adj_array, num_v, Adj_mat);
	printAdjArray(Adj_array, num_v, v_names);

	printf("\n");


	// Find Single Source Shortest-Path with BellmanFord
	int** SPs = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) SPs[i] = malloc(sizeof(int)*num_v);
	int** SP_pi = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) SP_pi[i] = malloc(sizeof(int)*num_v);
	BellmanFord(SPs, SP_pi, Adj_array, num_v);
	printf("All Pairs Shortest-Path with Bellman Ford Algorithm\n");
	printSPs(SPs, num_v, v_names);
	printf("\n");

	// Find Single Source Shortest-Path with Dijkstra
	Dijkstra(SPs, SP_pi, Adj_array, num_v);
	printf("All Pairs Shortest-Path with Dijkstra Algorithm\n");
	printSPs(SPs, num_v, v_names);
	printf("\n");


	// DFS
	struct DFS_vertex* DFS_vertices;
	DFS_vertices = malloc(sizeof(struct DFS_vertex)*num_v);
	int source_index = 0;
	printf("With %c as the source vertex,\n", v_names[source_index]);
	DFS(Adj_array, num_v, DFS_vertices, source_index, NULL);
	printDFS(DFS_vertices, num_v, v_names);

	printf("\n");

	// Transpose of adjacency matrix
	int** Adj_mat_trans = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) Adj_mat_trans[i] = malloc(sizeof(int)*num_v);
	getTransAdjMat(Adj_mat, Adj_mat_trans, num_v);
	printf("Transpose of Adjacency Matrix\n");
	printMat(Adj_mat_trans, num_v);

	printf("\n");

	// getting the array of adjacency list from the transpose of the adjacency matrix
	struct adj_list** Adj_array_trans;
	Adj_array_trans = malloc(sizeof(struct adj_list*)*num_v);
	getAdjArray(Adj_array_trans, num_v, Adj_mat_trans);
	printAdjArray(Adj_array_trans, num_v, v_names);

	printf("\n");

	// sorting the indices of vertices according to the finishi time
	int* finishTimes = malloc(sizeof(int)*num_v);
	int* originalIndex = malloc(sizeof(int)*num_v);
	for (int i = 0; i<num_v; i++)
	{
		finishTimes[i] = DFS_vertices[i].fin;
		originalIndex[i] = i;
	}
	//quickSort(finishTimes, 0, num_v - 1, originalIndex);
	//heapSort(finishTimes, originalIndex, num_v);
	countingSort(finishTimes, originalIndex, num_v);
	for (int i = 0; i<(int)((double)num_v / 2.0); i++)
	{
		swapItems(originalIndex, i, num_v - i - 1); // to get the reversed order
	}

	// Second DFS
	struct DFS_vertex* DFS_vertices_second;
	DFS_vertices_second = malloc(sizeof(struct DFS_vertex)*num_v);
	source_index = originalIndex[0];
	printf("With %c as the source vertex,\n", v_names[source_index]);
	DFS(Adj_array_trans, num_v, DFS_vertices_second, source_index, originalIndex);
	printDFS(DFS_vertices_second, num_v, v_names);

	printf("\n");
	printf("\n");

	// Find SSCs
	int** SSCs = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) SSCs[i] = malloc(sizeof(int)*num_v);
	int num_SSCs = findSSCs(DFS_vertices_second, num_v, SSCs);
	printSSCs(SSCs, num_SSCs, v_names);


	// freeing allocated memories
	for (int i = 0; i<num_v; i++)
	{
		free(Adj_mat[i]);
		free(Adj_mat_trans[i]);
		free(SSCs[i]);
	}
	free(Adj_mat);
	free(all_lists);
	free(Adj_array);
	free(DFS_vertices);
	free(Adj_mat_trans);
	free(Adj_array_trans);
	free(finishTimes);
	free(originalIndex);
	free(DFS_vertices_second);
	free(SSCs);
	free(SPs);
	free(SP_pi);

	return 0;
}

int getVeticesNames(char** v_names, char* filename)
{
	FILE* fp;
	fp = fopen(filename, "r");

	int num_v = 0;

	char buffer[256];
	fgets(buffer, 256, fp); // first line
	for (int i = 0; i<strlen(buffer) - 1; i++)
	{
		if (buffer[i] != '\t')
		{
			num_v++;
			*v_names = realloc(*v_names, sizeof(char)*num_v);
			(*v_names)[num_v - 1] = buffer[i];
		}
	}

	fclose(fp);

	return num_v;
}

void getAdjMat(int** Adj_mat, int num_v, char* filename)
{
	FILE* fp;
	fp = fopen(filename, "r");

	char buffer[256];
	fgets(buffer, 256, fp); // skip the first line

	for (int i = 0; i<num_v; i++)
	{
		for (int j = -1; j<num_v; j++)
		{
			fscanf(fp, "%s", buffer);
			if ((j != -1) && (strcmp(buffer, "INF") != 0))
				Adj_mat[i][j] = atoi(buffer);
			else if (strcmp(buffer, "INF") == 0)
				Adj_mat[i][j] = INFINITY; // just setting the maximum value of integer
		}
	}

	fclose(fp);
	return;
}

void printMat(int** Adj_mat, int num_v)
{
	for (int i = 0; i<num_v; i++)
	{
		for (int j = 0; j<num_v; j++)
		{
			if (Adj_mat[i][j] == INFINITY)
				printf("INF\t");
			else
				printf("%d\t", Adj_mat[i][j]);
		}
		printf("\n");
	}
}

void getAdjArray(struct adj_list** Adj_array, int num_v, int** Adj_mat)
{
	for (int i = 0; i<num_v; i++)
	{
		num_all_lists++;
		if (num_all_lists > max_num_lists) exit(1);
		all_lists[num_all_lists - 1].v_index = i;
		all_lists[num_all_lists - 1].w = 0;
		all_lists[num_all_lists - 1].next = NULL;
		Adj_array[i] = &all_lists[num_all_lists - 1];

		for (int j = 0; j<num_v; j++)
		{
			if (i != j)
			{
				if (Adj_mat[i][j] != INFINITY)
				{
					num_all_lists++;
					if (num_all_lists > max_num_lists) exit(1);
					all_lists[num_all_lists - 1].v_index = j;
					all_lists[num_all_lists - 1].w = Adj_mat[i][j];
					all_lists[num_all_lists - 1].next = NULL;
					all_lists[num_all_lists - 2].next = &all_lists[num_all_lists - 1];
				}
			}
		}
	}

	return;
}

void printAdjArray(struct adj_list** Adj_array, int num_v, char* v_names)
{
	// print the array of adjacency list (Adj_array)
	printf("The array of adjacency list\n");

	for (int i = 0; i<num_v; i++)
	{
		printf("%c: ", v_names[i]);
		struct adj_list* cur_list;
		cur_list = Adj_array[i]->next;
		while (cur_list != NULL)
		{
			printf("%c, %d; ", v_names[cur_list->v_index], cur_list->w);
			cur_list = cur_list->next;
		}
		printf("\n");
	}
}

int DFS_VISIT(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int u, int time)
{
	// change color to gray
	DFS_vertices[u].color = GRAY;
	time++;
	// update discovery time
	DFS_vertices[u].dis = time;
	struct adj_list* cur_list = Adj_array[u]->next;
	while (cur_list != NULL)
	{
		if (DFS_vertices[cur_list->v_index].color == WHITE)
		{
			DFS_vertices[cur_list->v_index].pi = u;
			time = DFS_VISIT(Adj_array, num_v, DFS_vertices, cur_list->v_index, time);
		}
		cur_list = cur_list->next;
	}
	// change color to black
	DFS_vertices[u].color = BLACK;
	time++;
	// update finish time
	DFS_vertices[u].fin = time;

	return time;
}

void DFS(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int source_index, int* indexOrder)
{
	for (int u = 0; u<num_v; u++)
	{
		//each vertex's color has to be initialized to WHITE;
		DFS_vertices[u].color = WHITE;
		DFS_vertices[u].pi = -1; //NULL;
	}
	int time = 0;

	// do DFS_VISIT for the source vertex first
	time = DFS_VISIT(Adj_array, num_v, DFS_vertices, source_index, time);

	for (int i = 0; i<num_v; i++)
	{
		if (indexOrder == NULL)
		{
			//if i-th vertex' color is white
			if (DFS_vertices[i].color == WHITE) {
				time = DFS_VISIT(Adj_array, num_v, DFS_vertices, i, time);
			}
		}
		else
		{
			//if (indexOrder[i])-th vertex' color is white
			if (DFS_vertices[indexOrder[i]].color == WHITE) {
				time = DFS_VISIT(Adj_array, num_v, DFS_vertices, indexOrder[i], time);
			}
		}

	}
}

void printDFS(struct DFS_vertex* DFS_vertices, int num_v, char* v_names)
{
	printf("DFS result\n");
	for (int u = 0; u<num_v; u++)
	{
		if (DFS_vertices[u].pi == -1)
			printf("%c: d=%d, f=%d, pi=%d(root)\n", v_names[u], DFS_vertices[u].dis, DFS_vertices[u].fin, DFS_vertices[u].pi);
		else
			printf("%c: d=%d, f=%d, pi=%d(%c)\n", v_names[u], DFS_vertices[u].dis, DFS_vertices[u].fin, DFS_vertices[u].pi, v_names[DFS_vertices[u].pi]);
	}
}

void getTransAdjMat(int** Adj_mat, int** Adj_mat_trans, int num_v)
{
	for (int i = 0; i<num_v; i++)
	{
		for (int j = 0; j<num_v; j++)
		{
			Adj_mat_trans[i][j] = Adj_mat[j][i];
		}
	}

	return;
}

void swapItems(int* Array, int i, int j)
{
	int temp = Array[i];
	Array[i] = Array[j];
	Array[j] = temp;
}

int partition(int* Array, int p, int r, int* originalIndex)
{
	int x = Array[r];
	int i = p - 1;
	for (int j = p; j<r; j++)
	{
		if (Array[j] <= x)
		{
			i++;
			swapItems(Array, i, j);
			swapItems(originalIndex, i, j);
		}
	}
	swapItems(Array, i + 1, r);
	swapItems(originalIndex, i + 1, r);

	return i + 1;
}

void quickSort(int* Array, int p, int r, int* originalIndex)
{
	if (p < r)
	{
		int q = partition(Array, p, r, originalIndex);
		quickSort(Array, p, q - 1, originalIndex);
		quickSort(Array, q + 1, r, originalIndex);
	}
}

void DownAdjust(int * Array, int * originalIndex, int mid, int arraySize) {
	int parent = mid;
	int left = parent * 2 + 1;
	int right = parent * 2 + 2;
	int largest = parent;

	if (left < arraySize && Array[left] > Array[largest]) {
		largest = left;
	}
	if (right < arraySize && Array[right] > Array[largest]) {
		largest = right;
	}
	if (parent != largest) {
		swapItems(Array, parent, largest);
		swapItems(originalIndex, parent, largest);
		DownAdjust(Array, originalIndex, largest, arraySize);
	}
}

void heapSort(int * Array, int * originalIndex, int arraySize) {
	int i;

	for (i = (arraySize / 2) - 1; i >= 0; i--) {
		DownAdjust(Array, originalIndex, i, arraySize);
	}

	while (arraySize > 1){
		swapItems(Array, 0, arraySize-1);
		swapItems(originalIndex, 0, arraySize-1);
		arraySize--;
		DownAdjust(Array, originalIndex, 0, arraySize);
	}
}

void countingSort(int * Array, int * originalIndex, int arraySize) {
	int* work_data = malloc(sizeof(int) * arraySize);
	int* work_data_orig = malloc(sizeof(int) * arraySize);
	int* count = malloc(sizeof(int) * (arraySize * 2 + 2));

	int i;

	for (i = 0; i <= arraySize * 2; i++)
		count[i] = 0;

	for (i = 0; i < arraySize; i++)
		count[Array[i]]++;

	for (i = 0; i < arraySize * 2; i++)
		count[i + 1] += count[i];

	for (i = arraySize - 1; i >= 0; i--) {
		work_data[--count[Array[i]]] = Array[i];
		work_data_orig[count[Array[i]]] = originalIndex[i];
	}

	for (i = 0; i < arraySize; i++) {
		Array[i] = work_data[i];
		originalIndex[i] = work_data_orig[i];
	}

	free(count);
	free(work_data);
	free(work_data_orig);
}

int findDecendants(struct DFS_vertex* DFS_vertices, int num_v, int rootIndex, int*decendants)
{
	int num_decendants = 0;

	for (int u = 0; u<num_v; u++)
	{
		if ((DFS_vertices[u].dis > DFS_vertices[rootIndex].dis) && (DFS_vertices[u].fin < DFS_vertices[rootIndex].fin))
		{
			num_decendants++;
			decendants[num_decendants - 1] = u;
		}
	}

	return num_decendants;
}

int findSSCs(struct DFS_vertex* DFS_vertices, int num_v, int**SSCs)
{
	// assuming SSCs has (num_v x num_v) dimension (maximum dimension)
	// use findDecendants function
	// how to indicate there is no further element? (hint: look at printSSCs() function)

	int num_SSCs = 0;
	int* used_v = malloc(sizeof(int)*num_v);
	for (int u = 0; u < num_v; u++) used_v[u] = 0;

	for (int u = 0; u < num_v; u++) {
		int i = 0;
		int* decendants = malloc(sizeof(int)*num_v);
		int num_decendants = findDecendants(DFS_vertices, num_v, u, decendants);

		if (used_v[u] == 0) {
			SSCs[num_SSCs][0] = u;
			used_v[u] = 1;

			for (i = 1; i < num_decendants + 1; i++) {
				SSCs[num_SSCs][i] = decendants[i-1];
				used_v[decendants[i-1]] = 1;
			}

			SSCs[num_SSCs][i] = -1;
			num_SSCs++;
		}
		free(decendants);
	}

	free(used_v);
	return num_SSCs;
}

void printSSCs(int** SSCs, int num_SSCs, char* v_names)
{
	printf("SSCs\n");
	for (int i = 0; i<num_SSCs; i++)
	{
		int index = 0;
		printf("SSC%d: ", i + 1);
		while (SSCs[i][index] >= 0) // if it meats negative value, it means there is no further element.
		{
			printf("%c, ", v_names[SSCs[i][index]]);
			index++;
		}
		printf("\n");
	}
}

void init_SPs(int**SPs, int**SP_pi, int arraySize, int s) {
	int j;
	for (j = 0; j < arraySize; j++) {
		if (s == j) {
			SPs[s][j] = 0;
			SP_pi[s][j] = -1;
		}
		else {
			SPs[s][j] = INFINITY;
			SP_pi[s][j] = -1;
		}
	}
}

void Relax(int**SPs, int**SP_pi, int u, int v, int w, int s)  {
	//printf("%d > %d + %d\n", SPs[s][v], SPs[s][u], w);
	if (SPs[s][u] != INFINITY) {
		if (SPs[s][v] > SPs[s][u] + w) {
			SPs[s][v] = SPs[s][u] + w;
			SP_pi[s][v] = u;
		}
	}
}

void BellmanFord(int**SPs, int**SP_pi, struct adj_list** Adj_array, int arraySize) {
	int i, j, s;
	for (s = 0; s < arraySize; s++) {
		init_SPs(SPs, SP_pi, arraySize, s);
		for (i = 0; i < arraySize - 1; i++) {
			for (j = 0; j < arraySize; j++) {
				//printf("%d\n", j);
				struct adj_list* cur_list = Adj_array[j]->next;
				while (cur_list != NULL) {
					//printf("%d: %d\n", cur_list->v_index, cur_list->w);
					Relax(SPs, SP_pi, j, cur_list->v_index, cur_list->w, s);
					cur_list = cur_list->next;
				}
			}
			//printf("next\n");
		}
	}
}

int isNotEmpty(int Q[], int arraySize) {
	int flag = 0;
	int i;
	for (i = 0; i < arraySize; i++) {
		if (Q[i] == 1)
			flag = 1;
	}
	return flag;
}

int ExtractMin(int**SPs, int Q[], int arraySize, int s) {
	int min = INFINITY;
	int min_index = -1;
	int i;
	for (i = 0; i < arraySize; i++) {
		if (Q[i] == 1 && min > SPs[s][i]){ 
			min = SPs[s][i];
			min_index = i;
		}
	}

	return min_index;
}

void Dijkstra(int**SPs, int**SP_pi, struct adj_list** Adj_array, int arraySize) {
	int i, j, s;
	for (s = 0; s < arraySize; s++) {
		init_SPs(SPs, SP_pi, arraySize, s);
		int S[7] = {0};
		int Q[7] = {1, 1, 1, 1, 1, 1, 1};
		//printf("s = %d\n", s);
		while (isNotEmpty(Q, arraySize) == 1) {
			int u = ExtractMin(SPs, Q, arraySize, s);
			if (u == -1) break;
			//printf("u: %d\n", u);
			Q[u] = 0;
			S[u] = 1;
			//printf("Q: ");
			//for(j = 0; j < arraySize; j++) {
			//	printf("%d, ", Q[j]);
			//}
			//printf("\n");
			struct adj_list* cur_list = Adj_array[u]->next;
			while (cur_list != NULL) {
				//printf("%d: %d\n", cur_list->v_index, cur_list->w);
				Relax(SPs, SP_pi, u, cur_list->v_index, cur_list->w, s);
				cur_list = cur_list->next;
			}
		}
		//printf("\n");
	}
}

void printSPs(int**SPs, int arraySize, char* v_names) {
	int i, j;
	printf("\t");
	for (i = 0; i < arraySize; i++) {
		printf("%c\t", v_names[i]);
	}
	printf("\n");
	for (i = 0; i < arraySize; i++) {
		printf("%c\t", v_names[i]);
		for (j = 0; j < arraySize; j++) {
			if (SPs[i][j] == INFINITY)
				printf("INF\t");
			else
				printf("%d\t", SPs[i][j]);
		}
		printf("\n");
	}
}
