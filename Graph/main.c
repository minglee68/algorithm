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
int findDecendants(struct DFS_vertex* DFS_vertices, int num_v, int rootIndex, int*decendants);
int findSSCs(struct DFS_vertex* DFS_vertices, int num_v, int**SSCs);
void printSSCs(int** SSCs, int num_SSCs, char* v_names);

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
	quickSort(finishTimes, 0, num_v - 1, originalIndex);
	// try heapSort() instead of quickSort
	// try countingSort() instead of quickSort
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
