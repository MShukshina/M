Lab_1.cpp
#include "stdafx.h"
#include <ctime>
#include <iostream>
#define _CRTDBG_MAP_ALLOC 
#include <stdlib.h> 
#include <crtdbg.h> 

using namespace std;

struct TNode
{
	int Data;
	TNode* Right;
	TNode* Left;
};
typedef TNode* PNode;

PNode InsertInTree(int x, PNode Root)//add one item to a tree
{
	if (Root == NULL)
	{
		PNode q;
		q = new TNode;
		q->Data = x;
		q->Left = NULL;
		q->Right = NULL;
		return q;
	}
	if (x < Root->Data)
		Root->Left = InsertInTree(x, Root->Left);
	else if (x > Root->Data)
		Root->Right = InsertInTree(x, Root->Right);
	return Root;
}

void CreatTree(PNode &Root)//to create a tree
{
	int n = 5, x;
	for (int i = 0; i<n; i++)
	{
		x = rand()%51;
		Root = InsertInTree(x, Root);
	}
}

PNode FindInTree(int x, PNode Root)//find the item in a tree
{
	if (Root == NULL)
		return NULL;
	if (x < Root->Data)
		return FindInTree(x, Root->Left);
	if (x > Root->Data)
		return FindInTree(x, Root->Right);
	return Root;
}

void PrintTree( PNode Root, int k)//print a tree
{
	if (Root != NULL)
	{
		PrintTree(Root->Right, k + 1);
		for (int i = 1; i <= k; i++)
			cout << "      ";
		cout << Root->Data << endl;
		cout << endl;
		PrintTree(Root->Left, k + 1);
	}
}

void Delete(PNode &Root, PNode &q) //to remove a vertex in the tree-1
{
	if (Root -> Right != NULL)
		Delete(Root -> Right, q);
	else
	{
		q = Root;
		Root = Root -> Left;
	}
}

void DeleteElem(PNode &Root, int x) //to remove a vertex in a tree
{
	if (Root == NULL) return;
	if (x < Root -> Data) DeleteElem(Root -> Left, x);
	else if (x > Root->Data) DeleteElem(Root->Right, x);
	else 
	{
		PNode q = Root;
		if (Root -> Left == NULL)
			Root = Root -> Right;
		else if (Root -> Right == NULL)
			Root = Root -> Left;
		else
		{
			Delete(Root -> Left, q);
			Root -> Data = q -> Data;
		}
		delete q;
	}
}

void DeleteTree(PNode Root)//to clear the memory
{
	if (Root != NULL)
	{
		DeleteTree(Root->Left);
		DeleteTree(Root->Right);
		delete Root;
	}
}

int main()
{
	srand(time(NULL));
	PNode Root = NULL;
	CreatTree(Root);
	cout << "Tree:  " << endl;
	PrintTree(Root, 0);
	DeleteElem(Root, 5);
	cout << "New Tree:  " << endl;
	PrintTree(Root, 0);
	DeleteTree(Root);
	_CrtDumpMemoryLeaks();
	system("Pause");
  return 0;
}
