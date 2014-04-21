/*
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <math.h>

using namespace std;
*/
#include "atom.h"

#ifndef LINKED_LIST_H
#define LINKED_LIST_H

//#include "protein.h"

//#include "aa.h"

template <class Item>
struct node
{
 Item data;
 node *next;
};

template<class Item>
void insert_an_node(node<Item>*& head_ptr,node<Item>*& previous_ptr,Item& entry);

template<class Item>
void display_node_list(node<Item> *start_ptr);

template<class Item>
int linked_list_length(node<Item> *head_ptr);

template<class Item>
double linked_list_mean(node<Item> *head_ptr);

template<class Item>
double linked_list_sd(node<Item> *head_ptr);

#include "linked_list.template"
#endif
