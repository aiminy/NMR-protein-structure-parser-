//#include "protein.h"
#include "linked_list.h"

//this function insert a atom node, if there is no atom node create anew one
//if there already is a atom node, then insert a new atom node after this atom node
template <class Item>  
void insert_an_node(node<Item>*& head_ptr,node<Item>*& previous_ptr,Item*& entry)
{
   node<Item> *insert_ptr;

   if(previous_ptr==NULL)
  {
    previous_ptr=new node<Item>;
    previous_ptr->data=entry;
    head_ptr=previous_ptr;
  }
 else
   {
     insert_ptr=new node<Item>;
     insert_ptr->data=entry;
     previous_ptr->next=insert_ptr;
     previous_ptr=insert_ptr;
   }

     previous_ptr->next=NULL;
   
}

//this function print atom node linked list
template <class Item>
void display_node_list(node<Item> *start_ptr)
{

  ofstream new_pdb("output_list.txt",ios::out);

  if(start_ptr == NULL)
   {cout<<"the list is empty ! "<<endl;}
  else
   {
    while(start_ptr!=NULL)
   {
    new_pdb<<start_ptr->data;
    start_ptr=start_ptr->next;
   }
    }

   new_pdb.close();
}

template<class Item>
int linked_list_length(node<Item>*& head_ptr)
{
  node<Item>* cursor_ptr;
  int length;

  length=0;
  for(cursor_ptr=head_ptr;cursor_ptr!=NULL;cursor_ptr=cursor_ptr->next)
  {
   length=length+1;
  }

  return length;


}

