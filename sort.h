#pragma once

#define LOADSIZE 4

class pSort {

public:
   typedef enum {
   	BEST,
   	QUICK,
   	MERGE,
   	RADIX
   }  SortType;

   typedef struct {
      int key;
      char payload[LOADSIZE];
   } dataType;

   void init();
   void print(dataType *data, int low, int high);
   void close();
   void sort(dataType *data, int ndata, SortType sorter=BEST); 
   void swap(dataType *data1, dataType *data2);
   int partition (dataType *data, int low, int high); 
   void quick_sort(dataType *data, int low, int high);
   int getMaximum(dataType *data, int low, int high); 
   void unit_sort(dataType *data, int low, int high, int exp); 
   void radix_sort(dataType *data, int low, int high);
   void merge(dataType *data, int low, int center, int high);
   void merge_sort(dataType *data, int low, int high);
   void heapify(dataType *data, int n, int i); 
   void heap_sort(dataType *data, int n); 
   void Rearrange(dataType *data, int ndata, int size, int index, int offset); 


};
