#include <iostream>
#include "sort.h"
#include <mpi.h>
#include <unistd.h>

#define LOADSIZE	4
using namespace std;
//********************QUICK SORT*******************//
void pSort::swap(dataType *data1, dataType *data2)
{ 
     pSort::dataType *t = new pSort::dataType[1];
     *t = *data1; 
     *data1 = *data2; 
     *data2 = *t; 
} 

int pSort::partition (dataType *data, int low, int high) 
{ 
    int pivot,i;
    i = (low - 1); 
    pivot = (data+high)->key; 
 
    for (int j = low; j <= high- 1; j++) 
    { 
        
        if ((data+j)->key <= pivot) 
        { 
            i++;    
            swap((data+i), (data+j)); 
        } 
    } 
    swap((data+i+1), (data+high)); 
    return (i + 1); 
} 

void pSort::quick_sort(dataType *data, int low, int high) 
{ 
// arr= 1 2 5 3 7 10 2 1 


    if (low < high) 
    {        
        int part = partition(data, low, high); 
        quick_sort(data, low, part - 1); 
        quick_sort(data, part + 1, high); 
    } 
} 
/*********************RADIX SORT***********************/
int pSort::getMaximum(dataType *data, int low, int high) 
{ 
    int mx = (data+low)->key; 
    for (int i = low; i < high; i++) 
        if ((data+i)->key > mx) 
            mx = (data+i)->key; 
    return mx; 
} 

void pSort::unit_sort(dataType *data, int low, int high, int exp) 
{ 
   
    pSort::dataType *output = new pSort::dataType[high-low];
    int i, count[10] = { 0 }; 
  
    
    for (i = low; i < high; i++) 
        count[((data+i)->key  / exp) % 10]++; 
    for (i = 1; i < 10; i++) 
        count[i] += count[i - 1]; 
  
    for (i = high-1; i >= low; i--) { 
        (output+(count[((data+i)->key / exp) % 10] - 1))->key = (data+i)->key; 

        *(output+(count[((data+i)->key / exp) % 10] - 1))->payload = *(data+i)->payload; 
        count[((data+i)->key / exp) % 10]--; 
    } 
  
    count[0]=0;
    for (i = low; i < high; i++) {
        (data+i)->key = (output+count[0])->key;
        *(data+i)->payload = *(output+count[0])->payload; 
        count[0]++;
    }
}

void pSort::radix_sort(dataType *data, int low, int high)
{
   int m = getMaximum(data, low,high); 
   int i;   
    for ( i = 1; m / i > 0; i *= 10) 
        unit_sort(data, low,high, i); 
    pSort::dataType *TEMP = new pSort::dataType[high]; 
    m=0;
    for(i = high-1; i >=0; i--){
	if((data+i)->key<0){
        *(TEMP+m) = *(data+i);
	m++;
		}
	}

    for(i = 0; i < high; i++){
	if((data+i)->key>=0){
        *(TEMP+m) = *(data+i);
	m++;
		}
	}
    for(i = 0; i < high; i++){
        *(data+i) = *(TEMP+i);
	}


} 

/**************************MERGE SORT*****************/
void pSort::merge(dataType *data, int low, int center, int high)
{
//arr = [1,2,3,4,5,6,7,6,5,4,3,2,1]
    int l1,l2,i,j,k;
    l1 = center - low+1;
    l2 = high  - center;
   
    pSort::dataType *LEFT = new pSort::dataType[l1]; 
    pSort::dataType *RIGHT = new pSort::dataType[l2]; 
    
    for(int i = 0; i < l1; i++)
        *(LEFT+i) = *(data+low + i);
   
        
    for(int j = 0; j < l2; j++)
        *(RIGHT+j) = *(data+center + j+1);

    i = 0;  
    j = 0;     
    k = low;
     

    while (i < l1 && j < l2)
    {
        if ((LEFT+i)->key <= (RIGHT+j)->key) 
        {
            *(data+k)=*(LEFT+i);
            i++;
        }
        else
        {
            *(data+k)=*(RIGHT+j);
            j++;
        }
        k++;
    }

    while (i < l1) 
    {
        *(data+k)=*(LEFT+i);
        i++;
        k++;
    }

    while (j < l2)
    {
        *(data+k)=*(RIGHT+j);
        j++;
        k++;
    }

}

void  pSort::merge_sort(dataType *data, int low, int high)
{
//arr = [2];
// arr1 = [1,3,5,5] arr2= [2,4,6,8]
//arr= [1,2,3,4,5,5,6,8]
    if (low < high)
    {        
       
        int m = low + (high - low) / 2;
 
       
        merge_sort(data, low, m);
        merge_sort(data, m + 1, high);
 
        merge(data, low, m, high);
    }
}

/*************************BEST SORT*******************/
void pSort::heapify(dataType *data, int n, int i) 
{ 
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;
  

    if (l < n && (data+l)->key > (data+largest)->key) 
        largest = l; 
  

    if (r < n && (data+r)->key > (data+largest)->key) 
        largest = r; 
  

    if (largest != i) { 
        swap(data+i, (data+largest)); 

        heapify(data, n, largest); 
    } 
} 

void pSort::heap_sort(dataType *data, int n) 
{ 
    int i;
    for ( i = n / 2 - 1; i>=0; i--){
        if(i < 0)
        break;
        heapify(data, n, i);
    } 
         
  
    
    for (int i = n - 1; i >= 0; i--) {        
        swap(data, data+i);        
        heapify(data, i, 0); 
    } 
}

void pSort::Rearrange(dataType *data, int ndata, int size, int index, int offset){
int length= ndata/size+2;
pSort::dataType *Temp = new pSort::dataType[length];
int k,j,i=0;
for (i=0;i<length;i++){
     *(Temp+i) = *(data + index+ i);  
 }

 i=0;
 while(offset){
     *(data+i)= *(Temp);
     i=i+1;
     offset=offset-1;
}
if(i!=0){
    k=1;
}else{
    k=0;
}

 while(i<ndata){
     j=0;
     for(j=0;j<size && i<ndata; j++){
         *(data+i)= *(Temp+k);
         i++;
     }
     k=k+1;   
 }

}

void pSort::init()
   {   
   MPI_Init(NULL, NULL);
   }

void pSort::print(dataType *data, int low, int high)
{
      for (int i=low;i<high;i++)
      {
        cout<<(data+i)->key<<" payload is:   ";
        for(int j=0;j<LOADSIZE;j++)
        cout<<(data+i)->payload[j];
        cout<<endl;
      }

}
void pSort::gatherSort(dataType *data,int ndata, MPI_Datatype Type, int rank, int size){
pSort::dataType *recv_data = new pSort::dataType[ndata *size];

//if (rank==0)
//print(data, 0, ndata);

//cout<<"data before broadcast"<<endl;
MPI_Allgather(
    data,
    ndata,
    Type,
    recv_data,
    ndata,
    Type,
    MPI_COMM_WORLD);

usleep(100);
//if (rank==0){
int lower=0;
int center=-1;
int high=ndata-1;
while(1){
    center = center+ndata;
    high=high+ndata;
    //cout<<"center n n data is "<<center<<"  "<<high<<endl;
    merge(recv_data, lower, center, high);
    if (high==(size)*ndata-1){
        //cout<<"break should be at cen n high"<<center<<"  "<<high<<endl;
        break;

    }
        
}
//}
usleep(100);
int index=0;
for (int i=ndata*rank;i<(ndata*rank)+ndata;i++)
    {
    *(data+index)=*(recv_data+i);
    index++;       
    }

 delete recv_data;




}
void pSort::merge2(dataType *data, dataType *recv_data, int ndata, int rndata, bool low){

       //int rank, size;
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   //MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
   //cout<<"merge2 is called by "<<rank<<" low is "<<low<<endl;
   int rank, size, lower;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
int i,j,k;
 pSort::dataType *TEMP = new pSort::dataType[ndata]; 
  //copy of data
    for( i = 0; i < ndata; i++)
        *(TEMP+i) = *(data + i);


    /*    if (rank == 2){
            cout<<"data n rcv data"<<endl;
         usleep(1000);
     //cout <<interations<< "  in loop sorting from " << rank << endl;
     for (int i = 0; i < ndata; i++) {
        cout << data[i].key << "(" << data[i].payload[0] << "," <<
                                            data[i].payload[1] << "," <<
                                            data[i].payload[2] << "," <<
                                            data[i].payload[3] << ")" << endl;
    }cout<< endl<< endl;
 }
         if (rank == 2){
         usleep(1000);
     //cout <<interations<< "  in loop sorting from " << rank << endl;
     for (int i = 0; i < rndata; i++) {
        cout << recv_data[i].key << "(" << recv_data[i].payload[0] << "," <<
                                            recv_data[i].payload[1] << "," <<
                                            recv_data[i].payload[2] << "," <<
                                            recv_data[i].payload[3] << ")" <<  endl;
    }cout<< endl<< endl;
 }*/
    if (low==true){
    i = 0;  
    j = 0;     
    k = 0;
    while (k < ndata) //need to modify in case of uneven length
    { if( j < rndata){
        if ((TEMP+i)->key <= (recv_data+j)->key) 
        {
            *(data+k)=*(TEMP+i);
            i++;
        }
        else
        {
            *(data+k)=*(recv_data+j);
            j++;
        }
    }else{
            //cout<<"low merging "<<endl;
            *(data+k)=*(TEMP+i);
            i++;
    }
        k++;
    }
    }else
    {
        //cout<<"we are here now "<<rank<<endl;
    i = ndata-1;  
    j = rndata-1;     
    k = ndata-1;
    while (k >= 0) //need to modify in case of uneven length
    { 
        if(!(j<0)){
        //print(TEMP, 0, ndata);
        if ((TEMP+i)->key >= (recv_data+j)->key) 
        {
            *(data+k)=*(TEMP+i);
            i--;
        }
        else
        {
            *(data+k)=*(recv_data+j);
            j--;
        }
    }else{
       //cout<<"high merge "<<j<<endl;
        *(data+k)=*(TEMP+i);
        i--;
    }
        k--;
    }
    }
    
  delete TEMP; 
        /*if (rank == 2){
            cout<<"after merging"<<endl;
         usleep(1000);
     //cout <<interations<< "  in loop sorting from " << rank << endl;
     for (int i = 0; i < ndata; i++) {
        cout << data[i].key << "(" << data[i].payload[0] << "," <<
                                            data[i].payload[1] << "," <<
                                            data[i].payload[2] << "," <<
                                            data[i].payload[3] << ")" << endl;
    }cout<< "........................................"<<endl<< endl;
 }*/
}


/* odd even merge*/
void pSort::double_merge(dataType *data, MPI_Datatype Type,  int ndata,int rank,int size)
{
bool step=true;
MPI_Status status;
//calculate max diff;
int diff;
int min_size,max_size,recv_size;
if(rank==0){
    
    min_size=ndata;
    max_size=ndata;
    for (int i=1;i<size;i++){
        MPI_Recv (&recv_size, 1, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
        if(recv_size>max_size)
        max_size=recv_size;
        
        if(recv_size<min_size)
        min_size=recv_size;
        
    }
}else{
   MPI_Send(&ndata, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);    
   //MPI_Recv (&diff, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
        
}
MPI_Barrier(MPI_COMM_WORLD);
if(rank==0){

  for (int i=1;i<size;i++){
        diff=max_size-min_size+1;
        //cout<<"send "<<i<<endl;
        MPI_Send(&diff, 1, MPI_INT, i, 1, MPI_COMM_WORLD);      
    }
}else{
     MPI_Recv (&diff, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
     // cout<<"recv "<<rank<<endl;
        
}
if(rank==0)
cout<<diff<<"/////////////"<<endl;
int interations, recv_ndata;
interations=size*diff*2;
//if (size%2==0)
//interations--;
while(interations>0){
 /*       if (rank == 3){
         usleep(1000);
     cout <<interations<< "  in loop sorting from " << rank << endl;
     for (int i = 0; i < ndata; i++) {
        cout << data[i].key << "(" << data[i].payload[0] << "," <<
                                            data[i].payload[1] << "," <<
                                            data[i].payload[2] << "," <<
                                            data[i].payload[3] << ")" << endl;
    }
 }*/
if(step==true){
if(rank%2==0){
    if(rank+1<size){ 
        //cout<<"rank from if "<<rank<<endl;    
        //mpi send(from rank to rank+1)
        MPI_Send (data, ndata, Type, rank+1, 1, MPI_COMM_WORLD);
        //mpi recv(from rank+1 to rank)

        MPI_Probe(rank+1, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, Type, &recv_ndata);
        //recv_ndata=recv_ndata/2;
        //cout<<"...........if "<<recv_ndata<<endl;
        pSort::dataType *recv_data = new pSort::dataType[recv_ndata];
        MPI_Recv (recv_data, recv_ndata, Type, rank+1, 1, MPI_COMM_WORLD, &status);
        //cout<<step<<"if interation: "<<interations<<" send from "<<rank<<" to "<<rank+1<<endl;
        //if(rank==0)
        //print(recv_data, 0, recv_ndata);
        //cout<<"recv ndata is "<<recv_ndata<<endl;
        //usleep(100000000);
        //get the lower values here
        //cout<<rank<<"calling merge low in iteration"<<interations<<endl;
        merge2(data, recv_data, ndata, recv_ndata, true);
        delete recv_data;
    }
    

}else{    
    //mpi recv(from rank-1 to rank)
    //cout<<"rank from else "<<rank<<endl; 
    MPI_Probe(rank-1, 1, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, Type, &recv_ndata);
    //recv_ndata=recv_ndata/2;
    //cout<<".........else..."<<recv_ndata<<endl;
    pSort::dataType *recv_data = new pSort::dataType[recv_ndata*3];
    MPI_Recv (recv_data, recv_ndata, Type, rank-1, 1, MPI_COMM_WORLD, &status);
    //cout<<step<<"else interation: "<<interations<<" send from "<<rank<<" to "<<rank-1<<endl;
    //print(recv_data, 0, ndata);
    //mpi send(from rank to rank-1)
    MPI_Send (data, ndata, Type, rank-1, 1, MPI_COMM_WORLD);
    //get the high values here
    //cout<<rank<<"calling merge HIGH in iteration"<<interations<<endl;
    merge2(data, recv_data, ndata, recv_ndata, false);
    delete recv_data;
    }
MPI_Barrier(MPI_COMM_WORLD);
interations--;
step=false;    
}
//break;/////////////////////////////REMOVE THIS
if(interations<=0)
break;

if(step==false){
if(rank%2==0){
    if(rank-1>=0){
        //cout<<step<<"interation: "<<interations<<" send from "<<rank<<" to "<<rank-1<<endl;
        //mpi send(from rank to rank-1)
        MPI_Send (data, ndata, Type, rank-1, 1, MPI_COMM_WORLD);
        //mpi recv(from rank-1 to rank)
        MPI_Probe(rank-1, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, Type, &recv_ndata);
        //recv_ndata=recv_ndata/2;
        pSort::dataType *recv_data = new pSort::dataType[recv_ndata];
        MPI_Recv (recv_data, recv_ndata, Type, rank-1, 1, MPI_COMM_WORLD, &status);
        //get the high values here
        //cout<<rank<<"calling merge HIGH in iteration"<<interations<<endl;
        merge2(data, recv_data, ndata, recv_ndata, false);
        delete recv_data;
        //print(recv_data, 0, ndata);     
        }
    
}else{
    if(rank+1<size){    
        //cout<<step<<"interation: "<<interations<<" send from "<<rank<<" to "<<rank+1<<endl;
        // mpi recv(from rank+1 to rank)
        MPI_Probe(rank+1, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, Type, &recv_ndata);
        //recv_ndata=recv_ndata/2;
        pSort::dataType *recv_data = new pSort::dataType[recv_ndata];
        MPI_Recv (recv_data, recv_ndata, Type, rank+1, 1, MPI_COMM_WORLD, &status);
        //print(recv_data, 0, ndata);
        //mpi send(from rank to rank+1)
        MPI_Send (data, ndata, Type, rank+1, 1, MPI_COMM_WORLD);
        //get the lower values here
        //cout<<rank<<"calling merge low in iteration"<<interations<<endl;
        merge2(data, recv_data, ndata, recv_ndata, true);
        delete recv_data;
    }
}
MPI_Barrier(MPI_COMM_WORLD);
interations--;
step=true;    
}
//interations--;
}
}





void pSort::close()
   {
   MPI_Finalize();
   }
void pSort::sort(dataType *data, int ndata, SortType sorter)
{
 
    
   //cout<<ndata<<endl;
   // printf("printing data in sort");
   // printf(" ndata is %d    ", ndata);
   // cout<< (data+1)->key;
   int rank, size, lower, upper,n_record,i,index, offset,count;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank); 


 // n_record=ndata/size;
 
  // if (rank==0){
    //cout<<"size "<<size<<" ndata "<<ndata<<" n_record "<<n_record<<endl;   
    //}

   // if(rank==1){
   //     delete data+(ndata-2);
  //      ndata--;
   // }
    

//    lower=n_record*rank;
//    upper=lower+n_record;  
//    if(rank==size-1)
//    {     
//       upper=ndata;
//    }
   
   lower=0;
   upper=ndata;

   index=-1;
   offset=0;
   count=1;
   //ndata=9;
   while(count<=ndata*rank){
       //cout<<"count in while is "<<count<<endl;
       count=count+size;
       index=index+1;
   }
   if ((ndata*rank)%count==0 && rank!=0){
       //cout<<"count is equal to size*rank"<<count<<endl;
       index=index-1;
   }
   //cout<<"count is "<<count<<endl;
   while(count!=ndata*rank){
      // cout<<"offset is getting changed "<<count<<" "<<ndata*rank <<endl;
       count=count-1;
       offset=offset+1;
   }
   if (offset==size){
       offset=0;
   }

//cout<<"rank is "<<rank<<" index is: "<<index<<" offset is: "<<offset<<" count is: "<<count<<endl;
// call a function according to offset and index




 //  if (rank==0)
 {
      //cout<<"rank "<<rank <<" lower "<<lower<<" upper "<<upper<<endl;

     
        //printing
        //print(data,lower,upper);
        switch(sorter)
        {
        case BEST:
        //cout<<"BEST is selected"<<endl;
        heap_sort(data+lower, upper-lower); 
        break;
        case QUICK:
        //cout<<"QUIK is selected"<<endl;
        quick_sort(data,lower,upper-1);
        break;
        case MERGE:
        //cout<<"MERGE is selected"<<endl;
        merge_sort(data,lower,upper-1);
        break;
        case RADIX:
        cout<<"RADIX is selected"<<endl;
        radix_sort(data,lower,upper);
        break;
        default:
        cout<<"no valid choise about type of sort function"<<endl;
        }
 
      
        //cout<<"print after sorting"<<endl;
        //Rearrange(data,  ndata,  size,  index,  offset);

        //print(data,lower,upper);
      
   }


//try to get all elments on all process and try to merge subparts:
//define the datatype
MPI_Datatype Type;
pSort::dataType *f_buf = new pSort::dataType[ndata *size];
f_buf[0].key =0;

int blocklengs[2]; 
MPI_Aint indices[2]; 
MPI_Datatype old_types[2];

blocklengs[0]=1;
blocklengs[1]=4;

old_types[0]=MPI_INT;
old_types[1]=MPI_CHAR;

MPI_Get_address(&((*f_buf).key), &indices[0]);
MPI_Get_address(&((*f_buf).payload), &indices[1]);

indices[1]=indices[1]-indices[0];
indices[0]=0;

MPI_Type_create_struct(2,blocklengs,indices,old_types,&Type);
MPI_Type_commit(&Type);

//if (rank==0)
//print(data, 0, ndata);
/*
if (rank==0){
print(data, 0, ndata);
cout<<"......................................"<<endl;}
usleep(1000000);

if (rank==1){
print(data, 0, ndata);
cout<<"***************************************"<<endl;}
usleep(1000000);

if (rank==2){
print(data, 0, ndata);
cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"<<endl;}
usleep(1000000);
*/
//gatherSort(data, ndata,  Type,  rank,  size);
double_merge(data , Type, ndata, rank, size);

/*
cout<<"datatype created successfully"<<endl;

if (rank==0){
print(data, 0, ndata);
cout<<"......................................"<<endl;}
usleep(1000000);

if (rank==1){
print(data, 0, ndata);
cout<<"***************************************"<<endl;}
usleep(1000000);

if (rank==2){
print(data, 0, ndata);
cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"<<endl;}
usleep(1000000);
*/
if (rank == 0){
            cout<<"just before return"<<endl;
         usleep(1000);
     //cout <<interations<< "  in loop sorting from " << rank << endl;
     for (int i = 0; i < ndata; i++) {
        cout << data[i].key << "(" << data[i].payload[0] << "," <<
                                            data[i].payload[1] << "," <<
                                            data[i].payload[2] << "," <<
                                            data[i].payload[3] << ")" << endl;
    }cout<< "........................................"<<endl<< endl;
 }

}










/*
645046143 payload is:   GZTR
639961816 payload is:   Ov4C
118849868 payload is:   T2PL
746963374 payload is:   vbVz
106884453 payload is:   hrIj
676511442 payload is:   wWqO
1470674618 payload is:   uevN
2032362113 payload is:   ilyk
889696710 payload is:   KDwi
362344524 payload is:   Ib5l
......................................
645046143 payload is:   GZTR
639961816 payload is:   Ov4C
118849868 payload is:   T2PL
746963374 payload is:   vbVz
106884453 payload is:   hrIj
676511442 payload is:   wWqO
1470674618 payload is:   uevN
2032362113 payload is:   ilyk
889696710 payload is:   KDwi
362344524 payload is:   Ib5l
***************************************
645046143 payload is:   GZTR
639961816 payload is:   Ov4C
118849868 payload is:   T2PL
746963374 payload is:   vbVz
106884453 payload is:   hrIj
676511442 payload is:   wWqO
1470674618 payload is:   uevN
2032362113 payload is:   ilyk
889696710 payload is:   KDwi
362344524 payload is:   Ib5l
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
print after sorting
print after sorting
print after sorting
106884453 payload is:   hrIj
106884453 payload is:   hrIj
106884453 payload is:   hrIj
106884453 payload is:   hrIj
118849868 payload is:   T2PL
118849868 payload is:   T2PL
118849868 payload is:   T2PL
118849868 payload is:   T2PL
362344524 payload is:   Ib5l
362344524 payload is:   Ib5l
......................................
106884453 payload is:   hrIj
106884453 payload is:   hrIj
118849868 payload is:   T2PL
118849868 payload is:   T2PL
362344524 payload is:   Ib5l
362344524 payload is:   Ib5l
639961816 payload is:   Ov4C
639961816 payload is:   Ov4C
645046143 payload is:   GZTR
645046143 payload is:   GZTR
***************************************
106884453 payload is:   hrIj
118849868 payload is:   T2PL
362344524 payload is:   Ib5l
639961816 payload is:   Ov4C
645046143 payload is:   GZTR
676511442 payload is:   wWqO
746963374 payload is:   vbVz
889696710 payload is:   KDwi
1470674618 payload is:   uevN
2032362113 payload is:   ilyk
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
