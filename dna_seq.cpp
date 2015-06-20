#include<iostream>
#include<cstdlib>
#include<string>
#include<ctime>
using namespace std;

// to create random DNA sequences..

int main()
{
string out[21];
string array[4] = "ATCG";
int num, i;
int store[21]; 
for (i=0; i<21; i++)
 {
   srand( time(NULL) );
   bool check;
   do{
   num = rand() % (4 + 1);  
   store[i] = num; 
  check = true;
  for (int j=0; j<i; j++){
  if (num == store[j])
   { check = false;
     break;
   }
 }
  }while(!check);
// store[i] = num;
 out[i] = array[num];
 cout<<out[i];
}
return 0;
}  
