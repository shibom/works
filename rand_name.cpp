#include<iostream>
#include<cstdlib>
#include<string>
#include<fstream>
#include<ctime>
using namespace std;

int main()
{
  ifstream infile;
  ofstream outfile;
 infile.open("dark1us_index.txt");
 outfile.open("dark_8K.txt");
 const int MAX = 70000;
 const int LIMIT = 8000;
 string innames[MAX];
 string outnames[LIMIT];
 int a = 0;
 int rand_num;
 int store[LIMIT];
 while(infile.good())
 {
   getline(infile, innames[a]);
   a++;
 }
 
 for(int i = 0; i<LIMIT; i++)
 {
   srand( time(NULL));
   bool check;
   do{
   rand_num = rand() % (MAX + 1);
   check = true;
   for(int j = 0; j<i; j++){
    if (rand_num == store[j])
    {
      check=false;
      break; 
    }
   }
   }while(!check);
   store[i] = rand_num;
   outnames[i] = innames[rand_num];
 }
 
 for(int i = 0; i<LIMIT; i++){
 outfile<< outnames[i] << endl;
 }
 outfile.close();
 infile.close();
 return 0;
 }
