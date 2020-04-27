// A population evolving as per the "Wright-Fisher process", with two types/alleles (A and a) present at t = 0

#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<fstream>
#include<iomanip>
#include<math.h>

using namespace std;

int main()
{
ofstream out;
out.open("fixation_prob_pop size 10, s = 0.txt");
//out.open("/home/jyoti/Desktop/C_Codes/text files/WF_pop size: 10, s = 0.3.txt");

out<<setw(20)<<"N"<<setw(20)<<"s"<<setw(20)<<"2*N*s"<<setw(30)<<"intial #A/N "<<setw(30)<<" fixation prob. (simulations) "<<"\n";

int i, j, Pop_size, Init_A, t_max, t, n, n_exp, fix_A, df;
float ran, ran_A, ran_a, l_bound, u_bound, sum, denom, sum_fit, mu_A, mu_a, s; 

char type_inter[1100], types[10000];
int tot_A[11000];
float w[10000], f[10000], w_inter[10000];

srand (time(NULL));

//***************************************************************************

//** Can also be taken as an input from users

t_max = 1000;                    // time (or number of generations) given for population to evolve 
Pop_size = 10;                // population size
Init_A = 8;                   // number of allele/type A one wants to start with, i.e., present at t = 0

n_exp = 10000;                  // number of times the same population was allowed to evolve, keeping initial conditions same

s = 0.0;                     // selection coefficient

mu_A = 0.0000;                  // mutation rate (A -> a) 
mu_a = 0.0000;                  // mutation rate (a -> A) 

//***************************************************************************

fix_A = 0;
f[0]=0;
tot_A[0] = Init_A; 

//out<<setw(15)<<"#populations "<<setw(15)<<" t "<<setw(15)<<" freq. of type A(t) "<<endl;

//***************************************************************************

for(n=1; n<(n_exp + 1); n++)
{

//out<<setw(15)<<n<<setw(15)<<0<<setw(15)<<(0.5)<<endl;

for(i=1; i<=(Pop_size); i++)
{

//********** 'A' allele present in 'Init_A' copies in a population of size 'Pop_size' at t = 0 ********

if(i<=Init_A)
{
types[i] = 'A';    
w[i] = 1 + s;                           // w[i] are the absolute fitnesses of the individuals
}
else
{
types[i] = 'a'; 
w[i] = 1;
}
}

//***************************************************************************

for(t=1; t<(t_max+1); t++)
tot_A[t]=0;                             // tot_A[t]: total number of allele/type 'A' present after 't' number of generations

for(t=1; t<(t_max+1); t++)
{
denom = (Pop_size + s*tot_A[t-1]);
sum_fit=0;

for(i=1; i<=(Pop_size); i++)
{
sum_fit+=w[i];                          // ADDING up the w[i] 
f[i]=sum_fit;
}

for(i=1; i<=(Pop_size); i++)            // composition of the population as it evolves, shaped by different evolutionary forces!!
{
ran = (rand()/(1.0*RAND_MAX));
for(j=1; j<=(Pop_size); j++)
{
l_bound = f[j-1]/(denom); 
u_bound = f[j]/(denom); 

if(((l_bound < ran)&&(ran < u_bound))||(ran==l_bound))
{
w_inter[i] = w[j];
type_inter[i] =  types[j];
if(type_inter[i]=='A')
tot_A[t]+=1;
break;
}
}  
} 

//**************     Adding mutations         ***************************

if((mu_A!=0)||(mu_a!=0))           // 'mu_A' and 'mu_a' are the mutation rates of alleles/types 'A' and 'a' respectively 
{
for(i=1; i<=(Pop_size); i++)
{
ran_A = (rand()/(1.0*RAND_MAX));
ran_a = (rand()/(1.0*RAND_MAX));

if(type_inter[i]=='a')  
{
if((ran_a < mu_a)&&(ran_A > mu_A))  // a -> A
{     
w_inter[i] = 1 + s;
type_inter[i] = 'A';
tot_A[t]+=1;             
}

}

if(type_inter[i]=='A')  
{

if((ran_A < mu_A)&&(ran_a > mu_a))  // A -> a
{     
w_inter[i] = 1;
type_inter[i] = 'a';
tot_A[t]-=1;             
}

}

}
}

//**************     Adding mutations         ***************************

for(i=1; i<=(Pop_size); i++)
{
w[i] = w_inter[i];
types[i] = type_inter[i];
}

//out<<setw(15)<<n<<setw(15)<<t<<setw(15)<<(1.0*tot_A[t]/(Pop_size))<<endl;

if((tot_A[t]==0)||(tot_A[t]==(Pop_size)))                    // one of the two alleles has reached fixation
{
df+=1;
break;
}
} 

if((tot_A[t]==(Pop_size)))                                   // allele 'A' has reached fixation
fix_A+=1;
} 

out<<setw(20)<<Pop_size<<setw(20)<<s<<setw(20)<<(Pop_size)*(2*s)<<setw(30)<<((1.0*Init_A)/(Pop_size))<<setw(30)<<(1.0*fix_A)/(1.0*n_exp)<<"\n";
out.close();

cout<<df;
//out<<" # pop size = "<<Pop_size<<" selection coefficient = "<<s<<endl;
//cout<<" Fixation fixability "<<(1.0*fix_A)/(1.0*n_exp)<<"    "; 
} // main function

