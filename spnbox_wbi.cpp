#include <iostream>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <x86intrin.h>

#define readSize  0x400000000
#define blocknum 128
#define l 4
#define ln 4
#define R 16

#define ltime(a) (((a) << (1)) ^ ((((a) >> (31)) & 1) * (0x8d)))

using namespace std;

void nonlinearLayer(uint32_t *state, char* memblock){
    uint64_t tmp;
    for(int j=0;j<ln;j++){
        tmp=(uint64_t)4*state[j]; 
        memcpy(&state[j],&memblock[tmp],sizeof(int));
    }
} 

void linearLayer(uint32_t *state)  //theta
{
    uint32_t temp[4]={};

    temp[0] = state[0]; temp[1] = state[1]; temp[2] = state[2]; temp[3] = state[3];

    state[0] = temp[0] ^ ltime(temp[1]^temp[3]) ^ ltime(ltime(temp[2]^temp[3]));
    state[1] = temp[1] ^ ltime(temp[0]^temp[2]) ^ ltime(ltime(temp[0]^temp[3]));
    state[2] = temp[2] ^ ltime(temp[1]^temp[3]) ^ ltime(ltime(temp[0]^temp[1]));
    state[3] = temp[3] ^ ltime(temp[0]^temp[2]) ^ ltime(ltime(temp[1]^temp[2]));
} 

void affineLayer(uint32_t *state, int round) //sigma_r
{
    for(int j=0;j<ln;j++){state[j]^=((round-1)*ln+j+1);}
}

void enc(uint32_t *state, char* memblock){
    for(int i= 1;i<=R;i++){
        nonlinearLayer(state, memblock);
        linearLayer(state);
        affineLayer(state, i);
    }
}

void generatemessage(uint32_t input[blocknum][ln]){
	for(int i=0; i<blocknum; i++){
        for(int j=0;j<ln;j++){
            input[i][j]=(rand()&0xFFFFFFFF);
        }
	}
}

int main(){
    uint32_t input[blocknum][ln]={};
    uint32_t state[blocknum][ln]={};

    int iter=100000;
    
    uint64_t start_time,end_time,cyc=0;

    ifstream table("~/table/sbox32.bin", ios::in | ios::binary );

    if(table.is_open()){
        char* memblock=new char[readSize];
        table.read(memblock,readSize);
        table.close();

        generatemessage(input);

        for(int count=0;count<iter;count ++){
            for(int j=0;j<blocknum;j++){
                for(int i=0;i<ln;i++){state[j][i]=input[j][i];} 
            }
            start_time=_rdtsc();
            for(int j=0;j<blocknum;j++){
                enc(state[j], memblock);
            }
            end_time=_rdtsc();
            cyc+=end_time-start_time;
        }

        delete[] memblock;
    }

    uint64_t aver=cyc/iter;

    cout << "average cost " << aver << "CPU cycles for " << iter << " tests \n";
	uint64_t cpb = aver/(4*ln*blocknum);
	cout << "average CPB for decryption is " << cpb << "\n";

    return 0;
}