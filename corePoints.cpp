#include <cstring>
#include <iostream>
//-------------------

using namespace std;

int get_cores(int* image1d, int sizeOf1d){
    int availble = sizeOf1d/8;
    int border[256];
    int epslon=3, mindensty=1, minPts=1, Totalmindensty = 2;
    int k =0;
    int histoImage[256];

    for(int i=0;i<256;i++){
        histoImage[i] = 0;
        border[i] = 0;
    }

    for(int i=0; i < sizeOf1d; i++){

        histoImage[image1d[i]]++;
    }

    for(int i=0; i<256; i++){

        int sumofdensties = 0;
        int sumOfpoints = 0;
        for(int j=-epslon ; j<=epslon; j++){

            if(i + j < 0)continue;
            if(i + j > 255)continue;

            if(histoImage[i + j] >= mindensty){
                sumOfpoints++;
            }
            sumofdensties += histoImage[i + j];
        }
        if(sumOfpoints > minPts){
            border[i] = sumofdensties;
        }
    }


    for(int i=0; i<256; i++)
        if(border[i] >= Totalmindensty){
            k++;
            cout<<border[i]<<"   " << i<<endl;
        }
    
    



    return k;
}

 int main(){

	int imagedata[16]={68, 178,  84,  68, 125, 148, 205,  58, 100,  33, 100,6,39,111,93,75};
    int s=0;
    // for(int i=0;i<16;i++){
        
    //     s+=imagedata[i];

    // }
    // cout<<s/255<<endl;
	cout<<get_cores(imagedata, 16)<<endl;





    return 0;
 }

 