#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <mpi.h>
#include <set>
#include <ctime>// include this header 
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>


#define eps 1e-6
using namespace std;
using namespace msclr::interop;
int rnk;
int sz;
int K = 6;
int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}
void initialize_Centroids(double* &Centroid,int * image,int W , int H)
{
	
	set<int>st;
	for (int i = 0; i < K; i++)
	{
		int tmp = image[(rand() % (W * H))];
		if (st.find(tmp) == st.end())
		{
			Centroid[i] = tmp;
			st.insert(tmp);
		}
		else i--;
	}
}

void pre_scatterv(int* &send_counts_scatterv,int* &displs_scatterv , int total)
{
	
	
	int basic = total / sz;
	int remain = total % sz;
	int sum = 0;
	int tmp;
	for (int i=0;i<sz;i++)
	{
		displs_scatterv[i] = sum;
		tmp = min(1, remain);
		remain -= tmp;
		send_counts_scatterv[i] = basic + tmp;
		sum += basic + tmp;
	}
}



int main()
{

	/*
	
	
	*/
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
	//All global variables 
	int ImageWidth = 4, ImageHeight = 4;

	int start_s, stop_s, TotalTime = 0;
	int* imageData=NULL;
	//initial values of centroids
	double* Centroid = new double[K];;
	double* new_Centroid = new double[K];;

	int* send_counts_scatterv = new int[sz];;
	int* displs_scatterv = new int[sz];


	//All Local Variables 


	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";

	imagePath = marshal_as<System::String^>(img);
	imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
	if (rnk == 0)
	{
		start_s = clock();
		//first step intialize centroids with random numbers 
		initialize_Centroids(Centroid, imageData, ImageWidth, ImageHeight);
		//filling scatterv attributes 
		pre_scatterv(send_counts_scatterv, displs_scatterv, ImageHeight * ImageWidth);
	}

	//giving meta scatterv data and image pixels to all proccessors 
	MPI_Bcast(send_counts_scatterv, sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs_scatterv, sz, MPI_INT, 0, MPI_COMM_WORLD);
	int g = 0;
	int pixel_cnt = send_counts_scatterv[rnk];
	int* pixels = new int[pixel_cnt];
	MPI_Scatterv(imageData, send_counts_scatterv, displs_scatterv, MPI_INT, pixels, pixel_cnt, MPI_INT, 0, MPI_COMM_WORLD);

		//some important data for every iteration 
		int* sum_of_distances_per_centroid = new int[K];
		int* count_pixel_per_centroid = new int[K];
		int* centroid_id_per_pixel = new int[pixel_cnt];
		int* cluster_size = new int[K];
		int* cluster_inner_sum = new int[K];
		
	while (true) 
	{
		MPI_Bcast(Centroid, K, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rnk == 0)
		{
			for (int i = 0; i < K; i++)
				cout << Centroid[i] << " ";
			cout << endl;

		}

		memset(sum_of_distances_per_centroid, 0, sizeof(int) * K);
		memset(count_pixel_per_centroid, 0, sizeof(int) * K);
		for (int i=0;i<pixel_cnt;i++)
		{
			int mn = 1000000;
			int id = 0;
			for (int j=0;j<K;j++)
			{
				if (mn>abs(Centroid[j]-pixels[i]))
				{
					mn = Centroid[j] - pixels[i];
					id = j;
				}
			}
			sum_of_distances_per_centroid[id] += pixels[i];
			count_pixel_per_centroid[id]++;
			centroid_id_per_pixel[i] = id;
			//if (centroid_id_per_pixel[i] > 2 || centroid_id_per_pixel[i] < 0)cout << "7a7a" << endl;
		}
		
		MPI_Reduce(sum_of_distances_per_centroid, cluster_inner_sum, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(count_pixel_per_centroid, cluster_size, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rnk == 0)
		{
			int same = 1;
			for (int i=0;i<K;i++)
			{
				new_Centroid[i] = cluster_inner_sum[i] / (double)(cluster_size[i]);
				new_Centroid[i] = min(255, max(new_Centroid[i], 0));
				if (abs(new_Centroid[i] - Centroid[i]) > eps)same = 0;
				Centroid[i] = new_Centroid[i];
			}
			
			for (int i = 1; i < sz; i++)
				MPI_Send(&same, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

			if (same)
			{
				break;
			}
				
				
		}
		else
		{
			int same;
			MPI_Status st;
			MPI_Recv(&same, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&st);
			if (same)break;
		}

		
	}
	int* cluster = new int[ImageHeight * ImageWidth];

	MPI_Gatherv(centroid_id_per_pixel, pixel_cnt, MPI_INT, cluster, send_counts_scatterv, displs_scatterv, MPI_INT, 0, MPI_COMM_WORLD);
	if (rnk == 0)
	{
		

		for (int i = 0; i < ImageHeight * ImageWidth; i++)
		{
			imageData[i] = Centroid[cluster[i]];
		}

		
		
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		createImage(imageData, ImageWidth, ImageHeight, 1);
		cout << "time: " << TotalTime << endl;
	}
	free(imageData);
	MPI_Finalize();
	return 0;

}



