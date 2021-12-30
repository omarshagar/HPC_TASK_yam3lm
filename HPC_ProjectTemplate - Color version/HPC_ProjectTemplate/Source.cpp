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
int K =4;
void inputImage(int* w, int* h, System::String^ imagePath ,int*&Red ,int*&Green ,int*&Blue) //put the size of image in w & h
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
	 Red = new int[BM.Height * BM.Width];
	 Green = new int[BM.Height * BM.Width];
	 Blue = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);
			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

		}

	}
}


void createImage(int *R,int *G ,int *B, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (R[i * width + j] < 0)
			{
				R[i * width + j] = 0;
			}
			if (R[i * width + j] > 255)
			{
				R[i * width + j] = 255;
			}
			if (G[i * width + j] < 0)
			{
				G[i * width + j] = 0;
			}
			if (G[i * width + j] > 255)
			{
				G[i * width + j] = 255;
			}
			if (B[i * width + j] < 0)
			{
				B[i * width + j] = 0;
			}
			if (B[i * width + j] > 255)
			{
				B[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(R[i * MyNewImage.Width + j], G[i * MyNewImage.Width + j],B[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}
int convert(int rank, int col)
{
	return rank * 3 + col;
}
void initialize_Centroids(double*& Centroid, int* R , int *G ,int *B, int W, int H)
{
	set<pair<double , pair<double ,double >>>st;
	for (int i = 0; i < K; i++)
	{
		int h = (rand() % (W * H));
		double  tmpR = R[h];
		double  tmpG = G[h];
		double  tmpB = B[h];
		if (st.find({ tmpR,{tmpG,tmpB} }) == st.end())
		{
			Centroid[convert(i,0)] = tmpR;
			Centroid[convert(i, 1)] = tmpG;
			Centroid[convert(i, 2)] = tmpB;
			st.insert({ tmpR,{tmpG,tmpB} });
		}
		else i--;
	}
}



void pre_scatterv(int*& send_counts_scatterv, int*& displs_scatterv, int total)
{


	int basic = total / sz;
	int remain = total % sz;
	int sum = 0;
	int tmp;
	for (int i = 0; i < sz; i++)
	{
		displs_scatterv[i] = sum;
		tmp = min(1, remain);
		remain -= tmp;
		send_counts_scatterv[i] = basic + tmp;
		sum += basic + tmp;
	}
}


double Distance(double * Centroids , int * R , int *G , int *B , int rnk , int idx  )
{
	double RR = Centroids[convert(rnk, 0)] - R[idx];
	double BB = Centroids[convert(rnk, 1)] - B[idx];
	double GG = Centroids[convert(rnk, 2)] - G[idx];
	RR *= RR;
	BB *= BB;
	GG *= GG;
	return sqrt(RR + BB + GG);
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
	int* imageData = NULL;
	int* R= NULL;
	int* G= NULL;
	int* B= NULL;
	//initial values of centroids
	double* Centroid = new double[K*3];;
	double* new_Centroid = new double[K*3];;

	int* send_counts_scatterv = new int[sz];;
	int* displs_scatterv = new int[sz];


	//All Local Variables 


	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test2.jpg";

	imagePath = marshal_as<System::String^>(img);
	inputImage(&ImageWidth, &ImageHeight, imagePath,R,G,B);
	if (rnk == 0)
	{
		start_s = clock();
		//first step intialize centroids with random numbers 
		initialize_Centroids(Centroid, R,G,B, ImageWidth, ImageHeight);
		//filling scatterv attributes 
		pre_scatterv(send_counts_scatterv, displs_scatterv, ImageHeight * ImageWidth);
	}

	//giving meta scatterv data and image pixels to all proccessors 
	MPI_Bcast(send_counts_scatterv, sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs_scatterv, sz, MPI_INT, 0, MPI_COMM_WORLD);
	int g = 0;
	int pixel_cnt = send_counts_scatterv[rnk];
	int* Rpixels = new int[pixel_cnt];
	int* Gpixels = new int[pixel_cnt];
	int* Bpixels = new int[pixel_cnt];
	MPI_Scatterv(R, send_counts_scatterv, displs_scatterv, MPI_INT, Rpixels, pixel_cnt, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(G, send_counts_scatterv, displs_scatterv, MPI_INT, Gpixels, pixel_cnt, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(B, send_counts_scatterv, displs_scatterv, MPI_INT, Bpixels, pixel_cnt, MPI_INT, 0, MPI_COMM_WORLD);

	//some important data for every iteration 
	int* sum_of_distances_per_centroid = new int[3*K];
	int* count_pixel_per_centroid = new int[K];
	int* centroid_id_per_pixel = new int[pixel_cnt];
	int* cluster_size = new int[K];
	int* cluster_inner_sum = new int [3*K];
	int ep = 70;
	while (ep--)
	{
		MPI_Bcast(Centroid, 3*K, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (rnk == 0)
		{
			for (int i = 0; i <3* K; i++)
				cout << Centroid[i] << " ";
			cout << endl;

		}

		memset(sum_of_distances_per_centroid, 0, sizeof(int) * 3*K);
		memset(count_pixel_per_centroid, 0, sizeof(int) * K);
		for (int i = 0; i < pixel_cnt; i++)
		{
			double  mn = 1000000;
			int id = 0;
			for (int j = 0; j < K; j++)
			{
				double D = Distance(Centroid, Rpixels, Gpixels, Bpixels, j, i);
				if (mn >D)
				{
					mn = D;
					id = j;
				}
			}
			sum_of_distances_per_centroid[convert(id,0)] += Rpixels[i];
			sum_of_distances_per_centroid[convert(id, 1)] += Gpixels[i];
			sum_of_distances_per_centroid[convert(id, 2)] += Bpixels[i];
			count_pixel_per_centroid[id]++;
			centroid_id_per_pixel[i] = id;
			//if (centroid_id_per_pixel[i] > 2 || centroid_id_per_pixel[i] < 0)cout << "7a7a" << endl;
		}

		MPI_Reduce(sum_of_distances_per_centroid, cluster_inner_sum, 3*K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(count_pixel_per_centroid, cluster_size, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rnk == 0)
		{
			int same = 1;
			for (int i = 0; i < K; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int k = convert(i, j);
					//cout << cluster_inner_sum[k] << " " << cluster_size[i] << endl;
					new_Centroid[k] = cluster_inner_sum[k] / (double)(cluster_size[i]);
					new_Centroid[k] = min(255, max(0, new_Centroid[k]));
					if (abs(new_Centroid[k] - Centroid[k]) > eps)same = 0;
					Centroid[k] = new_Centroid[k];
				}
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
			MPI_Recv(&same, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
			if (same)break;
		}


	}


	int* cluster = new int[ImageHeight * ImageWidth];

	MPI_Gatherv(centroid_id_per_pixel, pixel_cnt, MPI_INT, cluster, send_counts_scatterv, displs_scatterv, MPI_INT, 0, MPI_COMM_WORLD);
	if (rnk == 0)
	{
		for (int i = 0; i < ImageHeight * ImageWidth; i++)
		{
			R[i] = Centroid[convert(cluster[i],0)];
			G[i] = Centroid[convert(cluster[i], 1)];
			B[i] = Centroid[convert(cluster[i], 2)];
		}
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		createImage(R,G,B, ImageWidth, ImageHeight, 1);
		cout << "time: " << TotalTime << endl;
	}
	free(imageData);
	MPI_Finalize();
	return 0;
}