#include<stdio.h>
#include<math.h>

#include"volume.h"

struct segment
{
	int begin;
	int stop;
	int duration;
	struct segment *prv;
	struct segment *next;
};

//#define frameSize 512
int frameSize = 512;
int overlap = 0;
int frameCount = 0;
int minSegment = 0;
double dminSegment = 0.05;
int volRatio = 10;

double **out;
double *volume;

double volMin;
double volMax;
double volTh;

int index=0;

double q[513];

struct segment *sp;
struct segment *current_sp = 0;

int myround(double number)
{
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

void swap(double *x, double *y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

double mean(double *data, int size)
{
        int i=0;
        double sum=0;

        for (i=1;i<=size;i++)
        {
                sum = sum + data[i] ;
        }

        return (double) sum/(double) size;
}

void QuickSort(double *A, int p, int r)
{
    int q;
    if(p < r)
    {
        q = Pation(A, p, r);
        QuickSort(A, p, q-1);
        QuickSort(A, q+1, r);
    }
}

int Pation(double *A, int p, int r )
{
    int x;
    int i,j;
    x = A[r];
    i = p-1;
    for(j=p; j < r; j++)
    {
        if(A[j] <= x)
        {
            i += 1;
            swap(&A[i], &A[j]);
        }       
    }
    swap(&A[i+1], &A[r]);
    return i+1;
}

double median(double *ptr, int size)
{
        int i=0;
        int middle = size/2 + 1;
        double average;

        //clone
        for (i=1; i<=frameSize; i++)
	 {
                q[i] = ptr[i];
	 }

        QuickSort(q, 1, frameSize);

        if (size % 2 ==0)
	{
                average = q[middle-1] + q[middle]/2;
	}
        else    
	{
                average = q[middle];
	}
        
        return average;
}

void buffer2(double *wave, int frameSize, int overlap, int wave_size)
{
	int i=0, j=0, k=0;
	int step = frameSize-overlap;
	int startIndex;

	frameCount = floor((wave_size-overlap)/step);
	printf("step=%d, frameCount=%d\n", step, frameCount);
	//create matrix
	out = (double **) malloc((frameCount+1)*sizeof(double *)); 

	for (i=1; i<=frameCount; i++)
	{
		out[i] = (double*) malloc((frameCount+1)*(frameSize+1)*sizeof(double));
		startIndex = (i-1)*step+1;
		for (k=1,j=startIndex; j<=startIndex+frameSize-1; j++)
		{
			out[i][k] = wave[j];
			k++;
		}

		if (k != frameSize + 1)
		{
			printf("fill zero %d, %d\n", k, frameSize);
			//fill zero
			for (j=k; j<=frameSize; j++)
				out[i][j] = 0;
		}
	}
}

void frame2volume(double **frameMat, int usePolyfit)
{

	int i=0, j=0;
	double frame[513];
	double sum=0;
	double med;

	volume = (double *)malloc((frameCount+1)*sizeof(double *));

	for (i=1; i<=frameCount; i++)
	{
		//med = median(frameMat[i], frameSize);
		for (j=1; j<=frameSize; j++)
		{
			frame[j] = frameMat[i][j] - med;
			volume[i] = volume[i] + fabs(frame[j]);
		}
	}

	//free(frame);
}

void segmentFind(double *v, int size, double volTh)
{
	int i=0;
	int start_flag=0;
	int duration = 0;
	struct segment *sg;

	printf("segmentFind\n");
	for (i=1; i<=size; i++)
	{
		if (v[i] > volTh && !start_flag)
		{
			start_flag=i;
		}
		else if (v[i] <= volTh && start_flag != 0)
		{
			duration = i - start_flag + 1;
			//skip
			if (duration <= minSegment) 
			{
				start_flag=0;
				continue;
			}
			else
			{
			     if (!current_sp)
			     {
				sp = (struct segment *) malloc(1*sizeof(struct segment));
				sp->prv=0;
				current_sp = sp;
				sp->next=0;
			     }
			     else
			     {
				sg = current_sp;
			
				//create new node
			       current_sp->next = 
 				    (struct segment *) malloc(1*sizeof(struct segment));
			   	current_sp = current_sp->next;
				current_sp->prv = sg;
				current_sp->next = 0;
			     }	
			     printf("begin %d, stop %d, duration %d\n", 
							start_flag, i-1 , duration);

			     current_sp->begin = start_flag;
			     current_sp->stop = i-1;
			     current_sp->duration = duration;
			     start_flag=0;
			}				
		}
	}	
}

double getVolumVolTh(double *v ,int size)
{
	double q[513];
	int i=1;

	for (i=1; i<=frameSize; i++)
	{
		q[i] = v[i];
	}

	QuickSort(q, 1, frameSize);
	
	 index = myround(frameCount/32);
        volMin = q[index];
        volMax = q[frameSize-index+1];

	printf("volMin=%f, volMax=%f\n", volMin, volMax);

	return (volMax-volMin)/volRatio+volMin;
}


void epdByVol(int fs)
{
	buffer2(voice, frameSize,  overlap, 16896);	
	minSegment= myround( dminSegment*fs / (frameSize-overlap) );
	frame2volume(out, 0);

	volTh = getVolumVolTh(volume, frameCount);
	printf("minSegment = %d, volTh = %f\n", minSegment, volTh);
	segmentFind(volume, frameCount, volTh);
}

int main(void)
{
	int i, j;
	struct segment *free_sp;
	epdByVol(8000);

	for (; sp->next != 0;)
        {
	    for(i=sp->begin; i<=sp->stop; i++)
	    {
		for (j=1; j<= frameSize; j++)
		{
			printf("%5.5f\n", out[i][j]);			
		}
	    }
	     free_sp = sp;
            sp = sp->next;
	     free(free_sp);
        }
	//frame: free
	for (i=1; i<=frameCount; i++)
        {
                free(out[i]);
        }
     
       free(out);
	free(volume);

	printf("finish.\n");
	return 0;
}
