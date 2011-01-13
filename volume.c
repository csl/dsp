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

int frameSize = 512;
int overlap = 0;
int frameCount = 0;
int minSegment = 0;
double dminSegment = 0.05;
int volRatio = 10;

double **out;
double *out_sub;
double *volume;

double volMin;
double volMax;
double volTh;

int index=0;

struct segment *sp;
struct segment *current_sp = 0;

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

	printf("%10.10f\n",sum);

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
        double *q = (double *)malloc((size+1)*sizeof(double *));
        int i=0;
        int middle = size/2 + 1;
        double average;

        //clone
        for (i=1; i<=size; i++)
                q[i] = ptr[i];

        QuickSort(q, 1, size);

	 //record
	 volMin = q[index];
	 volMax = q[size-index+1];

        if (size % 2 ==0)
	{
                average = q[middle-1] + q[middle]/2;
	}
        else    
	{
                average = q[middle];
	}

        free(q);
        
        return average;
}

void buffer2(double *wave, int frameSize, int overlap, int wave_size)
{
	int i=0, j=0, k=0;
	int step = frameSize-overlap;
	frameCount = floor((wave_size-overlap)/step);
	double **out_sub;

	printf("step=%d, frameCount=%d\n", step, frameCount);

	//create matrix
	out = (double **) malloc((frameCount+1)*sizeof(double *)); 

	for (i=1; i<=frameCount; i++)
	{
		int startIndex = (i-1)*step+1;
		out[i] = (double*) malloc((frameCount+1)*(frameSize)*sizeof(double));
		k=1;
		for (j=	startIndex; j<=startIndex+frameSize-1; j++)
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
	double *frame = (double *)malloc((frameSize+1)*sizeof(double *));
	double sum=0;

	volume = (double *)malloc((frameCount+1)*sizeof(double *));

	for (i=1; i<=frameCount; i++)
	{
		double med = median(frameMat[i], frameSize);

		for (j=1; j<=frameSize; j++)
		{
			frame[j] = frameMat[i][j] - med;
			volume[i] = volume[i] + fabs(frame[j]);
		}
	}
}

void segmentFind(double *v, int size, double volTh)
{
	int i=0;
	int start_flag=0;
	int duration = 0;

	sp = (struct segment *) malloc(1*sizeof(struct segment));
	current_sp = sp->next;

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
				current_sp = sp;
				sp->prv=0;
			     }
			     else
			     {
			        current_sp->next = 
 				    (struct segment *) malloc(1*sizeof(struct segment));
				struct segment *sg = current_sp;			
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
	double *q = (double *)malloc((size+1)*sizeof(double *));
	int i=1;

	for (i=1; i<=size; i++)
	{
		q[i] = v[i];
	}

	QuickSort(q, 1, size);
	
	index = round(frameCount/32);
        volMin = q[index];
        volMax = q[size-index+1];

	printf("volMin=%f, volMax=%f\n", volMin, volMax);

	return (volMax-volMin)/volRatio+volMin;
}


void epdByVol(int fs)
{
	buffer2(voice, frameSize,  overlap, 16896);	
	minSegment= round( dminSegment*fs / (frameSize-overlap) );
	frame2volume(out, 0);

	volTh = getVolumVolTh(volume, frameCount);
	printf("minSegment = %d, volTh = %f\n", minSegment, volTh);
	segmentFind(volume, frameCount, volTh);
}

int main(void)
{
	int i, j;
	struct segment *free_sp;

	FILE *fp = fopen("out.txt","w");

	epdByVol(8000);
	
	for (; sp->next != 0;)
        {
	    for(i=sp->begin; i<=sp->stop; i++)
	    {
		for (j=1; j<= frameSize; j++)
		{
			fprintf(fp, "%5.5f\n", out[i][j]);			
		}
	    }
            sp = sp->next;
        }

	fclose(fp);
	//free: segment
	for (; current_sp->prv != 0;)
        {
		free_sp = current_sp;
                current_sp = current_sp->prv;
		free(free_sp);
        }

	//frame: free
	for (i=1; i<=frameCount; i++)
        {
                free(out[i]);
        }
        free(out);
	free(volume);

	return 0;
}
