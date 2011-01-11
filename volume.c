#include<stdio.h>
#include<math.h>

double voice_data = {0.0}

int frameSize = 256;
int overlap = 128;

double **out;
double *volume;

void buffer2(double *wave, int frameSize, int overlap, int wave_size)
{
	int i=0, j=0, k=0;
	int step = frameSize-overlap;
	int frameCount = floor((wave_size-overlap)/step);

	//create matrix
	out = (double **)malloc((frameCount+1)*sizeof(double *) +  
					(frameCount+1)*(frameSize)*sizeof(double)); 

	for (i=1; i<=frameCount; i++)
	{
		int startIndex = (i-1)*step+1;
		k=0;
		for (j=	startIndex; j<=startIndex+frameSize-1; j++)
		{
			out[i] = subout;
			out[i][k] = wave[j];
			k++;
		}

		if (k-1 != frameSize)
		{
			//file zero
			for (j=k; j<=frameSize; j++)
				out[i][j] = 0;
		}
	}
}

void frame2volume(double **frameMat, int frameNum, int frameSize, int usePolyfit)
{
	int i=0, j=0;
	double *frame = (double *)malloc((frameSize)*sizeof(double *));
	double sum=0;

	volume = (double *)malloc((frameNum)*sizeof(double *));

	for (i=1; i<=frameNum; i++)
	{
		for (j=1; j<=frameSize; j++)
		{
			frame[j]=frameMat[i][j] - median(frame);
			volume[i] = volume[i] + fabs(frame[j]);
		}
	}
	

/*
volume=zeros(1, frameNum);
for i=1:frameNum
	frame=frameMat(:,i);
	frame=frame-median(frame);
	volume(i)=sum(abs(frame));
end
*/
}


double mean(double *data, int size)
{
	int i=0;
	double sum=0;
	for (i=0;i<size;i++)
	{
		sum = sum + data[i] ;
	}
	return sum/size;
}


void main(void)
{
		

}
