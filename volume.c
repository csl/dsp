#include<stdio.h>
#include<stdlib.h>
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
const int frameSize = 512;

//volume
int overlap = 0;
int frameCount = 0;
int minSegment = 0;
long dminSegment = 0.05;
int volRatio = 10;
int volRatio2 = 5;

//zcr
double zcrTh=0;
double zcrRatio=0.1;
int zcrShiftGain=4;
struct segment segmentp[10];
int segment_count = 1;

//data
short out[34][513];
short outdata[513 * 34];

long outdiff[34][513];
long outhoddiff[34];
long VH[34];

long *volume;
short *zcr;

long volMin;
long volMax;
long volTh;

int index=0;

short myround(double number)
{
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

void buffer2(int frameSize, int overlap, int wave_size)
{
	int i=0, j=0, k=0;
	int step = frameSize-overlap;
	int startIndex;
	long all=0;
	int max=0, min=0;

	frameCount = floor((wave_size-overlap)/step);
	printf("step=%d, frameCount=%d\n", step, frameCount);
	//create matrix
	//out = (double **) malloc((frameCount+1)*sizeof(double *)); 

	for (i=1; i<=frameCount; i++)
	{
		//out[i] = (double*) malloc((frameCount+1)*(frameSize+1)*sizeof(double));
		startIndex = (i-1)*step+1;
		all=0;
		max=0;
		min=0;
		for (k=1,j=startIndex; j<=startIndex+frameSize-1; j++)
		{
			out[i][k] = voice[j];
			all = all + fabs(out[i][k]);
            		if (out[i][k] > max)  max = out[i][k];
            		if (out[i][k] < min) min = out[i][k];
			k++;			
		}
		
//		printf("startIndex=%d - %d, all %f\n",startIndex, startIndex+frameSize-1, all);
		if (k != frameSize + 1)
		{
			printf("fill zero %d, %d\n", k, frameSize);
			//fill zero
			for (j=k; j<=frameSize; j++)
				out[i][j] = 0;
		}
		/*
        	//normalize		
		printf("(max-min)/2 = %d\n", (max-min)/2);
		for (j=1; j<=frameSize; j++)
		{	
			if (max-min/2 != 0)
				out[i][j] = out[i][j] / ((max-min)/2);
		}
		*/
	}
}

void frame2volume(int usePolyfit)
{
	int i=0, j=0;
	long all = 0;

	volume = (long *)malloc((frameCount+1)*sizeof(long *));

	for (i=1; i<=frameCount; i++)
	{
		for (all=0, j=1; j<=frameSize; j++)
		{
			all = (double) all + (double) fabs(out[i][j]);
		}	
		volume[i] = all;
	}
}

void segmentFind(long *v, int size, double volTh)
{
	int i=0;
	int start_flag=0;
	int duration = 0;

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
			     printf("begin %d, stop %d, duration %d\n", 
							start_flag, i-1 , duration);

			     segmentp[segment_count].begin = start_flag;
			     segmentp[segment_count].stop = i-1;
			     segmentp[segment_count].duration = duration;
			     start_flag=0;
				 segment_count++;
			 }
		}
	}	
}

double getVolumVolTh(long *v ,int size)
{
	int i=1;
	double min=0, max=0;
	for (i=1; i<=frameCount; i++)
	{
		if (v[i] > max)
				max=v[i];
		if (v[i] < min)
				min=v[i];
	}

    index = myround(frameCount/32);
    volMin = min;
    volMax = max;

	return (volMax-volMin)/volRatio+volMin;
}


void epdByVol(int fs)
{
	buffer2(frameSize,  overlap, 15872);	
	minSegment= myround( dminSegment*fs / (frameSize-overlap) );
	frame2volume(0);

	volTh = getVolumVolTh(volume, frameCount);
	printf("minSegment = %d, volTh = %f\n", minSegment, volTh);
	segmentFind(volume, frameCount, volTh);
}

double mean(short *data, int size)
{
        int i=0;
        long sum=0;

        for (i=1;i<=size;i++)
        {
                sum = sum + data[i];
        }

        return (double)sum/size;
}

void frame2zcr(int method, long shiftAmount)
{
	int i=1, j=1;
	long frame[513];
	long all= 0;
	short total=0;

	zcr = (short *) malloc((frameCount+1)*sizeof(short));
/*
	//normalize
	for(i=1; i<= frameCount; i++)
	{
		for (j=1; j<=frameSize; j++)
		{
			fall = mean(out[i], frameSize);
			out[i][j] = out[i][j] - myround(fall);
		}	
	}
*/
	switch (method)
	{
		case 1:
			for(i=1; i<= frameCount; i++)
			{
				for (j=1; j<=frameSize; j++)
				{
					frame[j]=out[i][j] - shiftAmount;
				}

				total=0;
				for (j=1; j<=frameSize-1; j++)
				{
					all = frame[j] * frame[j+1];
					//condition
					if (all < 0)
					{
						total++;
					}
				}
				zcr[i] = total;
				//printf("zcr[%d]= %d, %d\n", i, total, zcr[i]);
			}

			break;
		case 2:
			for(i=1; i<= frameCount; i++)
			{
				for (j=1; j<=frameSize; j++)
				{
					frame[j]=out[i][j] - shiftAmount;
				}

				total=0;
				for (j=1; j<=frameSize-1; j++)
				{
					all = frame[j] * frame[j+1];
					//condition
					if (all <= 0)
					{
						total++;
					}
				}
				zcr[i] = total;
			}
			
			break;
	}
}

int getVectorMinIndex(long *t, int size)
{
	double min;
	int index = 1;
	int i=0;

	min = 0;
	for (i=1; i<=size;i++)
	{
		if (t[i] < min )
		{
			index = min;
		}
	}

	return index;
}


void epdByVolzcr(int fs)
{
	int index;
	short max = 0;
	int value = 0;
	int i = 1;
	long shiftAmount;

	index = getVectorMinIndex(volume, frameCount);
	printf("min index = %d\n", index);
	//fetch max
	max=0;
	for (i=1; i<=frameSize; i++)
	{
		value = fabs(out[index][i]);
		if (value > max)
		{
			max = value;
		}
	}

	shiftAmount = zcrShiftGain * max;
	printf("shiftAmount = %d\n", shiftAmount);
	frame2zcr(1, shiftAmount);

	//fetch max
	max=0;
	for (i=1; i<=frameCount; i++)
	{
		//printf("zcr[%d] = %d\n", i, zcr[i]);
		if (zcr[i] > max)
		{
			max = zcr[i];
		}
	}

	//get zcrTh
	zcrTh = max * zcrRatio;

	//printf("zcrTh, %f = %d * %f\n", zcrTh, max, zcrRatio);
}

int diff(int diffCount, int select)
{
	int i = 0, j =0;

	if (select == 0)
	{
		for (i=1; i<=frameCount; i++)
		{
			for (j=1; j<diffCount; j++)
			{
				outdiff[i][j] = out[i][j+1] - out[i][j];
			}	
		}
	}
	else
	{
		for (i=1; i<=frameCount; i++)
		{
			for (j=1; j<diffCount; j++)
			{
				outdiff[i][j] = outdiff[i][j+1] - outdiff[i][j];
			}	
		}
	}

	return diffCount-1;
}

int diffabs(int diffCount)
{
	int i = 0, j = 0;

	for (i=1; i<=frameCount; i++)
    {
    	for (j=1; j<=diffCount; j++)
        {
        	outdiff[i][j] = abs(outdiff[i][j]);
        }
    }

	return 1;
}

int diffsum(int diffCount)
{
    int i = 0, j = 0;
	long sum=0;

	for (i=1; i<=frameCount; i++)
	{
		sum=0;
		for (j=1; j<=diffCount; j++)
		{
           	sum = sum + outdiff[i][j];
		}	
		outhoddiff[i] = sum;
	}

    return 1;
}

void highOrderDiff(int fs)
{
	int max = 0, value = 0;
	int i,j;
	int diffCount = 0;

	buffer2(frameSize, overlap, 16896);
	minSegment= myround( dminSegment*fs / (frameSize-overlap) );
	frame2volume(0);

	diffCount = frameSize;  //3 x 5 
	diffCount = diff(diffCount, 0);

	diffCount = diff(diffCount, 1);
	diffCount = diff(diffCount, 1);
	diffCount = diff(diffCount, 1);
	diffabs(diffCount);

	diffCount = diffsum(diffCount);
}

int main(void)
{
	int i,j;
	int fsize = 0;
	int start,end;
	int head = 0, tail = 0;
	int max=0, min=0;
	int vhthr = 0;
	int startflag=0;
	double w = 0.5;
	double r = 0.125;

	//volume method
	epdByVol(8000);
/*
    for(i=1; i<segment_count; i++)
    {

		if (i==1)

		{

 		   start = segmentp[i].begin;

		}

		else if (i == segment_count-1)

		{

			end = segmentp[i].stop;

		}



	printf("begin %d, stop %d, duration %d\n",
			segmentp[i].begin, segmentp[i].stop ,segmentp[i].duration);



	}

	printf("begin %d, stop %d\n",
			start, end);



	for (j=start; j<= end; j++)
	{

		for (i=1; i<=frameSize; i++)

		{

			outdata[fsize] = out[j][i];

			fsize++;

		}
	}
*/
	epdByVolzcr(8000);
    for(i=1; i<segment_count; i++)
    {
		head = segmentp[i].begin;

		while ( (head-1)>=1 && (volume[head-1]>=volTh) )
		{
			head = head - 1;
		}
		segmentp[i].begin = head;

		tail = segmentp[i].stop;

		while ((tail+1)<=frameCount && volume[tail+1]>=volTh)
		{
        	tail = tail + 1;
		}

		segmentp[i].stop = tail;
		printf("begin %d, stop %d, duration %d\n",
			segmentp[i].begin, segmentp[i].stop ,segmentp[i].duration);
	}

	//% ====== Expansion 2: Expand end points to include high zcr region
    for(i=1; i<segment_count; i++)
    {
		head = segmentp[i].begin;

		while ( (head-1)>=1 && (zcr[head-1]>=zcrTh) )
		{
			head = head - 1;
		}
		segmentp[i].begin = head;

		tail = segmentp[i].stop;

		while ((tail+1)<=frameCount && zcr[tail+1]>=zcrTh)
		{
        	tail = tail + 1;
		}

		segmentp[i].stop = tail;
		printf("begin %d, stop %d, duration %d\n",
			segmentp[i].begin, segmentp[i].stop ,segmentp[i].duration);

	}

	highOrderDiff(8000);
	//VH = w*VOL + (1-w)*HOD
	//n = 4, w = 0.5, and r = 0.125
	for (i=1; i<=frameCount; i++)
	{
		VH[i] = w * volume[i] + (1-w) *  outhoddiff[i];
		if (VH[i] > max)
		{
			max = VH[i];
		}
		else if (VH[i] < min)
		{
			min = VH[i];
		}
	}

	//VHmin+(VHmax-VHmin)*r. 
	vhthr = min + (max - min) * r;
	printf("\nmax %d, min %d, vhthr = %d\n", max, min, vhthr);
	
	for (i=1; i<=frameCount; i++)
	{
		if (VH[i]>vhthr && startflag == 0)
		{
			startflag = 1;
			printf("start: %d\n", i);
		}
		else if (startflag==1 && VH[i]<vhthr)
		{
			startflag = 0;
			printf("end: %d\n", i);
		}
	}

	free(volume);
	//free(zcr);
	printf("volume finish.\n");
	return 0;
}
