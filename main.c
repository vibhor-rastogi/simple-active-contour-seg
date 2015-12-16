#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
FILE *fpt;
unsigned char *image;
unsigned char *imageCopy;
unsigned char *imageX;
unsigned char *imageY;
unsigned char *imageG;
signed char *sobelX;
signed char *sobelY;
float *map;
float *map1;
float *map2;
float *map3;
int *data;//store contour point data
//int *diffData;//store difference between consecutive points
char header[320];
char filename[30];
char filename1[30];
int    ROWS,COLS,BYTES;
int    r1,c1,r2,c2,i,j,k,sumX,sumY,x,x1,y,y1,distance,avgDistance;
int    min,max,min1,max1,min2,max2,min3,max3,temp,temp1,contourPoints=0;
int    nextPos;
static int windowSize=13; // set window size here

printf("Enter input image name: ");
gets(filename);
printf("Enter contour points file name: ");
gets(filename1);

printf("Reading input image...\n");
if ((fpt=fopen(filename,"rb")) == NULL)
{
  printf("Unable to open %s for reading\n",filename);
  exit(0);
}
fscanf(fpt,"%s %d %d %d\n",header,&COLS,&ROWS,&BYTES);
if (strcmp(header,"P5") != 0  ||  BYTES != 255)
{
  printf("Not a greyscale 8-bit PPM image\n");
  exit(0);
}

image=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
imageCopy=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
imageX=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
imageY=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
imageG=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
sobelX=(signed char *)calloc(3*3,sizeof(unsigned char));
sobelY=(signed char *)calloc(3*3,sizeof(unsigned char));
map=(float *)calloc(windowSize*windowSize,sizeof(float));
map1=(float *)calloc(windowSize*windowSize,sizeof(float));
map2=(float *)calloc(windowSize*windowSize,sizeof(float));
map3=(float *)calloc(windowSize*windowSize,sizeof(float));

sobelX[0]=-1;
sobelX[1]=0;
sobelX[2]=1;
sobelX[3]=-2;
sobelX[4]=0;
sobelX[5]=2;
sobelX[6]=-1;
sobelX[7]=0;
sobelX[8]=1;

sobelY[0]=-1;
sobelY[1]=-2;
sobelY[2]=-1;
sobelY[3]=0;
sobelY[4]=0;
sobelY[5]=0;
sobelY[6]=1;
sobelY[7]=2;
sobelY[8]=1;

header[0]=fgetc(fpt);
fread(image,1,COLS*ROWS,fpt);
fclose(fpt);
printf("Done.\n\n");

for(i=0;i<ROWS*COLS;i++)
    imageCopy[i]=image[i];


for (r1=1; r1<ROWS-1; r1++)
  for (c1=1; c1<COLS-1; c1++)
    {
        sumX=0;
        sumY=0;
        for (r2=-1; r2<=1; r2++)
            for (c2=-1; c2<=1; c2++)
            {
                temp=image[(r1+r2)*COLS+(c1+c2)];

                //convolve with sobelX
                sumX+=(temp*sobelX[(r2+1)*3+c2+1]);

                //convolve with sobelY
                sumY+=(temp*sobelY[(r2+1)*3+c2+1]);
            }
        imageX[r1*COLS+c1]=(int)(sumX/9);
        imageY[r1*COLS+c1]=(int)(sumY/9);
        imageG[r1*COLS+c1]=(int)(imageX[r1*COLS+c1]*imageX[r1*COLS+c1])+(imageY[r1*COLS+c1]*imageY[r1*COLS+c1]);
    }

printf("\nWriting sobelX to disk...\n");
fpt=fopen("sobelX.ppm","w");
fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
fwrite(imageX,COLS*ROWS,1,fpt);
fclose(fpt);

printf("\nWriting sobelY to disk...\n");
fpt=fopen("sobelY.ppm","w");
fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
fwrite(imageY,COLS*ROWS,1,fpt);
fclose(fpt);

printf("\nWriting sobelG to disk...\n");
fpt=fopen("sobelG.ppm","w");
fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
fwrite(imageG,COLS*ROWS,1,fpt);
fclose(fpt);

printf("Reading contour input file...\n");
if ((fpt=fopen(filename1,"r")) == NULL)
{
  printf("Unable to open %s for reading\n",filename1);
  exit(0);
}
while(fscanf(fpt,"%d %d\n",&x,&y)==2)
{
    image[y*COLS+x]=255;
    //form a 7x7 '+' at this position
    for(i=1;i<=3;i++)
    {
        image[y*COLS+x-i]=255;
        image[y*COLS+x+i]=255;
        image[(y-i)*COLS+x]=255;
        image[(y+i)*COLS+x]=255;
    }
    contourPoints++;
}
fclose(fpt);

data=(int *)calloc(2*contourPoints,sizeof(int));
//diffData=(int *)calloc(contourPoints,sizeof(int));

printf("Reading contour input file...\n");
fpt=fopen(filename1,"r");
i=0;
distance=0;
while(fscanf(fpt,"%d %d\n",&x,&y)==2)
{
    //store contour points in a matrix in memory
    data[i++]=x;
    data[i++]=y;

    //calculate average distance between contour points
    if(i==2) //leave first contour point
        continue;
    temp=abs(data[i-2]-data[i-4]);
    temp=temp*temp;

    temp1=abs(data[i-1]-data[i-3]);
    temp1=temp1*temp1;
    distance+=(sqrt(temp+temp1));
}
temp=abs(data[i-2]-data[0]);//last and first contour point
temp=temp*temp;
temp1=abs(data[i-1]-data[1]);
temp1=temp1*temp1;
distance+=(sqrt(temp+temp1));
avgDistance=distance/contourPoints;
fclose(fpt);


//iterations begin
//map1 - first internal energy parameter in the window
//map2 - second internal energy parameter
//map3 - external energy parameter
k=0;
while(k<3)//set no. of iterations here
{
printf("\n\n\n");
temp=0;//index for current contour point
temp1=2;//index for next contour point

//operating on energies
for(i=0;i<contourPoints;i++,temp=temp+2,temp1=(temp1+2)%(contourPoints*2))
{
    //x and y are coordinates of current contour point
    x=data[temp];
    y=data[temp+1];
    //printf("\nx and y are %d, %d",x,y);
    //x1 and y1 are coordinates of next contour point
    x1=data[temp1];
    y1=data[temp1+1];
    //printf("\nx1 and y1 are %d, %d",x1,y1);
    j=0;
    for (r2=-(windowSize/2); r2<=(windowSize/2); r2++)
        for (c2=-(windowSize/2); c2<=(windowSize/2); c2++)
        {
            //maps store energy data of windows

            //first internal energy parameter
            map1[j]=(x+c2-x1)*(x+c2-x1)+(y+r2-y1)*(y+r2-y1);

            //second internal energy parameter
            map2[j]=avgDistance - sqrt((x+c2-x1)*(x+c2-x1)+(y+r2-y1)*(y+r2-y1));
            map2[j]=map2[j]*map2[j];

            //external energy parameter
            map3[j]=imageG[(y+r2)*COLS+(x+c2)];

            j++;
        }

    //normalize data
    min1=map1[0];
    max1=map1[0];
    min2=map2[0];
    max2=map2[0];
    min3=map3[0];
    max3=map3[0];
    for(j=0;j<windowSize*windowSize;j++)
    {
        if(min1>map1[j])
            min1=map1[j];//find min
        if(max1<map1[j])
            max1=map1[j];//find max

        if(min2>map2[j])
            min2=map2[j];//find min
        if(max2<map2[j])
            max2=map2[j];//find max

        if(min3>map3[j])
            min3=map3[j];//find min
        if(max3<map3[j])
            max3=map3[j];//find max
    }
    for(j=0;j<windowSize*windowSize;j++)
    {
        //normalize
        map1[j]=((float)(map1[j]-min1))/((float)(max1-min1));
        map2[j]=((float)(map2[j]-min2))/((float)(max2-min2));
        map3[j]=((float)(map3[j]-min3))/((float)(max3-min3));
        map3[j]*=-1;//negate external energy
    }

    //sum energies
    for(j=0;j<windowSize*windowSize;j++)
        map[j]=map1[j]+map2[j]+map3[j];

    //find minimum energy and move current contour point
    min=map[0];
    nextPos=(windowSize/2)*(windowSize/2);
    for(j=0;j<windowSize*windowSize;j++)
    {
        if(min>map[j])
        {
            min=map[j];//find min
            nextPos=j;
        }
    }

    //update contour point
    //printf("\nBefore - data[temp] is %d, data[temp+1] is %d,nextPos is %d",data[temp],data[temp+1],nextPos);
    data[temp]=data[temp]+((nextPos%windowSize)-(windowSize/2));//x coordinate
    data[temp+1]=data[temp+1]+((nextPos/windowSize)-(windowSize/2));//y coordinate
    //printf("\nAfter - data[temp] is %d, data[temp+1] is %d",data[temp],data[temp+1]);

}
k++;
}//end of iterations loop


fpt=fopen("initial.ppm","w");
fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
fwrite(image,COLS*ROWS,1,fpt);
fclose(fpt);

printf("\nInitial image saved..\n");

fpt=fopen("finalPoints.txt","w");
temp=0;
for(i=0;i<contourPoints;i++)
    fprintf(fpt,"%d %d\n",data[temp++],data[temp++]);
fclose(fpt);
printf("\nFinal points text file generated and saved as finalPoints.txt...\n");


printf("\nReading contour output file...\n");
fpt=fopen("finalPoints.txt","r");
while(fscanf(fpt,"%d %d\n",&y,&x)==2)
{
    imageCopy[y*COLS+x]=255;
    //form a 7x7 '+' at this position
    for(i=1;i<=3;i++)
    {
        imageCopy[y*COLS+x-i]=255;
        imageCopy[y*COLS+x+i]=255;
        imageCopy[(y-i)*COLS+x]=255;
        imageCopy[(y+i)*COLS+x]=255;
    }
}
fclose(fpt);

fpt=fopen("finalImage.ppm","w");
fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
fwrite(imageCopy,COLS*ROWS,1,fpt);
fclose(fpt);
printf("\nFinal image saved to disk as finalImage.ppm...\n");

return 0;
}
