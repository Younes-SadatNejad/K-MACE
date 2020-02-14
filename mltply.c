//the following spatial filter takes an image and a kernel as an input and performs the convolution to obtain the output.
//another parameter is sent since the spatial filter can be used for edge detection and as an all purpose spatial filter
//spatila filter requires normalization so that the output is from 0 to 1. however, if this is used for edge detection, normalization
//is not applicable because we want - and + values to retain both the postive and the negative edges
//To choos, 3rd parameter shoudl be '1' if we want to use it as spatial edge detection filter otherwise it shoudl just be zero
#include "mex.h"
#include <math.h>
//this function does a non-linear contrast transform with the use of a Look up table

//the LUT is a 256 element array that defines that characteristic of the transform that will be applied



void mltply(int width_img, int height_img, const double *input,  double *output)
{
   int counter = 0;
   double temp;
   for (int i=0; i < height_img-1;i++)
       for (int j=i+1; j < height_img;j++)
       {    
           temp = 0;
           for (int d=0;d<width_img;d++)
           {
                output[counter + (d*height_img*(height_img-1)/2)] = input[i + (d*height_img)]*input[j + (d*height_img)];
           }
           counter++;
       }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    
    if (nlhs != 1)  //if the number of our outputs is not equal to 1, we should declared an error, we are expecting only one output which is the output image
        mexErrMsgTxt("mali");
    
    if (nrhs != 1)  //we are expecting only 1 input which is the array
        mexErrMsgTxt("mali");
    
    
    const mxArray *img = prhs[0];   //the image is the first input of the function call
   
    
    if (!mxIsDouble(img))
        mexErrMsgTxt("input should be of type 'double'");
    
    /*extracting the dimensions of the image*/
    mwSize ndims_img = mxGetNumberOfDimensions(img);
    const mwSize *dims_img = mxGetDimensions(img);
    
    int height_img  = dims_img[0]; //we figure out the dimension of our output image since 
    int width_img    = dims_img[1];
    
    //by knowing this, we can figure out now big the array that we need to allcoate to pass the image to another function
    mxClassID input_type_img = mxGetClassID(img);
     
    
    mxArray *output =  mxCreateDoubleMatrix(height_img*(height_img-1)/2, width_img, mxREAL);
   // mxArray *output = mxCreateNumericArray(1, height_img*(height_img-1), input_type_img, mxREAL);
    plhs[0] = output;
       
    
    double *img_ptr = ( double *)mxGetData(img);
    double *out_ptr = ( double *)mxGetData(output); 
    
    //call the operation, send in the required paramters
    mltply(width_img, height_img, img_ptr,out_ptr);
}
