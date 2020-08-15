cd ./ID/
mex -O -largeArrayDims QRpartial.c -lmwlapack -lmwblas
mex -O -largeArrayDims QRdim_c.c -lmwlapack -lmwblas
%If there is error about C99, C11 mode: 
%mex -O -largeArrayDims QRpartial.c -lmwlapack -lmwblas CFLAGS='$CFLAGS -std=c99'
%mex -O -largeArrayDims QRdim_c.c -lmwlapack -lmwblas CFLAGS='$CFLAGS -std=c99'

cd ../
cd ./kernel_function/
mex -O -largeArrayDims rpy_mex.c
cd ../
