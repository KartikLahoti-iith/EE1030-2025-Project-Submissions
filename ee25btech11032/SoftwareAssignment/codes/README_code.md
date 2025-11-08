# Execution Guidlines

## The following explains 

1.How to install and use **stb_image.h** and **stb_image_write.h** for image processing.
2.How to compile and run the **rSVD** and **Jacobi** algorithms to compress image according to requirements

### Required Tools
- **GCC Compiler**
- **C IDE** or **terminal**

### Required Libraries
1. "stb_image.h"
```https://github.com/nothings/stb/blob/master/stb_image.h```
2. "stb_image_write.h"
```https://github.com/nothings/stb/blob/master/stb_image_write.h```
3. Standard C Libraries:
 - "stdio.h"
 - "stdlib.h"
 - "math.h"

## Installing stb_image and stb_image_write

### Go to your project folder
cd path/to/your/project

### Download stb_image.h
```wget https://raw.githubusercontent.com/nothings/stb/master/stb_image.h```

### Download stb_image_write.h
```wget https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h```

## To Compress GreyScale or Coloured Image to GreyScale Using rSVD

In the codes/c_main/RunRsvd.
    -On line 59 provide the path to desired image
    -On line 94 choose the desired k value
    -On line 137 provide the path to save compressed image

## To Compress GreyScale or Coloured Image to GreyScale Using rSVD

In the codes/c_main/RunRsvd.
    -On line 59 provide the path to desired image
    -On line 94 choose the desired k value
    -On line 137 provide the path to save compressed image
    -```gcc RunRsvd.c ../c_libs/randomizedSvd.c -o RunRsvd -lm```
    -```./RunRsvd```

## To Compress Coloured Image Using rSVD

In the codes/c_main/ColorSvd.
    -On line 56 provide the path to desired image
    -On line 75 choose the desired k value
    -On line 149 provide the path to save compressed image
    -```gcc ColorSvd.c ../c_libs/randomizedSvd.c -o ColorSvd -lm```
    -```./ColorSvd```

## To Compress GreyScale or Coloured Image to GreyScale Using Only Jacobi

In the codes/c_main/Run.
    -On line 59 provide the path to desired image
    -On line 91 choose the desired k value
    -On line 163 provide the path to save compressed image
    -```gcc RunJacobi.c ../c_libs/oneSideJacobi.c -o RunJacobi -lm```
    -```./RunJacobi```



