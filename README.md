## Dissertation Topic: Improve UCR Suite by Lower Resolution Techniques

*Work in progress. To be completed in May 2020.*

UCR Suite written in Python.
The original UCR Suite is written in C++ (using pointers). I translated it to Python (using linked lists).
UCR Suite: https://www.cs.ucr.edu/~eamonn/UCRsuite.html

- *Note 1: The Python script works. But (1) LB_Keogh2 of original UCR Suite is not included here, and (2) the function **distance.dynamicTimeWraping( )** is an approximation of DTW only. The true DTW, which use dynamic programming, is **distance.dynamicTimeWraping_true( )**. It is very slow since Python is not designed for dynamic programming.*

- *Note 2: Library Matplotlib is needed if you wish to visualize the code. No external library is needed if you just want to find k-nearest neighbours.*

I tested searching a query with length = 421 and *Sakoe-Chiba band* = 5 among 0.5 million of data. {electrocardiogram data which is used in the original paper https://www.cs.ucr.edu/~eamonn/SIGKDD_trillion.pdf. I extracted the first 0.5 million of data}

**Experiment Result:**

Query length = 421
- *Need to run LB_Kim : 100 %*
- *Need to run LB_Keogh : 92.36 %*
- *Need to run Dynamic Time Wraping: 0.63 %*

As you can see, *LB_Kim* only filter out 7.64 % of subsequences. It is not useful for long query.
In the next step, I will try to improve the speed of UCR Suite by using lower resolution techniques and replacing *LB_Kim*.
It is expected to be completed before March 2020.

*Regarding to Note 1, please note that the purpose of this project is to design a new lower bound to replace LB_Kim. True DTW and LB_Keogh2 are not the major concerns in this project. I will update LB_Keogh2, but not intend to rewrite True DTW to make it more efficent.*

-----------------------------------------------------------------------------------------------------

Further updates here:
https://github.com/ngyiuwai/PolyU-COMP5940-Final-Dissertation-Improve-UCR-Suite-by-Low-Resolution-Technique

A new lower bound function, *LB_LowResED* is created. It accelerates the calculation of Euclidean distance by low resolution technique. *[Note: LB_Keogh is Euclidean distance for upper/lower bound of query & subsequences.]*

*LB_LowResED* is a better choice comparing with *LB_Kim* for long query. It is a new cascading lower bound, which can improve UCR Suite.
