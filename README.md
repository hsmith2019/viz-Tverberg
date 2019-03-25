# viz-Tverberg
Visualizing Tverberg's theorem in the plane.

https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;&\text{\underline{Tverberg's&space;Thm:}}&space;\\&space;&\quad&space;\text{For&space;$&space;\ge&space;(d&plus;1)(k-1)&plus;1$&space;points&space;in&space;$\mathbb{R}^d$,&space;there&space;exists&space;a&space;partition&space;}&space;\\&space;&\quad&space;\text{into&space;$A_1,\dots,A_k$&space;such&space;that&space;$\textstyle&space;\bigcap_{i=1}^k&space;\mathrm{conv}(A_i)&space;\not=&space;\varnothing$.}&space;\end{align*}


What we plan to do:
Create program that takes as input n points within a rectangular grid, and returns a Tverberg partition of size k = floor((n-1)/(d+1)+1). This will be used to plot the convex hulls of each of the k partitioned sets in a different color, outlining their intersection in black.
