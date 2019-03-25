# viz-Tverberg
Visualizing Tverberg's theorem in the plane.

<img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;&\text{\underline{Tverberg's&space;Thm:}}&space;\\&space;&\quad&space;\text{For&space;$&space;\ge&space;(d&plus;1)(k-1)&plus;1$&space;points&space;in&space;$\mathbb{R}^d$,&space;there&space;exists&space;a&space;partition&space;}&space;\\&space;&\quad&space;\text{into&space;$A_1,\dots,A_k$&space;such&space;that&space;$\textstyle&space;\bigcap_{i=1}^k&space;\mathrm{conv}(A_i)&space;\not=&space;\varnothing$.}&space;\end{align*}" title="\begin{align*} &\text{\underline{Tverberg's Thm:}} \\ &\quad \text{For $ \ge (d+1)(k-1)+1$ points in $\mathbb{R}^d$, there exists a partition } \\ &\quad \text{into $A_1,\dots,A_k$ such that $\textstyle \bigcap_{i=1}^k \mathrm{conv}(A_i) \not= \varnothing$.} \end{align*}" />

What we plan to do:
Create program that takes as input n points within a rectangular grid, and returns a Tverberg partition of size k = floor((n-1)/(d+1)+1). This will be used to plot the convex hulls of each of the k partitioned sets in a different color, outlining their intersection in black.

https://doi.org/10.1112/jlms/s1-41.1.123
