TODO: 

<del>1- Add error handling.</del>

2- <del>Add infeasibility check while solving.</del>

3- <del>Implement proper memory de-allocation.</del>

4- <del>Might as well substitute .txt data input with .csv data instead. Maybe even work on better data parsing since it does more computational effort, in my opinion. </del>

5- <del> Add debug mode </del>.

6- <del>Using Gauss' pivoting technique for matrix inversion is an effort of  $O(n^3)$, can I do better? </del> <del> Update: After learning a bit more on this topic, it would appear that using this algorithm is suitable for this application</del>. 
<del>Update 6/12/2025: Just learned that I don't need the entire B matrix. Using the LU factorization technique. There is way too much to improve in this area. </del>

7-Compress memory usage, make execution more efficient.

8- Given how portable this program is, maybe I should implement Branch-and-Bound to consider integer points in the space.

9- Implement model pre-processing cuts to reduce space.

10- Change parsing so that constraints are parsed before the objective function, making it less expensive for the CPU to request and fill memory upon transforming.

11- Do something about storing nothing but non-zero values to save up memory. 


