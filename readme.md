# Linear mathematical models solver in pure C
The repo was previously named RSA. To eliminate all sorts of confusion with the RSA encryption algorithm, I will just call it Solver for now.

TODO: 

<del>1- Add error handling.</del>

2- <del>Add infeasibility check while solving.</del>

3- <del>Implement proper memory de-allocation.</del>

4- <del>Might as well substitute .txt data input with .csv data instead. Maybe even work on better data parsing since it does more computational effort, in my opinion. </del>

5- <del> Add debug mode </del>.

6- <del>Using Gauss' pivoting technique for matrix inversion is an effort of  $O(n^3)$, can I do better? </del> <del> Update: After learning a bit more on this topic, it would appear that using this algorithm is suitable for this application</del>. 
<del>Update 6/12/2025: Just learned that I don't need the entire B matrix. Using the LU factorization technique. There is way too much to improve in this area. </del>

7- Given how portable this program is, maybe I should implement Branch-and-Bound to consider integer points in the space.

8- Implement model pre-processing cuts to reduce space.

9- Change parsing so that constraints are parsed before the objective function, making it less expensive for the CPU to request and fill memory upon transforming.

# How to contribute 

Test this with as many LP models please. So far I gave it nothing but LPs from linear algebra/Operation research text books and internet forums. Be sure to report any issues if necessary.

# About
In this code I attempted to write a mathematical solver that utilizes the Revised Simplex Algorithm by George B. Dantzig (1953) <del>using the Big-M method</del>. The model input should be written (In general form) as such:

```
OBJECTIVE, Variables count, Constraints count
Coeffs
Constraints
```
Where `OBJECTIVE` is the objective function's direction which takes either `MINIMIZE` and `MAXIMIZE` as keywords. `Variables count` and `Constraints count` describe the number of variables and constraints in the mathematical formulation respectively. Note: An error will be thrown if you give the wrong `Variables count`, but if you give a lower value $m$ for the `Constraints count` and include more constraints then only the first $m$ constraints will be considered . `Coeffs` refers to the objective function coefficients for the variables all in one line, the next lines will be strictly for the constraints where each constraint will have its left hand sign contain nothing but variable coefficients preceding the constraint's symbol `<=`, `=` or `>=`. Make sure each entry is separated with a `,`. 
One last thing: Make sure <del>the right-hand side is positive and</del> that there is no empty line beneath the last constraint. Use 0s in the LHS of each constraints if a variable is irrelevant for the latter. This is because the program checks if each LHS has exactly $n$ number of variables, else the solver won't start.

### Example 1: 

```
MAXIMIZE,2,3 
9,7
10, 5, <=, 50
6, 6, <=, 36
4.5, 18, <=, 81
```
Which corresponds to:

$$
\begin{aligned}
\text{MAXIMIZE } & z = 9x_1 + 7x_2 \\
\text{subject to } 
& 10x_1 + 5x_2  \leq 50 \\
& 6x_2 + 6x_3 \leq 36 \\
& 4.5x_1 + 18x_2 \leq 81 \\
& x_1, x_2 \ge 0
\end{aligned}
$$

### Example 2:
```
MINIMIZE,3,3
-3,1,1
1,-2,1,<=,11
-4,1,2,>=,3
2,0,1,=,1
```
Which corresponds to:

$$
\begin{aligned}
\text{MINIMIZE } \quad 
& z = -3x_1 + x_2 + x_3 \\
\text{subject to } \quad
& x_1 - 2x_2 + x_3 \le 11 \\
& -4x_1 + x_2 + 2x_3 \ge 3 \\
& 2x_1 + x_3 = 1 \\
& x_1, x_2, x_3 \ge 0
\end{aligned}
$$


Support for bounded variable simplex mode is available. To enable it, simply add a new parameter in the input file right after the `constraints count`. After inputting the constraints, in a new line add the keyword `BOUNDS`, `bounded variables count`, where the latter should be identical to the number stated in the header. Each line is associated to a bounded variable and is interpreted as such: `Variable index, Lower bound, Upper bound`. The `Variable index` is 1-indexed. 

### Example 3:

```
MINIMIZE,5,5,5
4,6,3,5,2
2,1,1,0,3,>=,18
1,3,2,1,0,<=,25
0,2,1,4,1,>=,20
3,1,0,2,2,<=,30
1,1,1,1,1,=,15
BOUNDS,5
1,0,8
2,1,7
3,0,6
4,2,9
5,1,5
```

Which corresponds to: 

$$
\begin{aligned}
\text{Minimize } \quad 
z \;=\;& 4x_1 + 6x_2 + 3x_3 + 5x_4 + 2x_5 \\
\text{subject to} \quad 
2x_1 + x_2 + x_3 + 3x_5 \;&\geq\; 18 \\
x_1 + 3x_2 + 2x_3 + x_4 \;&\leq\; 25 \\
2x_2 + x_3 + 4x_4 + x_5 \;&\geq\; 20 \\
3x_1 + x_2 + 2x_4 + 2x_5 \;&\leq\; 30 \\
x_1 + x_2 + x_3 + x_4 + x_5 \;&=\; 15 \\
0 \leq x_1 \leq 8, \quad 1 \leq x_2 \leq 7, \;&\quad 0 \leq x_3 \leq 6 \\
2 \leq x_4 \leq 9, \quad 1 \leq x_5 \leq 5 \;&
\end{aligned}
$$




# Usage:
To compile it run this in the parent folder:
```
make
```

Implicitly, all variables are non-negative (of course) you won't need to consider this. Explicitly adding non-negativity domain definitions for each variable will still allow the program to work and output the correct answers but it will result in having extra memory usage and more runtime (Similarly, if you use `>=` type constraint, you will introduce some extra memory and time costs, please consider reformulating). The algorithm uses the double floating precision. In the file `RSA.c` I defined some macros to limit the number of variables and constraints in the input file, feel free to adjust these for bigger models. To run this program, simply pass the text file and objective arguments:

```./RSA "filepath" "-Debug"```


Where `filepath`is the `.csv` file path, `-Debug` is an optional flag that can be added as an argument, allowing you to see the solver operations step by step, the displayed indices are 0 indexed. The program will convert the problem to its canonical form then iteratively execute the algorithm until it terminates. The reason I am using the .csv file format is because of how portable it is + you can easily view whether the entries are valid using a Graphical CSV reader (e.g. OnlyOffice) to easily display whether data is missing.


*Note*: By default the solver runs the Two-phase method scheme, which is slower since it solves twice. If you wish to use the Big M method (Artificials get max(coeffs) * 2 as a value), comment-out the function Solve_BigM and comment Solve in lines `86-87` before compiling. The Debug mode will run using the Big M method. The regular solving will use Two-Pass. I have yet to implement an additional flag to switch between the two modes for the sake of consistency. 



