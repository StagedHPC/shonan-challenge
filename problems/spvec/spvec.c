Sparse vectors

Data structure for sparse matrices is an important problem in high performance 
computing.  Here, I reduce the problem from 2-dimensional array to 1 dimension, 
that is, a sparse vector.

A standard dense vector is represented by an array (double*) to store the 
elements, and an integer (in int) for the length.  The following routine is called
daxpy, because it computes "a x + y" ('p' stands for 'plus') in double (here 'd'
comes).

void daxpy(int n, double *z, double a, double *x, double *y) {
  int i;
  for (i=0; i< n; i++)
    z[i] = a * x[i] + y[i];
}

Performance improvements are expected if we generate a specialized routine for
the cases without alias in the arguments.  In addition, we could generate a
specialized code for a == 0.  It is easy for application programmers to write
a code like

if (a == 0)
  memcpy(z, y, sizeof(double) * n);
else
  daxpy(n, z, a, x, y);

It is not good to write such a code here and there, so one has to write a function 
to do that.  If the programming system provides a simple mechanism for such a 
specialization, it could be useful.  There are some questions about this 
specialization.  Is it allowed to pass a NULL pointer to x when a is zero?  What we 
should do if a is zero and x contains values such as NaN and Inf? 

Next sparse vector is introduced.  The following data structure is one possible representation of sparse vectors.

typedef spvec {
  int n;		// the length of the vector
  int nnz;	// the number of the non-zero elements
  int *idx;	// the indices of the non-zero elements
  double *val;	// the values of the non-zero elements
} spvec;

The semantics can be specified by the following code which convert sparse 
vector into a dense one.

void spv2vec(double *u, spvec v) {
  int i;
  for (i=0; i< v.n; i++)
    u[i] = 0.0;

  for (i=0; i< v.nnz; i++)
    u[v.idx[i]] = v.val[i];
}

In this routine, the indices of the input vector v are not required to be 
sorted.  But probably the performance is better if the indices are sorted, 
because higher locality is expected and the number of cache mishits could be
reduced.  

We should assume that there is no two elements of v.idx are the same.

If we assume that the indices are sorted, then the following code works, which 
removes many zero-clears, especially in case of relatively dense vectors.

void spv2vec_sorted(double *u, spvec v) {
  int i, j=0;
  for (i=0; i< v.nnz; i++) {
    while (j < v.idx[i])
      u[j++] = 0.0;
    u[j] = v.val[i];
  }
  while (j < v.n)
    u[j++] = 0.0;
}

However, the previous code (spv2vec) can be faster on actual machines than this code 
(spv2vec_sorted).  The two loops in spv2vec can be parallelized, and compilers may 
notice it and parallelize those loops.  Another point is that the inner loop of 
spv2vec_sorted can be short and may incur more pipeline stall and lower precision 
of hardware branch prediction.  So we do not know which (of spv2vec or spv2vec_sorted) 
is faster on specific hardware. 

1) Some daxpy routines with sparse vectors.

By using the above routines, we can introduce a sparse vector into one of the 
parameters, say, y.

void daxpspy(int n, double *z, double a, double *x, spv y) {
  assert(n == y.n);
  double u[y.n];
  spv2vec(u, y);
  daxpy(n, z, a, x, u);
}

The definition could be derived from overloading and automatic type conversion.
The question is how to optimize the above routine.  Or more precisely, how to
make it assure that it is optimized so that the temporal array u is removed.
If the temporal u is not removed, we HPC people will write down codes without
temporal variable by ourselves.

Similar routine is defined by changing the parameter x into a sparse vector.

void daspxpy(int n, double *z, double a, spv x, double *y) {
  assert(n == x.n);
  double u[x.n];
  spv2vec(u, x);
  daxpy(n, z, a, u, y);
}

The difference of this routine from the previous one is that for the specialized cases 
that z and y are the same vector, the computational complexity is O(nnz), where nnz 
is the number of the non-zero elements in x.

And another routine is defined by changing both input parameters into sparse
data structure.

void daspxpspy(int n, double *z, double a, spv x, spv y) {
  assert(n == y.n && n == x.n);
  double u[n], v[n];
  spv2vec(u, x);
  spv2vec(v, y);
  daxpy(n, z, a, u, v);
}

In this case, we can remove temporal array allocation, as in the previous
cases.

There are many other possibilities of data structure for sparse vectors than the
previously defined.  First, if the non-zero elements appear in groups of contiguous 
elements, then the sparse vector can be defined as a set of such groups of elements.

struct {
  int n; // length
  int ng; // number of groups
  int *size_g;  // sizes of groups
  int *iidx_g:  // initial index of groups
  double **val;  // values
} spv_g;

void spv_g2vec(double *z, spv_g x) {
  int i;
  for (i=0; i< x.n; i++)
    z[i] = 0;

  for (i=0; i< x.ng; i++) {
    int j;
    for (j=0; j< x.size_g[i]; j++)
      z[j + x.iidx_g[i]] = x.val[i][j];
  }
}

The advantage of this data structure is less memory requirements for index.
In the previous data structure, x.idx is an array of nnz elements, but here the length
of x.iidx_g is the number of groups.

Second, if the sizes of the groups are (multiples of) the same integer 
BLOCKSIZE, then the group size can be removed from the data structure.

struct {
  int n; // length
  int ng; // number of groups
  int *iidx_g:  // initial index of groups
  double **val;  // values
} spv_b;

void spv_b2vec(double *z, spv_g x) {
  int i;
  for (i=0; i< x.n; i++)
    z[i] = 0;

  for (i=0; i< x.ng; i++) {
    int j;
    for (j=0; j< BLOCKSIZE; j++)
      z[j + x.iidx_g[i]] = x.val[i][j];
  }
}

The advantage of this structure is that the inner loop (for j) can be
fully unrolled, and thus the routine spv_b2vec has no nested loops.

In these data structured, we HPC people may introduce some zeros in x.val[i][j],
if we find it to give a better performance.

Third, those data structures can be combined, since any vector can be
decomposed into sum of a few vectors:
v = v_b + v_s
where v_b is represented by spv_b and v_s is represented by spvec.
By such a combination, the performance of spv_b can be utilized, while the generality
of spvec is not lost.

2) Conversion from dense vector to sparse vector

To make a sparse vector from a dense vector, we need to count the number
of the non-zero elements.

spv vec2spv(int n, double *y) {
  int i, nnz=0;
  for (i=0; i< n; i++)
    if (y[i] != 0)
      n ++;

  spv x;
  x.n = n;
  x.nnz = nnz;
  x.idx = (int*) malloc(sizeof(int) * nnz);
  x.val = (double*) malloc(sizeof(double) * nnz);

  int k=0;
  for (i=0; i< n; i++)
    if (y[i] != 0) {
      x.idx[k] = i;
      x.val[k] = y[i];
      k++;
    }
}

Note that the above routine gives indices in an ascending order.

The next routine computes saxpy of sparse vectors resulting in a sparse vector.

spv spdaxpy(double a, spv x, spv y) {
  assert(x.n == y.n);
  double xx[x.n], yy[y.n], zz[x.n];
  spv2vec(xx, x);
  spv2vec(yy, y);
  daxpy(zz, a, xx, yy);
  return vec2spy(x.n, zz);
}

It is not trivial whether we can remove temporal data structure or not.  If
we are sure that the indices are sorted, then it is easy to count the number 
of non-zero elements of zz with O(1) storage in a way similar to the merge sort.

Otherwise: 
If the input vectors are nearly dense, then we can use temporal dense data structure
without sacrificing complexity.  But if the input vectors are very sparse, we could use
another algorithm.  Since clearly z.nnz <= x.nnz + y.nnz, z can be constructed in time 
O(x.nnz + y.nnz).

The following is an inner product of sparse vectors.

double spdot(double a, spv x, spv y) {
  assert(x.n == y.n);
  double xx[x.n], yy[y.n];
  spv2vec(xx, x);
  spv2vec(yy, y);
  return dot(x.n, xx, yy);
}

double dot(int n, double *x, double *y) {
  double d = 0.0;
  int i;
  for (i=0; i< n; i++)
    d += x[i] * y[i];
  return d;
}

Again, the optimum routine will be different when the indices are sorted or not.

How about elementwise multiplication of sparse vectors?

spv spvemul(double a, spv x, spv y) {
  assert(x.n == y.n);
  double xx[x.n], yy[y.n], zz[x.n];
  spv2vec(xx, x);
  spv2vec(yy, y);
  vemul(x.n, zz, xx, yy);
  return vec2spv(zz);
}

double vemul(int n, double *z, double *x, double *y) {
  int i;
  for (i=0; i< n; i++)
    z[i] = x[i] * y[i];
}

Again, the solution will be different when the indices are sorted or not.

Obviously, it is a good choice to restrict the indices to be sorted.
That makes many routines simpler.

3) Can we derive sparse vector data structure?

Consider the following code.

void addxtoy(int n, double *y, double *x) {
  int i;
  for (i=0; i< n; i++)
    if (x[i] != 0)
      y[i] += x[i];
}

This code alludes that it is enough to store only non-zero elements of x.
Also it says that we need indices and values.

If the compiler (language system) is asked to optimize the above routine with
a hint that many of x are zeros, how can it do that?