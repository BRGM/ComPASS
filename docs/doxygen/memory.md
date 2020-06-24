Sparse matrix : CSR storage    {#memorypage}
===========================

# CSR type

## General presentation


The compressed sparse row (CSR) format represents a matrix
by three (one-dimensional) arrays, that respectively contain
nonzero values, the extents of rows, and column indices.
This format allows fast row access and matrix-vector
multiplications. <br>

Example of sparse matrix:  <br>
<center>
\f$
\left( \begin{array}{ccccc}
0 & 0 & 0 & 0 & 0 \\
1 & 2 & 0 & 0 & 0 \\
0 & 0 & 3 & 0 & 0 \\
0 & 0 & 4 & 0 & 5 \\
0 & 0 & 0 & 6 & 0 \end{array} \right)
\f$
</center> <br>


## Use in ComPASS code


In the ComPASS code, CSR format is mostly used to store continuous
informations (one-dimensional arrays) of different sizes,
ie the CSR format is used to store List Of List.
For example to store the numero of the nodes forming an
injection well. In this case, the number of line in the matrix correspond to
the number of injection wells and each line contains the nodes.
This is why it is not necessary to store the column
indices (continuous vectors).  <br>

The standard type CSR declared in the ComPASS code contains three main
informations (Nb, Pt and Num) and one optional vector (Val).

<pre class="fragment">
  type CSR
    integer :: Nb
    integer, allocatable, dimension(:) :: Pt
    integer, allocatable, dimension(:) :: Num
    integer, allocatable, dimension(:) :: Val
  end type CSR
</pre>

<ul>
<li> <b><code>Nb</code></b> represent the number of lines; </li>
<li> <b><code>Pt</code></b> is the index of the beginning of the rows.
By convention <code> Pt(1)=0 </code> then the nodes of the
injection well <code>i</code> are accessed by
<code>Pt(i)+1:Pt(i+1)</code>;     </li>
<li> <b><code>Num</code></b> contains the nonzero values of the matrix. </li>
<li> <b><code>Val</code></b> is a additionnal vector. In this declaration
it is an <code>integer</code> but it can also be
<code>double precision</code>. </li>
</ul>



Example of matrix stored using CSR format in the ComPASS code:   <br>
<center>
\f$
\left( \begin{array}{cccc}
10 & 8 & 2 & 0 & 0 \\
15 & 14 & 13 & 12 & 11 \\
60 & 64 & 0 & 0 & 0 \\
1 & 9 & 0 & 0 & 0 \\
30 & 35 & 25 & 0 & 0 \end{array} \right)
\f$
</center>

This matrix will hence be stored as:
<pre class="fragment">
    MyName\%Nb = 5
    MyName\%Pt(:)  = (0,3,8,10,12,15)
    MyName\%Num(:) = (10,8,2,15,14,13,12,11,60,64,1,9,30,35,25)
</pre><br>


# Own and Ghost elements stored with CSR format

ComPASS is a parallel code, thus there is a
\ref distributionpage "distribution" of the elements
nodes, faces, cells and wells. Each processor has its
own and ghost elements and it is necessary to identify
them. More precisely, it is useful to
have a list of the own and ghost objects,
but also to know the numero of the processor for
which the objects are own.                                 <br>

The local objects are stored in a matrix where the
number of lines is <b>the total number of processors</b>.
The matrix is filled as follow:
<ul>
<li> the first line contains the <b>own</b> objects; </li>
<li> the second line contains the <b>ghost objects which
are own for the processor <code>p1</code></b>;   </li>
<li> the third line contains <b>the ghost objects which
are own for the processor <code>p2>p1</code></b>; </li>
<li> etcetera </li>
</ul>

As the number of elements are not the same in each line,
the CSR format is convenient. The optional vector Val
is filled in this case, it contains the numero of the
concerned processor (<code>MycommRank, p1, p2, ...</code>).
Val is declared as <code>integer</code> and its size
is <code>Nb</code>.

<pre class="fragment">
    MyName\%Nb = Ncpus
    MyName\%Pt(1)  = 0
    MyName\%Pt(2) = Number of own objects
    MyName\%Val(1) = MycommRank    ! my CPU number
    do i=2,Ncpus
        MyName\%Pt(i+1) = MyName\%Pt(i) + Number of ghost objects which are own for processor p
        MyName\%Val(i) = p
    enddo
    MyName\%Num(:) = Numero of all objects
</pre>

The number of lines is the total number of processors
(and not only the neighbours of the MycommRank), then
some lines of the matrix can be empty (no ghost objects
are own for these processors). In this case
<code>Pt(i+1)=Pt(i)</code>.        <br>

\image html Type-CSR1.png "Illustration of the CSR storage with own and ghost elements, view from processor P0."
\image latex Type-CSR1.pdf "Illustration of the CSR storage with own and ghost elements." width=10cm


<br>
