Jacobian and second member    {#resjac}
===========================





## Notations

Concerning the mesh:
 - \f${\cal M}\f$ : set of all cells of the mesh

 - \f${\cal F} \f$ : set of all faces
of the mesh, <b>remark
that the faces are not assumed to be planar</b>

 - \f${\cal F}_{K}\f$ : set of faces
of the cell \f$K\in{\cal M}\f$, thus
\f${\cal F} = \bigcup\limits_{K\in {\cal M}}{\cal F}_K\f$

 - \f${\cal F}_\Gamma \f$ : set of
all fracture faces of the mesh

 - \f${\cal E} \f$ : set of all edges of the mesh

 - \f${\cal V}\f$ : set of all vertices of the mesh

 - \f${\cal M}_s\f$ : set of cells containing
the node \f$s \in {\cal V}\f$,
respectively \f${\cal M}_f\f$ is the set of
cells containing the fracture face
\f$f \in {\cal F}_\Gamma \f$

 - \f${\cal F}_{\Gamma,K} = {\cal F}_K \cap {\cal F}_\Gamma\f$
: set of fracture faces of the cell \f$K\in{\cal M}\f$,
respectively \f${\cal F}_{\Gamma,s}\f$ is the set
of fracture faces containing the node
\f$s \in {\cal V}\f$

 - \f${\cal V}_{K}\f$ : set of vertices of
the cell \f$K\in{\cal M}\f$, respectively
\f${\cal V}_{f}\f$ is the set of vertices of
the fracture face \f$f \in {\cal F}_\Gamma\f$

 - \f${\bf x}_K\f$ : "center" of the
cell \f$K\in{\cal M}\f$ under the
assumption that \f$K\f$ is star-shaped
with respect to \f${\bf x}_K\f$



Concerning the components and the phases:
 - \f$ {\cal P} \f$ : set of possible phases
 - \f$ N_c\f$ : number of components
 - \f$ {\cal C} \f$ : set of possible components
 - \f$ {\cal Q} \f$ : ???? the unknown repesenting the set of
present phases
 - \f$ P_i = \{ \alpha \in {\cal P} | {\cal Q}_\alpha = 1 \} \subset {\cal P} \f$
: set of present phases of context \f$ {\cal Q} \f$.

Concerning the physical parameters:












# Discrete system

The discrete system couples the following set of equations:

### Discrete balance equations on the cells

The discrete balance equations on each cell \f$ K \in {\cal M} \f$
and for each component \f$ i \in {\cal C} \f$:
<center>\f$
\frac{n_i(X_K^n) - n_i(X_K^{n-1})}{\Delta t^n}        +   \sum\limits_{s \in {\cal V}_K} \sum\limits_{\alpha \in {\cal Q}_{K_{K,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,K,{\bf s}}^\alpha  +     \sum\limits_{f \in {\cal F}_{\Gamma,K}} \sum\limits_{\alpha \in {\cal Q}_{K_{K,f}^{\alpha,n}}^n \cap P_i} G_{i,K,f}^\alpha       = 0
\f$   </center>



### Discrete balance equations on the vertices

The discrete balance equations on each vertex
\f$ {\bf s} \in {\cal V}_{\it int} \cup {\cal V}_N \f$
and for each component \f$ i \in {\cal C} \f$:
<center>\f$
\frac{n_i(X_{\bf s}^n) - n_i(X_{\bf s}^{n-1})}{\Delta t^n}   +   \sum\limits_{K \in {\cal M}_{\bf s}} \sum\limits_{\alpha \in {\cal Q}_{K_{k,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,{\bf s},K}^\alpha +   \sum\limits_{f \in {\cal F}_{\Gamma,s}} \sum\limits_{\alpha \in {\cal Q}_{K_{f,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,{\bf s},f}^\alpha   = 0
\f$   </center>



### Discrete balance equations on the fracture faces
The discrete balance equations on each fracture face
\f$ f \in {\cal F}_\Gamma \f$
and for each component \f$ i \in {\cal C} \f$:
<center>\f$
\frac{n_i(X_f^n) - n_i(X_f^{n-1})}{\Delta t^n}      + \sum\limits_{s \in {\cal V}_f} \sum\limits_{\alpha \in {\cal Q}_{K_{f,s}^{\alpha,n}}^n \cap P_i} G_{i,f,s}^\alpha       + \sum\limits_{K \in {\cal M}_f} \sum\limits_{\alpha \in {\cal Q}_{K_{f,K}^{\alpha,n}}^n \cap P_i} G_{i,f,K}^\alpha   = 0
\f$   </center>

### Closure Laws
The local closure laws for all control volume
\f$
K \in {\cal M} \cup {\cal V}_{\it int} \cup {\cal V}_N \cup {\cal F}_\Gamma
\f$


### Flash
The flash for all control volume
\f$
K \in {\cal M} \cup {\cal V}_{\it int} \cup {\cal V}_N \cup {\cal F}_\Gamma
\f$






# Algorithm for the computation of the fluxes

In this model, the control volume are the Nodes,
the Cells and the Fracture Faces. Thus the fluxes
will be computed over all this objects. In the code,
due to the previous balance equations on each control
volume, we should compute the six double sums.
But for all control volume Dof1, Dof2, we have:
<center>\f$
G_{i,Dof1,Dof2}^\alpha = - G_{i,Dof2,Dof1}^\alpha
\f$</center>

Then we can divide by two the number of
sum that are necessary to compute.
It is not necessary to have one loop
per balance equations, two are enough:
with a loop over the cells and a loop over
the fracture faces, all the fluxes are computed:

<ul>
<li> for \f$ K \in {\cal M}    \f$   <br>
    <ul>
    <li> for \f$ s \in {\cal V}_K   \f$,
    compute
    \f$ \sum\limits_{\alpha \in {\cal Q}_{K_{K,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,K,{\bf s}}^\alpha \f$,
    then deduce
    \f$ \sum\limits_{\alpha \in {\cal Q}_{K_{K,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,{\bf s},K}^\alpha \f$.  </li>
     <li> for \f$ f \in {\cal F}_{\Gamma,K} \f$,
     compute
     \f$ \sum\limits_{\alpha \in {\cal Q}_{K_{K,f}^{\alpha,n}}^n \cap P_i} G_{i,K,f}^\alpha \f$,
     then deduce
     \f$ \sum\limits_{\alpha \in {\cal Q}_{K_{K,f}^{\alpha,n}}^n \cap P_i} G_{i,f,K}^\alpha \f$.  </li>
</li> </ul>
<li> for \f$ f \in {\cal F}_{\Gamma} \f$  <br>
    <ul>
    <li> for \f$ s \in {\cal V}_f   \f$,
    compute
    \f$ \sum\limits_{\alpha \in {\cal Q}_{K_{f,s}^{\alpha,n}}^n \cap P_i} G_{i,f,s}^\alpha \f$,
     then deduce
    \f$ \sum\limits_{\alpha \in {\cal Q}_{K_{f,s}^{\alpha,n}}^n \cap P_i} G_{i,s,f}^\alpha \f$.  </li>
</li></ul>
</ul>




## Implementation of the Residue

In the code the residue is decomposed in two parts:
 - the residue from terms of accumulation :
\f$   \frac{n_i(X^n) - n_i(X^{n-1})}{\Delta t^n}  \f$

 - the residu from conservation components.



Here we will focuse on the computation of the second
part of the residu, and because it appears in every
balance equation we will focuse more precisely on the sum :
<center>\f$
\sum\limits_{\alpha \in {\cal Q}_{K_{Dof1,Dof2}^{\alpha,n}}^n \cap P_i} G_{i,Dof1,Dof2}^\alpha
\f$   </center>


In the implementation the sum is cut in two parts
due to the following : <br>
<center>\f$
\{ \alpha \in {\cal Q}_{K_{Dof1,Dof2}^{\alpha,n}}^n \} = \{ (\alpha | \alpha \in {\cal Q}_{Dof1}^n ) \cap (\alpha | K_{Dof1,Dof2}^{\alpha,n}=K) \} \cup \{ (\alpha | \alpha \in {\cal Q}_{Dof2}^n) \cap (\alpha | K_{Dof1,Dof2}^{\alpha,n}=s) \}
\f$   </center>


The following algorithm fill the vectors
ResiduCell, ResiduNode and ResiduFrac such that:
<center>\f$

\begin{array}{ccc}

ResiduCell &=&  \sum\limits_{s \in {\cal V}_K} \sum\limits_{\alpha \in {\cal Q}_{K_{K,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,K,{\bf s}}^\alpha  &+&     \sum\limits_{f \in {\cal F}_{\Gamma,K}} \sum\limits_{\alpha \in {\cal Q}_{K_{K,f}^{\alpha,n}}^n \cap P_i} G_{i,K,f}^\alpha     \\
ResiduNode &=& \sum\limits_{K \in {\cal M}_{\bf s}} \sum\limits_{\alpha \in {\cal Q}_{K_{k,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,{\bf s},K}^\alpha &+&   \sum\limits_{f \in {\cal F}_{\Gamma,s}} \sum\limits_{\alpha \in {\cal Q}_{K_{f,{\bf s}}^{\alpha,n}}^n \cap P_i} G_{i,{\bf s},f}^\alpha        \\
ResiduFrac &=& \sum\limits_{s \in {\cal V}_f} \sum\limits_{\alpha \in {\cal Q}_{K_{f,s}^{\alpha,n}}^n \cap P_i} G_{i,f,s}^\alpha       &+& \sum\limits_{K \in {\cal M}_f} \sum\limits_{\alpha \in {\cal Q}_{K_{f,K}^{\alpha,n}}^n \cap P_i} G_{i,f,K}^\alpha

\end{array}
\f$   </center>


<pre class="fragment">
do k=1, NbCellLocal_Ncpus(commRank+1)                                       ! cell k
    do s=1, NbNodeCell(k)                                                   ! s is node in cell k
        nums = num_Node(s)
        Flux_ks(:) = 0.d0
!
        do m=1, NbPhasePresente_ctx( IncCell(k)\%ic)                         ! Q_k
            if( FluxDarcyKI(mph,s,k)>=0.d0) then                            ! K_{k,s}^{alpha}=k
                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then                                 ! cap P_i
                        Flux_ks(icp) = Flux_ks(icp) + G_{i,K,s}^alpha
                   endif
                enddo
            endif
        enddo
!
        do m=1, NbPhasePresente_ctx(IncNode(nums)\%ic)                       ! Q_s
            if(FluxDarcyKI(mph,s,k)<0.d0) then                              ! K_{k,s}^{alpha}=s
                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then                                 ! cap P_i
                        Flux_ks(icp) = Flux_ks(icp) + G_{i,K,s}^alpha
                   endif
                enddo
            endif
        enddo
!
        ResiduCell(1:NbComp,k)    = ResiduCell(1:NbComp,k)    + Flux_ks(:)  ! G_{i,K,s}^alpha
        ResiduNode(1:NbComp,nums) = ResiduNode(1:NbComp,nums) - Flux_ks(:)  ! - G_{i,K,s}^alpha
!
    end do ! end of node s in cell k
!
    do f=1, NbFracCell(k)                                                   ! f is frac in cell k
        numf = num_Frac(f)
        Flux_kf(:) = 0.d0
!
        do m=1, NbPhasePresente_ctx(IncCell(k)\%ic)                          ! Q_k
            if( FluxDarcyKI(mph,f,k)>=0.d0) then                            ! K_{k,f}^{alpha}=k
                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then                                 ! cap P_i
                        Flux_kf(icp) = Flux_kf(icp) + G_{i,K,f}^alpha
                   endif
                enddo
            endif
        enddo
!
        do m=1, NbPhasePresente_ctx(IncNode(numf)\%ic)                       ! Q_f
            if(FluxDarcyKI(mph,f,k)<0.d0) then                              ! K_{k,f}^{alpha}=f
                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then                                 ! cap P_i
                        Flux_kf(icp) = Flux_kf(icp) + G_{i,K,f}^alpha
                   endif
                enddo
            endif
        enddo
!
        ResiduCell(1:NbComp,k)    = ResiduCell(1:NbComp,k)    + Flux_kf(:)  ! G_{i,K,f}^alpha
        ResiduFrac(1:NbComp,numf) = ResiduFrac(1:NbComp,numf) - Flux_kf(:)  ! - G_{i,K,f}^alpha
!
    end do ! end of frac face f in cell k
end do ! cell k
!
!
!
do f=1, NbFracLocal_Ncpus(commRank+1)                                       ! fracture face f
    do s=1, NbNodeFrac(f)                                                   ! s is node in fracture face f
        nums = num_Node(s)
        Flux_fs(:) = 0.d0
!
        do m=1, NbPhasePresente_ctx( IncFrac(f)\%ic)                         ! Q_f
            if( FluxDarcyKI(mph,s,f)>=0.d0) then                            ! K_{f,s}^{alpha}=f
                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then                                 ! cap P_i
                        Flux_fs(icp) = Flux_fs(icp) + G_{i,f,s}^alpha
                   endif
                enddo
            endif
        enddo
!
        do m=1, NbPhasePresente_ctx(IncNode(nums)\%ic)                       ! Q_s
            if(FluxDarcyKI(mph,s,f)<0.d0) then                              ! K_{f,s}^{alpha}=s
                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then                                 ! cap P_i
                        Flux_fs(icp) = Flux_fs(icp) + G_{i,f,s}^alpha
                   endif
                enddo
            endif
        enddo
!
        ResiduFrac(1:NbComp,f)    = ResiduFrac(1:NbComp,f)    + Flux_fs(:)  ! G_{i,f,s}^alpha
        ResiduNode(1:NbComp,nums) = ResiduNode(1:NbComp,nums) - Flux_fs(:)  ! - G_{i,f,s}^alpha
!
    end do ! end of node s in cell k
enddo  ! frac face f
</pre>





# Jacobian

First, each processor fills a matrix with
all the control volume:

<center>\f$
\begin{array}{ccccccc}
          & node & frac & cell & well\ inj & well\ prod &   &               \\
          & a11  & a12  & a13  & a14     & a15      &   & node\ own      \\
JacBigA = & a21  & a22  & a23  & 0       & 0        &   & frac\ own      \\
          & a31  & a32  & a33  & 0       & 0        &   & cell\ (own\ and\ ghost)   \\
          & a41  & 0    & 0    & a44     & 0        &   & well\ inj\ own   \\
          & a51  & 0    & 0    & 0       & a55      &   & well\ prod\ own
\end{array}
\f$</center>


Then there is a regularization, then a Schur.
<center>\f$
\begin{array}{cccccc}
          & node & frac & well\ inj & well\ prod &   &               \\
          & A11  & A12  & A13     & A14      &   & node\ own      \\
JacA =    & A21  & A22  & 0       & 0        &   & frac\ own      \\
          & A31  & 0    & A33     & 0        &   & well\ inj\ own   \\
          & A41  & 0    & 0       & a44      &   & well\ prod\ own
\end{array}
\f$</center>





# Discretisation of multiphase compositional Darcy flows

VAG scheme.

notations ???

One main point of the VAG scheme is to use generalised fluxes
\f$
F_{K,s}(u) = - F_{s,K}(u)
\f$
between a cell

The phase pressures at each control volume
\f$ K \in \kappa \f$
are defined by
<center>
\f$
P_K^\alpha = P_K + P_{c,\alpha}(S_K)
\f$
</center>
for all phases \f$ \alpha \in {\cal P} \f$
(and not only for present phases  \f$ \alpha \in {\cal Q}_K \f$).  <br>

Then, for all \f$ s \in {\cal V}_K \f$, \f$ K \in {\cal M} \f$,
and for all phases \f$ \alpha \in {\cal Q}_K \cup {\cal Q}_s \f$,
.... we define the Darcy fluxes
<center>
\f$
V_{K,s}^\alpha(X_K,X_{{\cal V}_K}) = - V_{s,K}^\alpha(X_K,X_{{\cal V}_K})
 = ...
\f$   </center>









## Computation of the fluxes


The expression of the matrix (resp. the fracture) fluxes
is linear and local to the cell (resp. fracture face).
More precisely, the matrix fluxes are given by

<center>\f$
F_{K,\nu}(u_{\cal D}) = \sum\limits_{\nu'\in \Xi_K} T_K^{\nu,\nu'} (u_K - u_{\nu'}),
\f$ </center>
with a symmetric positive definite
transmissibility matrix
\f$ T_K = (T_K^{\nu,\nu'})_{(\nu,\nu')\in \Xi_K\times \Xi_K} \f$
depending only on the cell \f$K\f$ geometry and
on the cell matrix diffusion tensor.
The fracture fluxes are given by  <br>
<center>\f$
F_{\sigma,\bf s}(u_{\cal D}) = \sum\limits_{s\in {\cal V}_\sigma} T_\sigma^{\bf s,\bf s'} (u_\sigma - u_{\bf s'}),
\f$</center>
with a symmetric positive definite
transmissibility matrix
\f$ T_\sigma = (T_\sigma^{\bf s,\bf s'})_{(\bf s,\bf s')\in {\cal V}_\sigma\times {\cal V}_\sigma} \f$
depending only on the fracture face
\f$ \sigma \f$ geometry and on the
fracture face tangential diffusion tensor.





The Darcy fluxes taken into account the gravity term are defined by
<center>\f$
\begin{equation}
\left\{\begin{array}{r@{\,\,}c@{\,\,}ll}
&V^\alpha_{K,\nu}(X_{\cal D}) &=& F_{K,\nu}(P_{\cal D}) + \rho_{K,\nu}^\alpha g F_{K,\nu}(Z_{\cal D}), \\
&V^\alpha_{\sigma,\bf s}(X_{\cal D}) &=& F_{\sigma,\bf s}(P_{\cal D}) + \rho_{\sigma,\bf s}^\alpha g F_{\sigma,\bf s}(Z_{\cal D}),
\end{array}\right.
\end{equation}
\f$ </center>

where \f$Z_{\cal D}\f$ denotes the vector \f$(z_\nu)_{\nu \in {\cal M} \cup {\cal V} \cup {\cal F}_\Gamma}\f$,
and the phase mass density is defined by the average
<center>\f$
\rho_{\mu,\nu}^\alpha =
{\rho^\alpha(P_\mu,T_\mu,C^\alpha_\mu) + \rho^\alpha(P_\nu,T_\nu,C^\alpha_\nu) \over 2 }
\f$.</center>


It could also be possible to compute \f$\rho\f$ as follow:
<center>\f$
\rho_{\mu,\nu}^\alpha =
...\f$.</center>


For all \f$ s \in {\cal V}_K \f$, \f$ K \in {\cal M} \f$,
and \f$ \alpha \in {\cal Q}_K \cup {\cal Q}_s \f$,
we define the fluxes
<center>\f$
\begin{array}{r@{\,\,}c@{\,\,}ll}
& G_{i,K,s}^\alpha &=& - G_{i,s,K}^\alpha    \\
& &=& ( C_i^\alpha \frac{\zeta^\alpha k_{r_\alpha}}{\mu^\alpha}  )(X_{K_{k,s}^{\alpha,n}}^n) V_{K,{\bf s}}^\alpha(X_K^n,X_{{\cal V}_K}^n)
\end{array}
\f$   </center>
