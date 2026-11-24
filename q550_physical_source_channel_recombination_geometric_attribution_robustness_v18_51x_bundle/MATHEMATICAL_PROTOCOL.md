# Q550 — v18.51x | Mathematical Protocol

**Question.** Does the near-right angle of the v18.51w Green-channel responses survive physically transparent reassignments of the source-channel split, holding the *total* original κ-derivative of the 10D vector field fixed?

## Authenticated baseline

The v18.51t dynamics defined the 10D parameter-forcing channels

- \(b_{v_2}(\xi)=e_5 b_5(\xi)\),
- \(b_{\omega_2}(\xi)=e_7 b_7(\xi)\),
- \(b_{\rm other}=b-b_{v_2}-b_{\omega_2}\).

Indices are zero-based and refer to the historical 10D state. Their Green responses are the columns \(S_{v_2}\), \(S_{\omega_2}\), and \(S_{\rm other}\) of the independently integrated v18.51t transport, later resolved into radial and azimuthal components by v18.51w. The physically defined split is the **baseline**, not an arbitrary original choice.

## Sum-preserving channel redefinition

At every \(\xi\), apply the same source-label transformation to the *source vector fields and all their linear responses*:

\[
[\widetilde b_1\;\widetilde b_2]=[b_{v_2}\;b_{\omega_2}]M_\eta,
\qquad M_\eta=\begin{pmatrix}1-\eta&\eta\\\eta&1-\eta\end{pmatrix}.
\]

\(M_\eta(1,1)^T=(1,1)^T\) and \(\det M_\eta=1-2\eta\). Thus \(\widetilde b_1+\widetilde b_2=b_{v_2}+b_{\omega_2}\), and linearity gives

\[
[\widetilde S_1\;\widetilde S_2]=[S_{v_2}\;S_{\omega_2}]M_\eta,
\quad \widetilde S_1+\widetilde S_2+S_{\rm other}=S_{\rm total}.
\]

This is tested in the **original 10D source, independent forced 10D state response, Green 10D response, Cartesian 5D response, local radial--azimuthal response and 2D lab frame**. The sums should agree to the numerical accuracy of the stored results. Both historical seed families are included.

For second-plane vectors \(a,b\),

\[
\det(\widetilde a,\widetilde b)=\det(M_\eta)\det(a,b),
\quad \angle(\widetilde a,\widetilde b)=\operatorname{atan2}\big(\det(\widetilde a,\widetilde b),\widetilde a\cdot\widetilde b\big).
\]

The signed-area *number* is not channel-recombination-invariant; its sign is preserved only while the channel transformation preserves channel-basis orientation. The split becomes singular at \(\eta=1/2\), excluded from invertible-basis claims.

The projection-cancellation percentage is evaluated on the signed natural radial components of the three source-labelled response vectors:

\[
C_r=1-\frac{|r_1+r_2+r_{\rm other}|}{|r_1|+|r_2|+|r_{\rm other}|}.
\]

It is undefined at zero response; the origin is masked.

The critical mixing at a given nondegenerate section for vanishing radial projection of the *second redefined* channel is

\[
\eta_*(\xi)=\frac{-b_r(\xi)}{a_r(\xi)-b_r(\xi)}.
\]

## Scientifically admissible claims

The sum-preserving transformation represents a change of **source attribution/bookkeeping**. The two redefined channels no longer occupy exclusively the original physical acceleration-component supports. Therefore they are **not** experimentally independent actuators or proof that the historical physical split is incorrect. The v18.51t channel definition retains its concrete model-coordinate meaning. This audit tests whether geometry claimed as a *property of the total response* survives its source decomposition, not whether the original source physics is itself arbitrary.

Agreement of the recomposition formulas is primarily linear algebra, not an independent nonlinear replay or a constitutive-memory certificate. T/W cross-archive equality and matched source-vs-response recomposition ground that algebra in the authentic model data. To establish a definition-independent causal mechanism would require a physically justified admissible class of decompositions or independently controlled perturbations.
