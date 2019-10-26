Keywords used in the specification of objects 
=============================================


This is a list of keywords used to specify various parameters of Duo objects. 

\begin{itemize} 

\item{\verb!name!} is a text label which can be assigned to any object for 
reference in the output. The string must appear within quotation marks. 
\index{name} Examples: \begin{myverbatim} name "X 1Sigma+" name 
"<X1Sigma\|HSO\|A3Pi>" \end{myverbatim} 

\item{\verb!lambda!}  specifies the quantum number(s) $\Lambda$, i.e. 
projections of the electronic angular momentum onto the molecular axis, 
either for one (PECs) or two states (couplings). It must be an integral 
number and is allowed to be either positive or negative. %     the sign of 
$\Lambda$ is relevant when specifying couplings between degenerate states. 
\index{lambda} Examples: \begin{myverbatim} lambda 1 lambda 0 -1 
\end{myverbatim} The last example is relative to a coupling-type object and 
the two numbers refer to the bra and ket states. 

\item{\verb!sigma!} specifies the quantum number(s) $\Sigma$, i.e. the 
    projections of the total spin onto the molecular axis, either for one (diagonal) or two
    states (couplings). These values should be real ($-S\le \Sigma \le S$) and can be half-integral,
    where $S$ is the total spin. \verb!sigma! is currently required for the spin-orbit couplings only.
\index{sigma} Example: \begin{myverbatim} sigma 0.5 1.5 \end{myverbatim} 
where two numbers refer to the bra and ket states. 

\item{\verb!mult!} (alias: \verb!multiplicity!) specifies the multiplicity of 
    the electronic state(s), given by $(2S + 1)$, where $S$ is the total spin.
    It must be an integer number and is an alternative to the \verb!spin! keyword. \index{mult}
    \index{multiplicity}
Examples: \begin{myverbatim} mult 3 mult 1 3 \end{myverbatim} The last 
example is relative to a coupling-type object and the two numbers refer to 
the bra and ket states. 

\item{\verb!spin!} The total spin of the electronic state(s), an integer or 
half-integer number. \index{spin} Example: \begin{myverbatim} spin 1.0 spin 
0.5 1.5 \end{myverbatim} The last example is relative to a  coupling-type 
object and the two numbers refer to the bra and ket states. 

\item{\verb!symmetry!}  This keyword tells Duo if the electronic state has 
    gerade \verb!g! or ungerade \verb!u! symmetry (only for homonuclear diatomics)
    and whether it has positive (\verb!+!) or negative \verb!-! parity (only for
    $\Sigma$ states, i.e. states with $\Lambda=0$, for which it is mandatory).
     \index{Symmetries of objects}  %ma!g!u!+!-}
    \index{gerade} \index{ungerade} \index{+} \index{-}
Examples: \begin{myverbatim} symmetry + symmetry + u symmetry g 
\end{myverbatim} The \verb!g!/\verb!u! or \verb!+!/\verb!-! can appear in any 
order. 

\item{\verb!type!} selects the parametrised analytical function used for 
representing the objects 
                   or selects the interpolation type to be used. The function types supported by Duo
                   are listed in Section~\ref{s:functions}. \index{type}
Examples: %\item \verb!type ЕМО! \begin{myverbatim} type grid type polynomial 
type morse \end{myverbatim} In the examples above \verb!grid! selects 
numerical interpolation of values given on a grid, \verb!polynomial! selects 
a polynomial expansion and \verb!morse! selects a polynomial expansion in the 
Morse variable. See Section~\ref{s:functions} for details. 


\item{ \verb!Interpolationtype!} \index{Interpolationtype} is used only for 
\texttt{type grid} and specifies 
    the method used for the numerical interpolation of the numerical values.
    The currently implemented interpolation methods are
    \verb!Cubicsplines! and \verb!Quinticsplines! (default).
Example: \begin{myverbatim} Interpolationtype Cubicsplines Interpolationtype 
Quinticsplines \end{myverbatim} 

\item{\verb!factor!} \index{factor}  This optional keyword permits to rescale 
any object by 
    an arbitrary multiplication factor. At the moment the accepted values are any real number,
    the imaginary unit $i$, the square root of two, written as \verb!sqrt(2)!, or products
    of these quantities. To write a product simply leave a space between the factors, but do not
    use the \texttt{*} sign. All factor can have a $\pm$ sign.
    The default value for \verb!factor! is 1. This keyword is useful, for example,
    to temporarily zero a certain object without removing it from the input file.
Examples: \begin{myverbatim} factor 1.5 factor -sqrt(2) factor -2 sqrt(2) i 
\end{myverbatim} In the last example the factor is read in as $-2 \sqrt{2} 
i$. Note that imaginary factors make sense only in some cases for some 
coupling terms (in particular, spin-orbit) in the Cartesian-representation, 
see Section~\ref{s:representations}. 


\item{ \verb!units!} \index{units} This keyword selects the units of measure 
used for the %grid representation of 
    the object in question. Supported units are: \verb!angstroms!
    (default) and \verb!bohr! for the bond lengths; \verb!cm-1! (default),
    \verb!hartree! (aliases are \verb!au!, \verb!a.u.!, and \verb!Eh!), and \verb!eV! (electronvolts)
    for energies;
%     for \verb!poten!, \verb!spin-orbit!, \verb!spin-spin!, \verb!spin-rot!, 
\verb!diabatic!; 
   \verb!debye! (default) and \verb!ea0! (i.e., atomic units) for dipoles; units can appear in any order. %\verb!dipole!.
%     Note that the objects \verb!L+! are \verb!bobrot! unitless. %LL I think 
L+ has units hbar (same as angular momentum) Example: \begin{myverbatim} 
units angstrom cm-1 (default for poten, spin-orbit, lambda-doubling etc) 
units bohr cm-1 units debye  (default) units ae0 bohr \end{myverbatim} 

\item{\verb!ASSIGN_V_BY_COUNT!} The vibrational quantum number $v$ is 
assigned by counting the rovibronic states of the same `State', $\Lambda$, 
$\Sigma$ arranged by increasing energy. The corresponding `State', $\Lambda$, 
$\Sigma$ labels are defined using the largest-contribution approach (the 
quantum labels corresponding to the basis set contribution with the largest 
expansion coefficient).   The keyword should appear anywhere in the body of 
the input file. The default is to use the largest-contribution  approach also 
to assign the vibrational quantum number (no \verb!ASSIGN_V_BY_COUNT!). 



\item{\verb!values!} \index{values} This keyword starts the subsection 
containing the numerical 
    values defining the object. %specifying the object's grid (\verb!type grid!) or analytical (any
%     other \verb!type!s) representations. For one of the \verb!type!'s 
corresponding to an analytical function (Section~\ref{s:functions}), the 
input between \verb!values! and \verb!end! contains the values of the 
parameters of the function. The input consists in two columns separated by 
spaces containing \emph{(i)} a string label identifying the parameter and 
\emph{(ii)} the value of the parameter (a real number). 

In case of \verb!fitting! (see Section~\ref{sec:fitting}) a third column 
should also be provided; the parameters which are permitted to vary during 
fitting must have in the third column the string \texttt{fit} or, 
alternatively, the letter \texttt{f} or the number 1. Any other string or 
number (for example, the string \texttt{nofit} or the number 0) % containing 
% providing the switch $ = 0/1$ identifying if the parameters should be 
varied in the fit (1) % or not (0) (ignored when fitting is not activated, 
see Section~\ref{sec:fitting}). implies the parameter should be kept at its 
initial value. In the case of fitting, the keyword \verb!link! can be also 
appear at the end of each the line; this keyword permits to cross-reference 
values from different objects and is explained below in this section. 

In the case of objects of type \texttt{grid} only two columns are normally 
needed, %, the input between \verb!values! and \verb!end! should contain two 
a first containing the bond length $r_i$ and a second with the value of the 
object. Only in the case of object of the \verb!abinitio! (\verb!reference!) 
type and specified as \texttt{grid} a third column should be present 
specifying the fitting weights (see Section~\ref{sec:fitting}). 

\item{\verb!<x|Lz|y>!}, \verb!<z|Lz|xy>! (aliases \verb!<a|Lz|b>! and 
\verb!<1|Lz|2>!) \index{$\langle x \vert L_z \vert y \rangle$} This keyword 
     is sometimes needed when specifying coupling curves between electronic states
     with $|\Lambda| > 0$ in order to resolve ambiguities in the definition of the
     degenerate components of each electronic state, see Section~\ref{s:representations}.
%     switching to the alternative \red{MOLPRO?}, spherical harmonics 
\red{???} representation, which %     is also the value of the matrix element 
of the $\hat{L}_z$ operator computed for %     the two component spherical 
harmonic, degenerate functions $|x\rangle$ and %     $|y\rangle$ for the 
$\Pi$ states or $|z\rangle$ and $|xy\rangle$ for the $\Delta$ %     states 
etc. These values are used in order to transform the spherical harmonics %     
representation (as used e.g. by MOLPRO and Gaussian \red{XXX}) to the %     
$\Lambda$-representation as internally used by \duo. The corresponding 
\verb!<x|Lz|y>! values for both %     coupled states must be provided. %      
\red{This feature is under construction but should be finished by the paper 
submission.} 
      This keyword specifies the matrix element of the $\hat{L}_z$ operator between the degenerate components
      of the electronic wave function. %Specifically,
%       The spacial part of electronic states with $|\Lambda| >0 $ is doubly 
degenerate  and 
      Quantum chemistry programs such as Molpro choose the degenerate components so that they transform
      like the $x$ or $y$ functions (for states with odd $|\Lambda|$, i.e. $\Pi$, $\Phi$, $\cdots$, corresponding to
      symmetry species $b_1$ and $b_2$ in the C$_{2v}$ point group) or like $z$ and $xy$ (for states with even $|\Lambda|$,
      i.e. $\Delta$, $\Gamma$, $\cdots$, corresponding to  symmetry species $a_1$ and $a_2$ in the C$_{2v}$ point group).
      In this keyword we specify matrix elements of the type $\langle \Pi_x | \hat{L}_z| \Pi_y \rangle$ or
      $\langle \Delta_z | \hat{L}_z| \Delta_{xy} \rangle$  for the bra and ket states. % of the coupling curve in question.
Examples: \begin{myverbatim} <x|Lz|y>   i  -i <z|Lz|xy> -2i  i 
\end{myverbatim} These matrix elements are pure imaginary number in the form 
$\pm |\Lambda | i$. It is the overall $\pm$ sign which Duo needs and cannot 
be otherwise guessed. As shown in the examples above, each factor should be 
written in the form $\pm |\Lambda | i$ without any space or \texttt{*} sign. 


\item{ \verb!link!}  \index{link} This special keyword is used in fitting %     
to link (cross-reference) a set of parameters %(analytical representation 
only) 
    to force a set of parameters %(analytical representation only)
    (which may be relative to a different object) to have the same value.
    For example, in a typical situation one may want to fit a set of PECs and to constrain their
    dissociation (asymptotic) energy to the same value (because they are expected from theory to share the same
    dissociation channel).
%     if dissociation limits of states 1 and 2 are the same, the second state 
parameter \verb!De! %     can take the value of \verb!De! of the 1st state. 
It is useful during the fittings %     (see the addons section.) 
    After the keyword \texttt{link} one should provide three numbers $i_1$, $i_2$, $i_3$ defining the parameter ID, where
    $i_1$ identifies the object type (e.g. \texttt{poten}, \texttt{spin-orbit}, \texttt{spin-rot} etc.), $i_2$ is the object number within the type $i_1$ and $i_3$ is the parameter number as it appears after \texttt{values}. The ID numbers $i_1$,$i_2$, $i_3$ are specified in the fitting outputs in the form \texttt{[i,j,k]}. Example of the input:
\begin{myverbatim} DE            0.50960000000000E+05   fit     link   1   1   
3 \end{myverbatim} Example of the corresponding output \begin{myverbatim} DE            
0.50960000000000E+05   [ 1   1   3 ] \end{myverbatim} 


\item{\verb!morphing!} This keyword is used for fitting and switches on the 
morphing method, see Ref.~\cite{YuLoTe15}. \index{morphing} 
%Section~\ref{s:morphing} as well as 

\item{\verb!ZPE!}  allows to explicitly input the zero-point energy (ZPE) of 
the molecule (in \cm). This affects the value printed, as by default 
   Duo  prints energy of rovibronic levels by subtracting the ZPE. if not specified, the lowest energy of the first $J$-block (independent of parity) will be used as appear on the line \verb!Jlist!.

\item{\verb!fit_factor!} \index{fit-factor} This factor ($d_{\lambda}$) is 
used as a part of the reference \ai\ curves of the \texttt{abinitio} type 
which (when given) is applied to the corresponding weights assigned to the 
corresponding values of this object, (see Section 4.3 of \cite{YuLoTe15}). It 
is different from {\verb!fit_factor!} defined within the {\verb!Fitting!} 
section. 

Example: \begin{myverbatim} abinitio poten 1 name "A 1Pi" type   grid lambda 
1 mult   1 units bohr cm-1 fit_factor  1e1 values 2.00	32841.37010	0.01 2.20	
17837.88960	0.10 2.40	8785.33147	0.70 2.60	3648.35520	1.00 2.70	
2107.10737	1.00 2.80	1073.95670	1.00 2.90	442.52180	1.00 3.00	
114.94960	1.00 3.10	0.00000	    1.00 3.20	48.46120	1.00 3.30	
213.34240	1.00 3.40	455.16980	1.00 3.50	739.61170	1.00 3.60	
1038.82620	1.00 3.70	1332.46170	1.00 4.00	2059.31119	1.00 4.50	
2619.19233	0.30 5.00	2682.84741	0.30 6.00	2554.34992	0.30 8.00	
2524.31106	0.30 10.00	2561.48269	1.00 12.00	2575.09861	1.00 end 
\end{myverbatim} 

\end{itemize}
