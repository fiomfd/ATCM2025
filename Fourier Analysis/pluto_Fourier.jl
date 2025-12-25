### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 581f8ff0-d1fc-11ec-3cab-5f335ec587ad
begin
	using PlutoUI 
	using Plots
	using Colors
	using LaTeXStrings

	# Small patch to make images look more crisp:
	# https://github.com/JuliaImages/ImageShow.jl/pull/50
	Base.showable(::MIME"text/html", ::AbstractMatrix{<:Colorant}) = false

	nbsp = html"&nbsp"
end

# ╔═╡ 25063436-79e0-4f01-ab92-f01defb4e948
md"""
#### Fourier Analysis
##### [Hiroyuki Chihara](https://fiomfd.github.io/) (University of the Ryukyus, Okinawa Island, Japan)

##### Fourier series
1. Definition of Fourier series
2. Real Form of the Fourier series
3. Example of Fourier series
4. Dirichelet kernel
5. Fourier series of periodic continuous functions
6. Fourier series of triangular function
7. A continuous function whose fourier series does not converges even pointwisely
8. Step function and sawtooth function
9. Decomposition of piecewise smooth functions
10. Weyl's equidistribution theorem
11. Continuous but nowhere differentiable functions

##### Fourier Series and Fourier transform
12. Shannon's sampling theorem and aliasing

"""

# ╔═╡ 79521e83-f32a-48aa-bdf1-27f3f95ecb9e
md"""
## Fourier series
"""

# ╔═╡ c299148e-4b85-4bf3-a514-20d1dff0dc59
md"""
#### 1. Definition of Fourier series
Throughout this note we consider piecewise continuous $1$-periodic functions $f(x)$ on $\mathbb{R}$, that is, 

-  $f(x+1)=f(x)$ for all $x\in\mathbb{R}$.

- There exist numbers $0=a_0<a_1<\dotsb<a_m=1$ such that $f(x)$ is continuous on each open interval $(a_{j-1},a_j$, $j=1,\dotsc,m$, and limits exist:
$f(a_0+0), \quad f(a_1-0), \quad f(a_1+0), \quad \dotsc \quad f(a_m-0).$

We remark that if we restrict $f(x)$ on each closed interval $[a_{j-1}+k,a_{j}+k]$ for $j=1,\dotsc,m$ and $k\in\mathbb{Z}$, $f(x)$ can be regarded as a continuous function. We denote by $\mathcal{L}$ the set of all piecewise continuous $1$-periodic functions on $\mathbb{R}$. 

All the $e^{2\pi i nx}$, $n\in\mathbb{Z}$ are $1$-periodic smooth functions. In the spirit of Joseph Fourier (1768-1830), we consider whether $f(x)\in\mathcal{L}$ can be expressed by a trigonometric series of the form

$\sum_{n=-\infty}^{\infty}c_ne^{2\pi i nx}.$

If $f(x)=\displaystyle\sum_{n=-\infty}^{\infty}c_ne^{2\pi i nx}$ holds, then we have by formal computation 

$\begin{aligned}
  \int_0^1
  f(x)
  e^{-2\pi i mx}
  dx
& =
  \int_0^1
  \sum_{n=-\infty}^{\infty}
  c_n
  e^{2\pi i (n-m)x}
  dx
\\
& =
  \sum_{n=-\infty}^{\infty}
  c_n
  \int_0^1
  e^{2\pi i (n-m)x} 
  dx
\\
& =
  \sum_{n=-\infty}^{\infty}
  c_n
  \delta_{mn}
  =
  c_m.
\end{aligned}$

Note that the condition $f(x)\in\mathcal{L}$ guarantees the existence of the left hand side of the above. The formal trigonometric series defined by 

$\sum_{n=-\infty}^{\infty}c_ne^{2\pi i nx},
\quad
c_n=\int_0^1f(x)e^{-2\pi i nx}dx$

is said to be the Fourier series of $f(x)$, and we denote this by 

$f(x)\sim\sum_{n=-\infty}^{\infty}c_ne^{2\pi i nx}.$

 $c_n$ is called the $n$-th Fourier coefficient of $f(x)$.

"""

# ╔═╡ 6f6479b6-265a-4474-8a90-3df0ff98cfac
md"""
#### 2. Real Form of the Fourier series
If $f(x)$ is real valued, then the Fourier series can be written by the real-valued trigonometric series. In this case $\overline{c_n}=c_{-n}$ holds, and $c_0$ is a real number. If we set for $n=0,1,2,\dotsc$ 

$\overline{c_n}=c_{-n}=\frac{a_n+ib_n}{2},$

$a_n
=
2
\int_0^1
f(x)
\cos(2\pi nx)
dx,$

$b_n
=
2
\int_0^1
f(x)
\sin(2\pi nx)
dx.$

then $b_0=0$, and 

$f(x)
\sim
\frac{a_0}{2}
+
\sum_{n=1}^\infty
\bigl(
a_n\cos(2\pi nx)
+
b_n\sin(2\pi n x)
\bigr).$

It is convenient to compute 

$\frac{a_0}{2}
=
c_0
=
\int_0^1
f(x)
dx,$

$a_n+ib_n
=
2c_{-n}
=
2
\int_0^1
f(x)
e^{2\pi i nx}
dx,
\quad
n=1,2,3,\dotsc.$


For $\lambda=\mu+\nu\in\mathbb{C}$, $\mu,\nu\in\mathbb{R}$, we have 

$\begin{aligned}
  e^{\lambda x}
& =
  e^{\mu x + i\nu x}
  =
  e^{\mu x}
  \{\cos(\nu x) + i\sin(\nu x)\},
\\
  \frac{d}{dx}
  e^{\lambda x}
& =
  \mu
  e^{\mu x}
  \{\cos(\nu x) + i\sin(\nu x)\}
  +
  \nu
  e^{\mu x}
  \{-\sin(\nu x) + i\cos(\nu x)\}
\\
& =
  \mu
  e^{\mu x}
  \{\cos(\nu x) + i\sin(\nu x)\}
  +
  i\nu
  e^{\mu x}
  \{i\sin(\nu x) + \cos(\nu x)\}
\\
& =
  (\mu+i\nu)
  e^{\mu x}
  \{\cos(\nu x) + i\sin(\nu x)\}
  =
  \lambda
  e^{\lambda x},
\end{aligned}$

which is useful for computing Fourier coefficients. 

"""

# ╔═╡ fd9270b6-bcb1-444a-b2e8-bc07efee3a24
md"""
#### 3. Example of Fourier series
Let $f(x)$ be a function defined by 

$f(x):=\min\{x-[x],1-x-[1-x]\}, \quad x\in\mathbb{R},$

$[s]:=\max\{m\in\mathbb{Z} : m \leqq s\}, \quad x\in\mathbb{R}.$

In other words, $f(x)$ is a $1$-periodic Lipshitz continuous function such that 

$f(x)
=
\begin{cases}
x, &\ x\in[0,1/2),
\\
1-x, &\ x\in[1/2,1). 
\end{cases}$

We shall confirm that the Fourier series of $f(x)$ is

$f(x)
\sim
\frac{1}{4}
-
\sum_{k=1}^\infty
\frac{2}{(2k-1)^2\pi^2}
\cos\bigl(2(2k-1)\pi x\bigr).$

Indeed

$\frac{a_0}{2}
=
\int_0^{1/2}xdx+\int_{1/2}^1(1-x)dx
=
\frac{2}{2^3}
=
\frac{1}{4},$

and for $n=1,2,3,\dotsc$, 

$\begin{aligned}
  a_n+ib_n
& =
  2
  \int_0^{1/2}xe^{2\pi i nx}dx
  +
  2
  \int_{1/2}^1(1-x)e^{2\pi i nx}dx
\\
& =
  \left[
  \frac{xe^{2\pi i nx}}{n\pi i}
  \right]_0^{1/2}
  +
  \left[
  \frac{(1-x)e^{2\pi i nx}}{n\pi i}
  \right]_{1/2}^1
  -
  \int_0^{1/2}
  \frac{e^{2\pi i nx}}{n\pi i}
  dx
  +
  \int_{1/2}^1
  \frac{e^{2\pi i nx}}{n\pi i}
  dx
\\
& =
  0
  -
  \left[
  \frac{e^{2\pi i nx}}{2(n\pi i)^2}
  \right]_0^{1/2}
  +
  \left[
  \frac{e^{2\pi i nx}}{2(n\pi i)^2}
  \right]_{1/2}^1
\\
& =
  \frac{e^{\pi i n}-1-1+e^{\pi i n}}{2n^2\pi^2}
  =
  \frac{\cos(\pi n)-1}{n^2\pi^2}
  =
  -
  \frac{1-(-1)^n}{n^2\pi^2}
\\
& =
  \begin{cases}
  -
  \dfrac{2}{(2k-1)^2\pi^2},
  &\ n=2k-1,
  \\
  0
  &\ n=2k,
  \end{cases}
  \quad
  k=1,2,3,\dotsc.
\end{aligned}$



"""

# ╔═╡ a94cb318-c3e9-4285-af58-154d6955f2d5
md"""
#### 4. Dirichelet kernel
Suppose $f(x)\in\mathcal{L}$ and $f(x)\sim\displaystyle\sum_{n=-\infty}^\infty c_ne^{2\pi i nx}$. Our basic problems are the following.

- Q1. Does the formal series $\displaystyle\sum_{n=-\infty}^\infty c_ne^{2\pi i nx}$ converge?

- Q2. If the formal series $\displaystyle\sum_{n=-\infty}^\infty c_ne^{2\pi i nx}$ converges, $f(x)=\displaystyle\sum_{n=-\infty}^\infty c_ne^{2\pi i nx}$ holds?

We introduce the partial $N$-the sum of the series

$S_N(s):=\sum_{n=-N}^{n=N}c_ne^{2\pi i nx}.$

We have 

$S_N(x)
=
\int_0^1D_N(x-y)f(y)dy
=
\int_{-1/2}^{1/2}D_N(t)f(x+t)dt,$

$D_N(t)
=
\sum_{n=-N}^{n=N}
e^{2\pi i nt}
=
1+
2\sum_{n=1}^N
\cos(2\pi n x)$

by the definition of the Fourier coefficients. $D_N(t)$ is called the $N$-th Dirichelet kernel. $D_N(t)$ is $1$-periodic smooth function and 

$\int_0^1D_N(t)dt=1.$

Moreover we obtain 

$D_N(t)
=
\begin{cases}
\dfrac{\sin\bigl(\pi(2N+1)t\bigr)}{\sin(\pi t)},
&\ t\not\in\mathbb{Z},
\\
2N+1,
&\ t\in\mathbb{Z}.
\end{cases}$
 
This expression is useful to study the convergence and divergence of Fourier series. Here we observe the behavior of the Dirichelet kernel in a one-period interval $[-1/2,1/2]$. Roughly speaking 

- The mass is concentrating on an interval $I_N:=[-1/(2N+1),1/(2N+1)]$ as $N\rightarrow\infty$.

-  $D_N(t)$ oscillates more violently in $[-1/2,1/2]\setminus{I_N}=[-1/2,-1/(2N+1))\cup(1/(2N+1),1/2]$ as $N\rightarrow\infty$.

"""

# ╔═╡ 3edf319a-48e4-48f0-82ab-f81e5f2c581c
begin
	M0=100; 
	t = range(-0.6, 0.6, length = 301);
	S0 = ones(M0+1,301);
	for m=2:M0+1
		for k=1:301
		    S0[m,k]=S0[m-1,k]+2*cos(2*pi*(m-1)*t[k]);
		end
	end
end

# ╔═╡ a1aa2ce7-db42-44d2-a9d8-4ae67976f76f
md"""
 $N$ = $(@bind N0 Slider(0:M0, show_value=true))
"""

# ╔═╡ cceffb37-cb1f-4436-aab6-17f325504276
begin
    plot(t,S0[N0+1,:],
		 grid=false,
		 linewidth=2,
		 ylim=(-M0/2,2*M0+M0/10),
         title="Dirichelet Kernel \$D_N(t)\$",
         xticks = ([-0.5 0 0.5;], [-0.5,0,0.5]),
		 yticks = ([0 100 200;], [0 100 200]),
		 xlabel="t",
         label=false,
	     legend=:false)
end

# ╔═╡ 897a6954-9e77-41d8-a88d-d43b4c9902c6
md"""
#### 5. Fourier series of periodic continuous functions 
We have the following results on the convergence and divergence of the Fourier series of continuous function. 

1. Suppose that $f(x)$ is $1$-periodic Hoelder continuous of degree $\alpha\in(0,1]$, that is, there exists a constant $L>0$ such that $\lvert{f(x)-f(y)}\rvert \leqq L\lvert{x-y}\rvert^\alpha$ for any $x,y\in\mathbb{R}$. Then
$\max_{x\in[0,1]}\lvert{S_N(x)-f(x)}\rvert \rightarrow 0 \quad (N\rightarrow\infty).$

2. Moreover if $\alpha>1/2$, then
$\sum_{n=-\infty}^\infty\lvert{c_n}\rvert<\infty.$

3. There exists a is -periodic continuous function $f(x)$ such that 
$S_N(0) \rightarrow \infty \quad (N\rightarrow\infty).$

Here we have some remarks on Hoelder continuity. The condition $\lvert{f(x)-f(y)}\rvert \leqq L\lvert{x-y}\rvert^\alpha$ is a restriction only for small $\lvert{x-y}\rvert$ since the absolute value of a periodic continuous function has a maximum. Roughly speaking, the exponent $\alpha$ expresses that $f(x)$ is differentiable of order $\alpha$. When $\alpha=1$, the Hoelder continuity is said to be the Lipshitz continuity. The Rademach theorem shows that the Lipshitz continuous function is differentiable almost everywhere with respect to the Lebesgue measure on $\mathbb{R}$. Roughly speaking, Part 1 and Part 3 assert that if $f(x)$ is a little bit smoother than continuous functions, its Fourier series uniformly converges to , and otherwise, its Fourier series does not necessarily converge even pointwisely. 

We show the outline of the proof of Part 1. Note that $S_N(x)-f(x)$ can be written as 

$S_N(x)-f(x)
=
\int_{-1/2}^{1/2}
D_N(t)\{f(x+t)-f(x)\}
dt.$

The Hoelder condition implies that for any small $\delta>0$

$S_N(x)-f(x)
=
\int_{\lvert{t}\rvert<\delta}
\mathcal{O}(\lvert{t}\rvert^{-1+\alpha})
dt
+
\int_{\delta<\lvert{t}\rvert<1/2}
\sin\bigl(\pi(2N+1)t\bigr)
\times\dotsb
dt.$

The second term of the right hand side is vanishing as $N\rightarrow\infty$ due to the oscillation and $\displaystyle\int_0^1\sin(2\pi t)dt=0$. Indeed for any $\varphi(t) \in C^1[0,1]$, 

$\begin{aligned}
  \int_0^1\sin(2\pi Nt)\varphi(t)dt
& =
  \left[-\frac{\cos(2\pi Nt)\varphi(t)}{2\pi N}\right]_0^1
  +
  \int_0^1\frac{\cos(2\pi Nt)\varphi^\prime(t)}{2\pi N}dt
\\
& =
  \mathcal{O}(N^{-1})
  \quad
  (N\rightarrow\infty).
\end{aligned}$

The differentiability of  is not required for this. In fact one can also prove 

$\int_0^1\sin(2\pi Nt)\varphi(t)dt \rightarrow 0 \quad (N\rightarrow\infty)$

for $\varphi(t) \in C[0,1]$ by using the uniform continuity of $\varphi(t)$ and  $\displaystyle\int_0^1\sin(2\pi t)dt=0$. So we have for any small $\delta>0$

$S_N(x)-f(x)
=\mathcal{O}(\delta^\alpha)+o(1)
\quad
(N\rightarrow\infty).$

"""

# ╔═╡ e956146c-d9f1-4c0f-8422-32c29e6d5b69
md"""
#### 6. Fourier series of triangular function
Consider a $1$-periodic function  defined by 

$f(x)
:=
\begin{cases}
x, &\ x\in[0,1/2),
\\
1-x, &\ x\in[1/2,1).
\end{cases}$

This is the same as the example in Section 3. Here we observe the uniform convergence of the Fourier series of $f(x)$, which is 

$\begin{aligned}
  f(x)
& =
  \frac{1}{4}
  -
  \sum_{k=1}^\infty
  \frac{2}{(2k-1)^2\pi^2}
  \cos\bigl(2(2k-1)\pi x\bigr),
\\
  S_{2K-1}[f](x)
& =
  \frac{1}{4}
  -
  \sum_{k=1}^K
  \frac{2}{(2k-1)^2\pi^2}
  \cos\bigl(2(2k-1)\pi x\bigr).
\end{aligned}$

"""

# ╔═╡ b13742b2-45ca-42ee-bdcd-649b036da15f
begin
    M1=20;
	x = range(-0.1, 1.1, length = 481);
	
	f = zeros(481);
	for l=1:481
	f[l] = min(x[l]-Float64(floor(x[l])),1-x[l]-Float64(floor(1-x[l])));
	end
	
	S1=ones(M1+1,481)/4;
    for k=2:M1+1
		for l=1:481
		    S1[k,l]=S1[k-1,l]-2*cos(2*(2*k-3)*pi*x[l])/(2*k-3)^2/pi^2;
		end
	end
end

# ╔═╡ c05b1f88-2750-4f76-8af0-10091651289a
md"""
 $K$ = $(@bind N1 Slider(0:M1, show_value=true))
"""

# ╔═╡ cb7d861a-8d69-4e08-a69d-7638de638e3f
begin
    plot(x,S1[N1+1,:],
		 grid=false,
		 linewidth=2,
		 xlim=(-0.15,1.15),
		 ylim=(-1/16,9/16),
		 title="Triangular Function and its Fourier Series",
		 xticks = ([0 1/2 1;], [0,1/2,1]),
		 yticks = ([0 1/2;], [0,1/2]),
		 xlabel="\$x\$",
		 label="\$S_{2K-1}[f](x)\$",
		 legend=:topright,
		 legendfont=font(12))
	plot!(x,f,
		  linewidth=2,
	      label="\$f(x)\$")
end

# ╔═╡ 729497d1-2a1e-4fd3-bef4-fd3280ffd673
md"""
#### 7. A continuous function whose fourier series does not converges even pointwisely

Define a function

$f(x)
=
\sum_{k=1}^\infty
\frac{2}{3^k}
\sin\{2\pi(2k+1)n(k)x\}
\sum_{l=1}^{n(k)}
\frac{\sin(2\pi lx)}{l},$

$n(k)=3^{3^k},
\quad\text{i.e.,}\quad
\log{n(k)}=3^k\log{3},$

which converges uniformly on $\mathbb{R}$, and is a $1$-periodic continuous function. 
Consider its Fourier series and the patrial sum

$f(x) \sim \sum_{n=-\infty}^\infty C_ne^{2\pi nx}, 
\quad
S_N(x):=\sum_{n=-N}^N C_ne^{2\pi nx}.$

It is well-known that the sequence $\{S_N(0)\}$ does not converge, more precisely, we have 

$S_{(2k+1)n(k)}(0)-S_{2kn(k)}(0)
=
\frac{1}{3^k}\sum_{j=1}^{n(k)-1}\frac{1}{j}
>
\log{3}.$

We tried to see the graph of $f$ by using 

$f_K(x)
:=
\sum_{k=1}^K
\frac{2}{3^k}
\sin\{2\pi(2k+1)n(k)x\}
\sum_{l=1}^{n(k)}
\frac{\sin(2\pi lx)}{l}.$

Unfortunately, we could not do it well since 

$n(1)=9, \quad n(2)=19683, \quad n(3)=7625597484987.$

We could compute $f_2$, and not $f_3$, and show the graph of $f_2$ below.
"""

# ╔═╡ adf3292f-4b05-4779-a316-dfe9662c055a
begin
	x5 = range(-0.6, 0.6, length = 1201);

	f5=zeros(length(x5))
	for p=1:length(x5)
		for k=1:2
			g5=0.0
			for l=1:3^(3^k)
			    g5 += sin(2*pi*l*x5[p])
		    end
		f5[p] += (2/(3^k))*sin(2*pi*(2*k+1)*3^(3^k)*x5[p])*g5
		end
	end
end

# ╔═╡ 9956a3d2-9ae9-4ef4-9f1c-55064001a025
begin
    plot(x5,f5,
		 grid=false,
		 linewidth=2,
		 xlim=(-0.6,0.6),
		 title="The continuous function \$f_2\$",
		 xticks = ([-0.5 0 0.5;], [-0.5,0,0.5]),
		 yticks = ([-20 0 20;], [-20,0,20]),
		 xlabel="\$x\$",
		 legend=false,
		 legendfont=font(12))
end

# ╔═╡ 596f1b0f-7a83-4f48-9316-17ef15d88cbe
md"""
#### 8. Step function and sawtooth function
Let $Y(x)$ be the heaviside function, that is, 

$Y(x)
=
\begin{cases}
1, &\ x\geqq0,
\\
0, &\ x<0.
\end{cases}$

The value $Y(0)=1$ does not matter essentially, and is not needed to be defined. Consider a step function $g(x)$ and a sawtooth function $h(x)$ defined by 

$g(x)
:=
\sum_{n=-\infty}^\infty
Y(x+n) \times Y(1/2-x-n),
\quad
h(x):=x-[x].$

They are $1$-periodic function such that for $x\in[0,1)$

$g(x)
=
\begin{cases}
1, &\ x\in[0,1/2],
\\
0, &\ x\in(1/2,1),
\end{cases}
\quad
h(x)=x.$

Their Fourier series are 

$\begin{aligned}
  g(x)
& \sim
  \frac{1}{2}
  +
  \sum_{k=1}^\infty
  \frac{2}{(2k-1)\pi}
  \sin\bigl(2(2k-1)\pi x\bigr),
\\
  S_{2K-1}[g](x)
& =
  \frac{1}{2}
  +
  \sum_{k=1}^K
  \frac{2}{(2k-1)\pi}
  \sin\bigl(2(2k-1)\pi x\bigr),
\\
  h(x)
& \sim
  \frac{1}{2}
  -
  \sum_{n=1}^\infty
  \frac{1}{n\pi}
  \sin(2\pi nx),
\\
  S_N[h](x)
& =
  \frac{1}{2}
  -
  \sum_{n=1}^N
  \frac{1}{n\pi}
  \sin(2\pi nx).
\end{aligned}$

On one hand, at the discontinuous points 

$\begin{aligned}
  S_{2K-1}[g](k)
  =
  \frac{1}{2}
& \rightarrow
  \frac{1}{2}
  =
  \frac{g(k-0)+g(k+0)}{2}
  \quad
  (K\rightarrow\infty),
\\
  S_{2K-1}[g](k+1/2)
  =
  \frac{1}{2}
& \rightarrow
  \frac{1}{2}
  =
  \frac{g(k+1/2-0)+g(k+1/2+0)}{2}
  \quad
  (K\rightarrow\infty),
\\
  S_N[h](k)
  =
  \frac{1}{2}
& \rightarrow
  \frac{1}{2}
  =
  \frac{h(k-0)+h(k+0)}{2}
  \quad
  (N\rightarrow\infty)
\end{aligned}$

for all $k\in\mathbb{Z}$. On the other hand, we can prove that for any small $\delta>0$

$\begin{aligned}
  \max_{x\in[\delta,1/2-\delta]\cup[1/2+\delta,1-\delta]}
  \lvert{S_{2K-1}[g](x)-g(x)}\rvert
& \rightarrow
  0
  \quad
  (K\rightarrow\infty),
\\
  \max_{x\in[\delta,1-\delta]}
  \lvert{S_N[h](x)-h(x)}\rvert
& \rightarrow
  0
  \quad
  (N\rightarrow\infty).
\end{aligned}$

The partial sums oscillate violately as $N$ increases near the discontinuous points. This is called the Gibbs Phenomenon.
"""

# ╔═╡ 1bdd1071-c4a0-474a-a394-cd5ff77d5481
begin
    M2=80;
	
	g = zeros(481)
	for l=41:241
	    g[l] = 1;
	end
	for l=441:481
		g[l]=1;
	end
	
	S2=ones(M2+1,481)/2;
    for k=2:M2+1
		for l=1:481
		    S2[k,l]=S2[k-1,l]+2*sin(2*(2*k-3)*pi*x[l])/(2*k-3)/pi;
		end
	end
end

# ╔═╡ e4f95f78-cd34-4ce6-96f4-65c0151c55cc
md"""
 $K$ = $(@bind N2 Slider(0:M2, show_value=true))
"""

# ╔═╡ 2c0a572e-1da4-40cb-a5d4-ff510175d7fb
begin
    plot(x,S2[N2+1,:],
		 grid=false,
		 linewidth=2,
		 xlim=(-0.15,1.15),
		 ylim=(-1/8,9/8),
		 title="Step Function and its Fourier Series",
		 xticks = ([0 1/2 1;], [0,1/2,1]),
		 yticks = ([0 1/2 1;], [0,1/2,1]),
		 xlabel=L"x",
		 label=L"S_{2K-1}[g](x)",
		 legend=:topright,
		 legendfont=font(12))
	plot!(x,g,
		  linewidth=2,
	      label=L"g(x)")
end

# ╔═╡ 83ee20aa-e518-442b-af61-d1473d2d3d1a
begin
    h = zeros(481);
	for l=1:481
		h[l]=x[l]-Float64(floor(x[l]));
	end
		
	S3=ones(M2+1,481)/2;
    for k=2:M2+1
		for l=1:481
		    S3[k,l]=S3[k-1,l]-sin(2*(k-1)*pi*x[l])/(k-1)/pi;
		end
	end
end

# ╔═╡ 7cabd643-b904-416f-9382-7a4703e19b22
md"""
 $N$= $(@bind N3 Slider(0:M2, show_value=true))
"""

# ╔═╡ 36c46073-1425-43fc-93e8-fea1b2e1f192
begin
    plot(x,S3[N3+1,:],
		 grid=false,
		 linewidth=2,
		 xlim=(-0.15,1.15),
		 ylim=(-1/8,9/8),
		 title="Sawtooth Function and its Fourier Series",
		 xticks = ([0 1/2 1;], [0,1/2,1]),
		 yticks = ([0 1/2 1;], [0,1/2,1]),
		 xlabel=L"x",
		 label=L"S_N[h](x)",
	     legend=:top,
		 legendfont=font(12))
	plot!(x,h,
		  linewidth=2,
	      label=L"h(x)")
end

# ╔═╡ 47454abd-5b71-465e-bb2a-a4e3c544d350
md"""
#### 9. Decomposition of piecewise smooth functions
Finally we consider piecewise smooth $1$-periodic functions $f(x)$ on $\mathbb{R}$, that is, 

-  $f(x+1)=f(x)$ for all $x\in\mathbb{R}$.

- There exist numbers $0=a_0<a_1<\dotsb<a_m=1$ such that $f(x)$ is continuously differentiable on each open interval $(a_{j-1},a_j$, $j=1,\dotsc,m$, and limits exist:
$f(a_0+0), \quad f(a_1-0), \quad f(a_1+0), \quad \dotsc \quad f(a_m-0),$ 
$f^\prime(a_0+0), \quad f^\prime(a_1-0), \quad f^\prime(a_1+0), \quad \dotsc \quad f^\prime(a_m-0).$ 

We remark that if we restrict  on each closed interval $[a_{j-1}+k,a_j+k]$ for $j=1,\dotsc,m$ and $k\in\mathbb{Z}$, $f(x)$ can be regarded as a continuously differentiable function. We can split $f(x)$ into a Lipshitz continuous function, some step functions, and a sawtooth function. We shall illustrate this for a simplest example. 

Let $f(x)$ be a piecewise smooth $1$-periodic functions on $\mathbb{R}$ with discontinuities at most $0,c,1$ in the interval $[0,1]$ with some $c\in(0,1)$. Note that $f(0\pm0)=f(1\pm0)$. We denote by $\chi_{[c,1]}(x)$ the characteristic function of the interval $[c,1]$. We define $f_0(x)$, $f_1(x)$, and $f_2(x)$ by 

$\begin{aligned}
  f_0(x)
& :=
  f(x)-f_1(x)-f_2(x),
\\
  f_1(x)
& :=
  \bigl(f(c+0)-f(c-0)\bigr)\chi_{[c,1]}(x),
\\
  f_2(x)
& :=
  -
  \bigl\{
  \bigl(f(0+0)-f(0-0)\bigr)
  +
  \bigl(f(c+0)-f(c-0)\bigr)
  \bigr\}
  h(x).
\end{aligned}$

Then $f_0(x)$ is Lipshitz continuous. So the piecewise smooth function $f(x)$ is splitted into the Lipshitz function part $f_0(x)$, the step function part $f_1(x)$, and the sawtooth function part $f_2(x)$. From the view point of Fourier series, $f_0(x)$ is harmless, and $f_1(x)+f_2(x)$ contains discontinuities where the Gibbs phenomenon occurs. 

We illustrate this decomposition visually. 

"""

# ╔═╡ 252f9cc5-5c34-49f1-96a7-bf5108f30fbd
begin
	Z=100*ones(481);
    f0=sin.(2*pi*x);
	f1=(ones(481)-g)/2;
	f2=h/3;
	F=f0+f1+f2;

	plot(x,F,
		 grid=false,
		 linewidth=2,
		 xlim=(-0.15,1.15),
		 ylim=(-1.1,1.7),
		 title="A piecewise smooth function",
		 xticks = ([0 1/2 1;], [0,L"c",1]),
		 xlabel=L"x",
		 label=L"f(x)",
		 legend=:topright,
		 legendfont=font(12))
end

# ╔═╡ d8503bba-e3f3-4f9b-ac4e-bc4f1ffd7fbf
begin
	plot(x,F,
		 grid=false,
		 linewidth=2,
		 xlim=(-0.15,1.15),
		 ylim=(-1.1,1.7),
		 title="Decomposition of the piecewise smooth function",
		 xticks = ([0 1/2 1;], [0,L"c",1]),
		 xlabel=L"x",
		 label=L"f(x)",
		 legend=:top,
		 legendfont=font(12))
	plot!(x,Z,
		  linewidth=2,
	      label=:false)
	plot!(x,Z,
		  linewidth=2,
	      label=:false)
	plot!(x,Z,
		  linewidth=2,
	      label=:false)
	plot!(x,f0,
		  linewidth=2,
		  color="magenta",
	      label=L"f_0(x)")
	plot!(x,f1,
		  linewidth=2,
		  color="green",
	      label=L"f_1(x)")
	plot!(x,f2,
		  linewidth=2,
		  color="orange",
	      label=L"f_2(x)")
end

# ╔═╡ ed6e4c4d-a686-4271-b964-0fe31c0361be
md"""
#### 10. Weyl's equidistribution theorem

Let $\lfloor{x}\rfloor$ be the floor function for $x\in\mathbb{R}$, that is, $\lfloor{x}\rfloor:=\max\{m \in \mathbb{Z} : m \leqq x\}$, and set 

$\langle{x}\rangle:=x-\lfloor{x}\rfloor.$

Then $0\leqq\langle{x}\rangle<1$ and $x\equiv\langle{x}\rangle$ mod $1$.

We say that a sequence $\{a_n\}_{n=1}^\infty \subset [0,1)$ is equidistributed in $[0,1)$ if 

$\frac{\sharp\{n=1,\dotsc,N : a<a_n<b\}}{N} \rightarrow b-a \quad (N \rightarrow \infty)$

**Weyl's equidistribution theorem**: For any $\gamma\in\mathbb{R}\setminus\mathbb{Q}$, $\{\langle{n\gamma}\rangle\}_{n=1}^\infty$ is equidistributed in $[0,1)$.

We observe Weyl's equidistribution theorem using a bijection 

$[0,1)\ni\langle{x}\rangle \mapsto e^{2\pi{i}\langle{x}\rangle}=e^{2\pi{i}x} \in \mathbb{S}^1:=\{z\in\mathbb{C} : \lvert{z}\rvert=1\}$.
"""

# ╔═╡ bb9a5642-c3ec-46e6-8dc5-f1100d3a1ed4
MWeyl=250;

# ╔═╡ 95c91522-88db-4f8f-9c2f-4e2e9cab84b6
md"""
 $N$ = $(@bind N10 Slider(1:MWeyl, show_value=true))

 $\gamma^2$ = $(@bind g2 Slider(2:15, show_value=true))
"""

# ╔═╡ e6c5f708-b224-42de-996c-e7ed7b3e9125
begin
	t10 = range(0, 1, length = 101);
	c101=cos.(2*pi*t);
	c102=sin.(2*pi*t);

    γ=sqrt(g2);
	x10=zeros(MWeyl);
	y10=zeros(MWeyl);

	for n=1:MWeyl
		x10[n]=cos(2*pi*γ*n);
		y10[n]=sin(2*pi*γ*n);
	end

end

# ╔═╡ 85863bee-e916-4d1e-98ab-52c6bc080fac
begin
    scatter(x10[1:N10],y10[1:N10],
		    title="Weyl's Equidistribution Theorem",
	        grid=false,
	        xlim=(-1.1,1.1),
	        ylim=(-1.1,1.1),
	        aspect_ratio=1.0,
	        xticks=false,
	        yticks=false,
	        label=false,
	        xaxis=false,
	        yaxis=false,
	        annotations=[(0,0.1,L"e^{2πiγn}, γ=\sqrt{%$g2}",20),(0,-0.15,L"n=1,...,%$N10",20)],   
		    annotationguide=:auto)
	plot!(c101,c102,label=false,linecolor=:magenta)
end

# ╔═╡ 75e78914-4f93-4f98-8137-6512ba96e868
md"""
#### 11. Continuous but nowhere differentiable functions

For any $\alpha \in (0,1)$ we set 

$f_\alpha(x):=\sum_{n=0}^\infty\frac{\exp(2\pi i2^nx)}{2^{n\alpha}}.$

The right hand side of the above converges uniformly since $2^\alpha>1$, and $f_\alpha$ becomes a $1$-periodic continuous function on $\mathbb{R}$. 
To show the graph of $f_\alpha$ below, set

$S_{\alpha,N}(x):=\sum_{n=0}^N\frac{\exp(2\pi i2^nx)}{2^{n\alpha}}.$

"""

# ╔═╡ 7bcc6226-cfcf-4e27-b70a-5e88bbf271aa
M6=20;

# ╔═╡ 76b7bc7b-d29e-43b8-bd38-0772479c5fda
md"""
 $N$ = $(@bind N6 Slider(0:M6, show_value=true, default=0)) 

 $\alpha$ = $(@bind α Slider(0.4:0.1:0.9, show_value=true, default=0.5))
"""

# ╔═╡ 90937af9-9230-4676-8b6c-34559a96d635
begin
	x6 = range(-0.1, 1.1, length = 121);
	
	f6 = zeros(length(x6),M6+1);
	g6 = zeros(length(x6),M6+1);
	for l=1:length(x6)
		f6[l,1]=cos(2*pi*x6[l]);
		g6[l,1]=sin(2*pi*x6[l]);
		for n=2:M6+1
		    f6[l,n] = f6[l,n-1] + cos(2*pi*(2^n)*x6[l])/(2^(α*n));
		    g6[l,n] = g6[l,n-1] + sin(2*pi*(2^n)*x6[l])/(2^(α*n));
		end
	end
	
end

# ╔═╡ 7f3c0c03-a149-49b6-bbe8-a57556dce69d
begin
    plot(x6,f6[:,N6+1],
		 grid=false,
		 linewidth=2,
		 ylim=(-2,3.5),
		 title="Continuous but nowhere differentiable functions",
		 xticks = ([0 1/2 1;], [0,1/2,1]),
		 #yticks = ([0 1/2 1;], [0,1/2,1]),
		 xlabel=L"x",
		 label="real part",
		 legend=:top,
		 legendfont=font(12))
	plot!(x6,g6[:,N6+1],linewidth=2,label="imaginary part",linecolor=:magenta)
end

# ╔═╡ 69f2d10c-38bf-43c5-a7dc-3ac70f5b9763
md"""
## Fourier series and Fourier transform
"""

# ╔═╡ 99a69d04-aed2-436c-85f9-5f0ca3fb1ace
md"""
#### 12. The Nyquist–Shannon sampling theorem and aliasing

We state the Nyquist-Shannon sampling theorem and observe an example. This is basically concerned with one-variable integrable functions whose Fourier transforms are compactly supported. Such functions are automatically real analytic and can be extended to entire functions on $\mathbb{C}$. We employ the definition of Fourier transform as follows. For an appropriate one-variable function $f(x)$, its Fourier transform is defined by 

$\hat{f}(\xi)
:=
\frac{1}{\sqrt{2\pi}}
\int_{-\infty}^\infty
e^{-ix\xi}
f(x)
dx,
\quad
\xi\in\mathbb{R}.$

We use the normalized sinc function defined by 

$\operatorname{sinc}(x)
:=
\frac{\sin(\pi x)}{\pi x}
\quad
(x\ne0),
\quad
\operatorname{sinc}(0)
:=1.$
"""

# ╔═╡ 445ccce6-7439-45ad-8724-98bb8a9ca4be
md"""
#### Theorem (the Nyquist-Shannon sampling theorem)

Let $B>0$. Suppose that $f(x)$ is an integrable smooth functions on $\mathbb{R}$ and $\hat{f}(\xi)=0$ for $\lvert\xi\rvert \geqq B$. Then for any $s\in(0,\pi/B]$, we have 

$f(x)
=
\sum_{n=-\infty}^\infty
f(sn)\cdot\operatorname{sinc}\left(\frac{x}{s}-n\right),$


$\int_{-\infty}^\infty
\lvert{f(x)}\rvert^2
dx
=
s
\sum_{n=-\infty}^\infty
\lvert{f(sn)}\rvert^2.$

The parameter $s>0$ is said to be the sampling period. 

"""

# ╔═╡ 21d71d4d-fa84-4dd7-bb86-dc5aa6b07cad
md"""
#### Observation: sampling and aliasing

Let $f(x)$ be a function defined by 

$f(x)
:=
\frac{2\bigl(1-\cos(\pi x)\bigr)}{(\pi x)^2}
\quad
(x\ne0),
\qquad
f(0):=1.$

Then 

$\hat{f}(\xi)
=
\frac{\sqrt{2}(\pi-\lvert\xi\rvert)}{\pi^{3/2}}
\quad
(\lvert\xi\rvert<\pi),
\qquad
\hat{f}(\xi)=0
\quad
(\lvert\xi\rvert\geqq\pi).$

If we set 

$F(x;s):=
\sum_{-\infty}^\infty
f(sn)
\cdot
\operatorname{sinc}\left(\frac{x}{s}-n\right),
\quad
s>0,$

then the Nyquist–Shannon Sampling Theorem implies that 

$f(\cdot)=F(\cdot;s), \quad s\in(0,1],$

$f(\cdot){\ne}F(\cdot;s), \quad s>1,$

The latter one is said to be aliasing. 

"""

# ╔═╡ 373db211-b7cf-4f64-ac04-5b97fac9186b
md"""
s = $(@bind s Slider(0.7:0.1:2, show_value=true))
"""

# ╔═╡ 2283f0b6-58db-4f16-8c69-85d7c5a3eca8
begin
	x11 = range(-10, 10, length = 201);
	F11=sinc.(x11/s);
	f11=ones(201);
	for k=1:100
		f11[k]=2*(1-cos(pi*x11[k]))/(pi*x11[k])^2;
		f11[k+101]=2*(1-cos(pi*x11[k+101]))/(pi*x11[k+101])^2;
	end
	
	for n=1:20
		a=2*(1-cos(pi*n*s))/(pi*n*s)^2;
        b=sinc.(x11/s-n*ones(201));
        c=sinc.(x11/s+n*ones(201));
        F11=F11+a*b+a*c;
	end
end

# ╔═╡ 048c94e7-31fa-41e4-9336-70b44320378b
begin
    plot([F11,f11],
		 grid=false,
		 linewidth=2,
		 ylim=(-0.3,1.2),
         title="Nyquist-Shannon sampling and aliasing, s=$s",
         xticks = ([0 51 101 151 201;], [-10,-5,0,5,10]),
		 xlabel=L"x",
         yticks = ([0 0.5 1;], [0,0.5,1]),
         label=[L"F(x;s)" L"f(x)"],
	     legend=:topright,
		 legendfont=font(12))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Colors = "~0.12.10"
LaTeXStrings = "~1.3.0"
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "42d370b6ae201137fc5f9c1e364cc913578b870f"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccc81ba5e42497f4e76553a5545665eed577a663"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.0+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "50c11ffab2a3d50192a228c313f05b5b5dc5acb2"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "386b47442468acfb1add94bf2d85365dea10cbab"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "064b532283c97daae49e544bb9cb413c26511f8c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.8"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5cb3c5d039f880c0b3075803c8bf45cb95ae1e91"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.51+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─581f8ff0-d1fc-11ec-3cab-5f335ec587ad
# ╟─25063436-79e0-4f01-ab92-f01defb4e948
# ╟─79521e83-f32a-48aa-bdf1-27f3f95ecb9e
# ╟─c299148e-4b85-4bf3-a514-20d1dff0dc59
# ╟─6f6479b6-265a-4474-8a90-3df0ff98cfac
# ╟─fd9270b6-bcb1-444a-b2e8-bc07efee3a24
# ╟─a94cb318-c3e9-4285-af58-154d6955f2d5
# ╟─3edf319a-48e4-48f0-82ab-f81e5f2c581c
# ╟─a1aa2ce7-db42-44d2-a9d8-4ae67976f76f
# ╟─cceffb37-cb1f-4436-aab6-17f325504276
# ╟─897a6954-9e77-41d8-a88d-d43b4c9902c6
# ╟─e956146c-d9f1-4c0f-8422-32c29e6d5b69
# ╟─b13742b2-45ca-42ee-bdcd-649b036da15f
# ╟─c05b1f88-2750-4f76-8af0-10091651289a
# ╟─cb7d861a-8d69-4e08-a69d-7638de638e3f
# ╟─729497d1-2a1e-4fd3-bef4-fd3280ffd673
# ╟─adf3292f-4b05-4779-a316-dfe9662c055a
# ╟─9956a3d2-9ae9-4ef4-9f1c-55064001a025
# ╟─596f1b0f-7a83-4f48-9316-17ef15d88cbe
# ╟─1bdd1071-c4a0-474a-a394-cd5ff77d5481
# ╟─e4f95f78-cd34-4ce6-96f4-65c0151c55cc
# ╟─2c0a572e-1da4-40cb-a5d4-ff510175d7fb
# ╟─83ee20aa-e518-442b-af61-d1473d2d3d1a
# ╟─7cabd643-b904-416f-9382-7a4703e19b22
# ╟─36c46073-1425-43fc-93e8-fea1b2e1f192
# ╟─47454abd-5b71-465e-bb2a-a4e3c544d350
# ╟─252f9cc5-5c34-49f1-96a7-bf5108f30fbd
# ╟─d8503bba-e3f3-4f9b-ac4e-bc4f1ffd7fbf
# ╟─ed6e4c4d-a686-4271-b964-0fe31c0361be
# ╟─bb9a5642-c3ec-46e6-8dc5-f1100d3a1ed4
# ╟─95c91522-88db-4f8f-9c2f-4e2e9cab84b6
# ╟─85863bee-e916-4d1e-98ab-52c6bc080fac
# ╟─e6c5f708-b224-42de-996c-e7ed7b3e9125
# ╟─75e78914-4f93-4f98-8137-6512ba96e868
# ╟─7bcc6226-cfcf-4e27-b70a-5e88bbf271aa
# ╟─76b7bc7b-d29e-43b8-bd38-0772479c5fda
# ╟─7f3c0c03-a149-49b6-bbe8-a57556dce69d
# ╟─90937af9-9230-4676-8b6c-34559a96d635
# ╟─69f2d10c-38bf-43c5-a7dc-3ac70f5b9763
# ╟─99a69d04-aed2-436c-85f9-5f0ca3fb1ace
# ╟─445ccce6-7439-45ad-8724-98bb8a9ca4be
# ╟─21d71d4d-fa84-4dd7-bb86-dc5aa6b07cad
# ╟─373db211-b7cf-4f64-ac04-5b97fac9186b
# ╟─048c94e7-31fa-41e4-9336-70b44320378b
# ╟─2283f0b6-58db-4f16-8c69-85d7c5a3eca8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
