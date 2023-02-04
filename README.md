# Fast Lower Bounds of Functions on Finite Abelian Groups

This is the a tool to compute the lower bound of function on finite abelian groups, which is based on SOS relaxation.  

## Dependencies
- Matlab 2016b
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)

## Usage
The Numerical experiments of this paper can be performed by "Numerical_experiments.m"

### construct a function on fintie abelian group  
$$f:G \mapsto \mathbb{Z},~G=\prod_{i=1}^{k}  \mathbb{C}_{n_i},$$
with 

$$f=\sum_{j=1}^{s}a_jx^{\alpha_j}, ~x_i \in   \mathbb{C}_{n_i}$$
Let n be a 1-dim array such that n(i)=n_i, then set
```
f=CZ(n);
for j=1:s
f(alpha_j)=a_j;
end
```
### compute fast lower bound  
If already known that the set $S$, such that $\operatorname{Im}(f)$ is a subset of $S$, to compute the lower bound with input parameters d and k, just by:
```
[Q,lb,Index]=FastLowerBound_Poly(f,d,S,k)
```
the output $Q$ is the Gram matrix with rows and columns indexed by Index,and lb is the lower bound.
### rounding
For a cnf/wcnf file, to use
```
[Q,f,Index,y_Gram,y_Moment,lb]=OurRounding(File,k);
```
to get the rounding reslut
y_Gram and y_Moment is the rounding results by null space of Gram matrix and rank-1 approx of Moment matrix.
