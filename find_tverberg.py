import cvxpy as cp
import numpy as np

r = 10 # number of partitions (colors)
d = 2 # space dimension
m = (r-1)*(d+1)+1 # number of points
np.random.seed(4)
P = np.random.randn(m,d) # point set
print_output = True


V = np.ones((m,d+1)); V[:,0:d] = P
W = np.zeros((r,r-1)); W[0:(r-1),:] = np.eye(r-1); W[-1,:] = -1
VxW = V.reshape( (1,m,1,d+1) ) * W.reshape( (r,1,r-1,1) )

S = np.transpose( VxW.reshape( (r,m,m-1) ), (1,0,2) ) # S.shape = (r,m-1,d)

# Now we find colorful Carathéodory set s from S using the Barany's Algorithm
# Barany, I. (1982). A generalization of Carathéodory’s theorem. Discrete 
# Mathematics, 40(23), 141152. https://doi.org/10.1016/0012-365X(82)90115-7

I = np.zeros((m,), np.int32);
min_sqr_norm = 1
sqr_norm_tol = 1e-12 # close enough to be considered unchanged
min_sqr_norm_tol = 1e-16 # close enough to be considered zero
residual_history = np.zeros((0,))
while min_sqr_norm > min_sqr_norm_tol:
	s = np.zeros( (m-1,m) ) # the colorful set
	for i in range(m):
	   s[:,i] = S[i,I[i],:]

	## here we solve a constrained least squares minimization to find the min-norm point:
	x = cp.Variable(m)
	objective = cp.Minimize(cp.sum_squares(s*x))
	constraints = [0 <= x, x <= 1, sum(x) == 1]
	prob = cp.Problem(objective, constraints)
	min_sqr_norm = prob.solve( eps_abs = 1.0e-10, eps_rel = 1.0e-10, verbose=False)
	coeff = x.value
	min_norm_pt = s.dot(coeff)

	# we then pivot out a point that does not increase the min_sqr_norm
	pivot_order = coeff.argsort()
	x2 = cp.Variable(m-1)

	for pivot in pivot_order:
		not_pivot = np.setdiff1d(range(m), pivot, assume_unique=True)
		s2 = s[:,not_pivot]
		objective2 = cp.Minimize(cp.sum_squares(s2*x2))
		constraints2 = [0 <= x2, x2 <= 1, sum(x2) == 1]
		prob2 = cp.Problem(objective2, constraints2)
		min_sqr_norm2 = prob2.solve( eps_abs = 1.0e-10, eps_rel = 1.0e-10, verbose=False)
		if  np.abs( min_sqr_norm - min_sqr_norm2 ) < sqr_norm_tol:
			prev_pt_id = I[pivot]
			t = S[pivot,:,:] 
			not_prev_pt_id = np.setdiff1d(range(r), prev_pt_id, assume_unique=True)
			pt_id = not_prev_pt_id[t.dot(min_norm_pt)[not_prev_pt_id].argmin()]
			I[pivot] = pt_id
			break
	if print_output:
		print("residual =",min_sqr_norm)
	residual_history = np.append(residual_history, min_sqr_norm)

# create incidence matrix from I
color_bool = np.zeros((r, m), dtype=bool)
for i in range(r):
    color_bool[i,I==i] = True;

# compute the center-point
center_point = np.array(P).transpose().dot(coeff)


# %%%%%%%
# output
print("center_point =",center_point)
print("incidece_mat =\n",color_bool+0)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
y = residual_history

ax.semilogy(y, 'o-')

ax.set(xlabel='iteration', ylabel='residual',
       title='residual history')
# ax.grid()
plt.show()