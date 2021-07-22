"""
Implements the string class with performs calculations
on the string object
"""
import numpy as np

class String():
    """
    Class which operates the string itself as a series of 
    replicas in CV space. Note that the string only maintains
    its current value of CVs and not the history, so calling 
    update will erase the previous values.
    """
    def __init__(self, 
        z,
        k,
        dt,
        smooth,
        grad_F=None, 
        M=None
    ):
        """
        Paramters
        ---------
        z (np.ndarray(R, N)): matrix of CVs (R, N) in shape
        R (int): number of replicas
        N (int): knumber of CVs
        k (np.ndarray(N)): array force constant used for each CV
        dt (int): forward Euler timestep
        smooth (float): smoothing parameter between 0 and 1
        grad_F (np.ndarray(R, N)): mean force for each replica (gradient)
        M (np.ndarray(R, N, N)): metric tensor at each replica of size N x N
        """
        self.z = z
        self.R, self.N = self.z.shape
        self.k = k
        self.dt = dt
        self.smooth = smooth
        self.grad_F = grad_F
        self.M = M

    def compute_grad_F(self, z_curr):
        """
        Compute grad(F) using eq. (41) in Maragliano et al. 2006

        Parameters
        ----------
        z_curr (np.ndarray(R, N))
        """
        self.grad_F = (z_curr - np.mean(cvs, axis=0)) * self.k
    
    def compute_M(self, z_curr, z_deriv):
        """
        Compute M tensor using eq. (42) in Maragliano et al. 2006

        Parameters
        ----------
        z_curr:
        z_deriv:
        """
        #Getting the derivatives, this is only of pwd
        dthdx=[]
        for i in range(num_CV):
            lines=np.genfromtxt(filename+str(i))
            dthdx.append(CVder[CVinfo[i][1]](lines,CVinfo[i]))
        #Getting tensor M, this should be universal.
        M=np.zeros((num_CV,num_CV))
        for i in range(num_CV):
            for j in range(num_CV):
                for ki,k_i in enumerate(CVinfo[i][2][:2]):
                    for kj,k_j in enumerate(CVinfo[j][2][:2]):
                        if k_i==k_j:
                            for r in range(3):
                                M[i,j]+=np.dot(dthdx[i][r+3*ki,:],dthdx[j][r+3*kj,:])/(len(dthdx[i][r+3*kj,:])*masses[int(k_i)])

    
    def forward_euler(self):
        """
        Evolve images forward one time step using the discretized eq. (46) in 
        Maragliano et al. 2006

        Parameters
        ---------

        """
        z_new = np.zeros_like(self.z)
        
        # first and last images are constant
        z_new[0, :] = self.z[0, :]
        z_new[-1, :] = self.z[-1, :]

        # compute M \cdot \nabla_z F
        M_dot_grad_F = self.M @ self.grad_F

        for m in range(1, self.R - 1):
            #checking sign of sigma
            if np.dot((zcurrent[m + 1, :] - zcurrent[m, :]), M_dot_grad_F) >= 0:
                sigma = 1
            else:
                sigma = -1

            # projection operator, eq. (47)
            dz = (zcurrent[m + sigma, :] - zcurrent[m, :])
            P = np.eye(self.N)
            for i in range(self.N):
                for j in range(self.N):
                    P[i,j] -= dz[i] * dz[j]/ np.dot(dz, dz)
                    # P[j,i]=P[i,j]

            for i,check in enumerate(znewm):
                if check<0:
                    znewm[i]=0
            z_new[m] = self.z - self.dt * P @ self.M[m] @ self.grad_F[m]
        self.z = z_new

    def reparameterize(self):
        """
        Smooth and reparameterize the string
        """
        zsmooth = np.zeros_like(self.z)
        zsmooth[0, :] = self.z[0, :]
        zsmooth[-1, :] = self.z[0, :]

        for m in range(1, self.R - 1):
            zself.smooth[m, :] = (1 - self.smooth) * self.z[m, :] + 
                0.5 * self.smooth * (self.z[m - 1] + self.z[m + 1])
        self.z = zsmooth
        
        # reparameterization
        L = np.zeros(self.R)
        for k in range(self.R):
            if k == 0:
                L[k] = np.linalg.norm(self.z[k, :] - self.z[k - 1, :]) + 1
            L[k] = np.linalg.norm(self.z[k, :] - self.z[k - 1, :]) + L[k - 1]

        segl=L[-1]/(self.R-1)
        k=1

        for m in range(1, self.R - 1):
            while m*segl>=L[k]:
                if k==self.R-1:
                    break
                k+=1
            zcurrent[m,:]=zsmooth[k-1,:]+(m*segl-L[k-1])*(zsmooth[k,:]-zsmooth[k-1,:])/np.linalg.norm(zsmooth[k,:]-zsmooth[k-1,:])

        np.save("zcurrent", zcurrent)
        np.save("gradF",gradF)
        np.save("../gradF",gradF)
        Dt=np.linalg.norm(zcurrent-z0)/(np.shape(z0)[0]*np.shape(z0)[1])
        Ddt=np.linalg.norm(zcurrent-zold)/(np.shape(z0)[0]*np.shape(z0)[1])
        f=open("../deviation.txt","a+")
        f.write(str(Dt)+" "+str(Ddt)+"\n")
        f.close()

    def update(self, cvs, derivs):
        self.compute_grad_F(cvs))
        self.compute_M(cvs, derivs)
        self.forward_euler()
        self.reparameterize()

    def write_to_file(self, path):
        """ Writes data about current grad_F, M, and parameterization to file"""
        np.save(f"{path}/gradF.npy", self.grad_F)
        np.save(f"{path}/M.npy", self.M)
        np.save(f"{path}/cvs.npy", self.z)


