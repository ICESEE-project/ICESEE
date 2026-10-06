"""non-cycled 2d test of LETKF solver with B or R localization"""
import numpy as np
from scipy.linalg import eigh, inv, pinvh
from scipy.special import gamma, kv
from argparse import ArgumentParser

# function definitions.

def cartdist(x1,y1,x2,y2,xmax,ymax):
    """cartesian distance on doubly periodic plane"""
    dx = np.abs(x1 - x2)
    dy = np.abs(y1 - y2)
    dx = np.where(dx > 0.5*xmax, xmax - dx, dx)
    dy = np.where(dy > 0.5*ymax, ymax - dy, dy)
    return np.sqrt(dx**2 + dy**2)

def gasp_cohn(r):
    """
    Gaspari-Cohn taper function.
    very close to exp(-(r/c)**2), where c = sqrt(0.15)
    r should be >0 and normalized so taper = 0 at r = 1
    """
    rr = 2.*r
    rr += 1.e-13 # avoid divide by zero warnings from numpy
    taper = np.where(r<=0.5, \
            ( ( ( -0.25*rr +0.5 )*rr +0.625 )*rr -5.0/3.0 )*rr**2 + 1.0,\
            np.zeros(r.shape,r.dtype))
    taper = np.where(np.logical_and(r>0.5,r<1.), \
            ( ( ( ( rr/12.0 -0.5 )*rr +0.625 )*rr +5.0/3.0 )*rr -5.0 )*rr \
               + 4.0 - 2.0 / (3.0 * rr), taper)
    return taper

def generalized_normal(r, beta):
    # https://en.wikipedia.org/wiki/Generalized_normal_distribution
    # beta=1 is exponential (laplace), beta=2 is gaussian.
    if beta < 1 or beta > 2:
        raise ValueError('1 <= beta <= 2 for generalized normal')
    return np.exp(-r**beta)

def rq(r, alpha):
    # rational quadratic cov function.
    # equivalent to a sum of gaussians with different length scales
    # length scale (l) parameter = 0.5
    return (1+r**2/alpha)**-alpha

def matern(r,v,l=1):
    # matern covariance function (v=0.5 is exponential, v->inf is gaussian)
    # overflow will result for values of v greater than about 35
    r[r == 0] = 1e-8
    part1 = 2 ** (1 - v) / gamma(v)
    part2 = (np.sqrt(2 * v) * r / l) ** v
    part3 = kv(v, np.sqrt(2 * v) * r / l)
    return part1 * part2 * part3

# CL args.
parser = ArgumentParser(description='test EnKF solvers (global solve with B loc vs local solve with R loc)')
parser.add_argument('--lscale', type=float, required=True, help='localization scale in grid points')
parser.add_argument('--covscale', type=float, required=False, default=0.1, help='covariance scale in grid points (as a fraction of domain size (ndim))')
parser.add_argument('--verbose', action='store_true', help='verbose output')
parser.add_argument('--localsolve', action='store_true', help='local analysis for B loc')
parser.add_argument('--cov_param', type=float, required=True, help='covariance parameter')
parser.add_argument('--nsamples', type=int, default=10, help='ensemble members')
parser.add_argument('--ntrials', type=int, default=100, help='number of trials')
parser.add_argument('--ndim', type=int, default=50, help='domain size (ndim x ndim square)')
parser.add_argument('--random_seed', type=int, default=0, help='random seed (default is to not set)')
parser.add_argument('--oberrvar', type=float, default=1.0, help='observation error variance')
args = parser.parse_args()
# update local namespace with CL args and values
locals().update(args.__dict__)
if verbose:
    print(args)

scale = ndim*covscale
if random_seed > 0: # set random seed for reproducibility
    np.random.seed(random_seed)

# specify covariance and localization matrices
# (two-dimensional periodic domain with grid ndim x ndim)
local = np.zeros((ndim**2,ndim**2),np.float64)
cov = np.zeros((ndim**2,ndim**2),np.float64) # B variance = 1
yg,xg = np.unravel_index(np.arange(ndim**2),(ndim,ndim))
nmid = ndim**2//2+ndim//2 # index of middle of domain
for n in range(ndim**2):
    dist = cartdist(xg,yg,xg[n],yg[n],ndim,ndim)
    #cov[:,n] = generalized_normal(dist/scale, cov_param)
    #cov[:,n] = rq(dist/scale, cov_param)
    cov[:,n] = matern(dist/scale, cov_param)
    local[:,n] = gasp_cohn(dist/lscale)
    n += 1

#cov2d = cov[nmid].reshape(ndim,ndim)
#import matplotlib.pyplot as plt
#plt.imshow(cov2d)
#plt.show()
#plt.figure()
#plt.plot(np.arange(ndim),cov2d[ndim//2])
#plt.show()
#raise SystemExit

# optimal gain matrix
Rm = oberrvar*np.eye(ndim**2)
kfopt = np.dot(cov, inv(cov + Rm))
paopt = np.dot((np.eye(ndim**2) - kfopt), cov)
if verbose:
    print('tr(pa)/tr(pb) = ',np.trace(paopt)/np.trace(cov))
kfopt_frobnorm = (kfopt**2).sum()

# compute eigenanalysis of true covariance matrix to sample.
evals, evecs = eigh(cov)

if verbose:
    for n in range(1,ndim**2):
        percentvar = evals[-n:].sum()/evals.sum()
        if percentvar > 0.99:
            nrank = n
            break
    print('rank of covariance matrix = %s' % nrank)
evals = evals.clip(min=np.finfo(evals.dtype).eps)
scaled_evecs = np.dot(evecs, np.diag(np.sqrt(evals)))/np.sqrt(nsamples-1)

# run trials with different ensembles
if verbose:
    kfensmean_bloc = np.zeros((ndim**2,ndim**2),np.float64)
    kfensmean_rloc = np.zeros((ndim**2,ndim**2),np.float64)
meankferr_bloc = 0.0
meankferr_rloc = 0.0
nlocal = 0; neig = 0
for ntrial in range(ntrials):

    # generate ensemble (x is an array of unit normal random numbers)
    x = np.random.normal(size=(ndim**2,nsamples))
    x = x - x.mean(axis=1)[:,np.newaxis] # zero mean
    # noise = np.random.multivariate_normal(np.zeros(ndim**2), cov)
    # x = x + noise[:,np.newaxis]
    y = np.dot(scaled_evecs, x)

    # compute kalman gain for global solution with B loc
    if not localsolve:
        cov_sample = local*np.dot(y,y.T)
        kfens_bloc = np.dot(cov_sample, inv(cov_sample + Rm))
    else:
        kfens_bloc = np.zeros((ndim**2,ndim**2),np.float64)

    # compute local solution for R loc
    kfens_rloc = np.zeros((ndim**2,ndim**2),np.float64)
    for n in range(ndim**2): # loop over analysis grid points

        # find local grid points (obs since H=I)
        dist = cartdist(xg,yg,xg[n],yg[n],ndim,ndim)
        indx = dist < np.abs(lscale)
        nmindist = np.argmin(dist[indx])
        ylocal = y[np.ix_(indx,np.ones(nsamples,np.bool_))].T

        # 'traditional' R localization
        if not localsolve:
            YbRinv = ylocal*local[indx,n]/oberrvar
            pa = np.eye(nsamples) + np.dot(YbRinv, ylocal.T)
            kfens_rloc[indx,n] = np.dot(y[n,:], np.dot(inv(pa), YbRinv))
        else:
            # R localization achieved by tapering ens perts a la
            # Sakov DOI 10.1007/s10596-010-9202-6)
            taper = np.sqrt(local[indx,n])
            ylocalloc = ylocal*taper # depends on distance between ob and analysis point
            YbRinv = ylocalloc/oberrvar
            pa = np.eye(nsamples) + np.dot(YbRinv, ylocalloc.T)
            # note the extra application of taper here
            kfens_rloc[indx,n] = taper*np.dot(y[n,:], np.dot(inv(pa), YbRinv))

            # B loc modulate ensemble with eigenvectors of 'local' localization matrix.
            if not ntrial:
                # compute and save modulation vectors.
                localloc = local[np.ix_(indx,indx)]
                if not nlocal: nlocal = localloc.shape[0]
                nlocal2 = localloc.shape[0]
                if nlocal != nlocal2:
                    raise ValueError('nlocal not constant')
                # symmetric square root of localization (truncated eigenvector expansion)
                evals, evecs = eigh(localloc)
                for nn in range(1,nlocal):
                    percentvar = evals[-nn:].sum()/evals.sum()
                    if percentvar > 0.99:
                        neigcount = nn
                        break
                if not neig: neig = neigcount
                if neigcount != neig:
                    raise ValueError('neig not constant')
                evecs_norm = (evecs*np.sqrt(evals/percentvar)).T
                if not n:
                    sqrtlocalloc = np.zeros((ndim**2,neig,nlocal),np.float64)
                sqrtlocalloc[n,...] = evecs_norm[nlocal-neig:nlocal,:]
            # modulated ensemble (permuted element-wise products of ylocal and sqrtlocalloc)
            #ylocal2 = np.multiply(np.tile(sqrtlocalloc[n],(nsamples,1)),np.tile(ylocal,(neig,1)))
            ylocal2 = np.multiply(np.repeat(sqrtlocalloc[n],nsamples,axis=0),np.tile(ylocal,(neig,1)))
            #ylocal2 = np.zeros((neig*nsamples,nlocal),ylocal.dtype); nsamp2 = 0
            #for j in range(neig):
            #    for nsamp in range(nsamples):
            #        ylocal2[nsamp2,:] = ylocal[nsamp,:]*sqrtlocalloc[n,neig-j-1,:]
            #        nsamp2 += 1
            YbRinv = ylocal2/oberrvar
            pa = np.eye(neig*nsamples) + np.dot(YbRinv, ylocal2.T)
            kfens_bloc[indx,n] = np.dot(ylocal2[:,nmindist], np.dot(inv(pa), YbRinv))
    
    # normalized Frobenius norm
    if verbose:
        kfensmean_rloc += kfens_rloc/ntrials
        kfensmean_bloc += kfens_bloc/ntrials
    diff_rloc = kfens_rloc-kfopt
    diff_bloc = kfens_bloc-kfopt
    kferr_rloc = (diff_rloc**2).sum()
    kferr_bloc = (diff_bloc**2).sum()
    meankferr_rloc += (kferr_rloc/kfopt_frobnorm)/ntrials
    meankferr_bloc += (kferr_bloc/kfopt_frobnorm)/ntrials

if verbose:
    import matplotlib.pyplot as plt
    plt.figure()
    nlocal = 2*int(lscale) - 1
    x = np.arange(ndim)-ndim//2
    x2 = np.linspace(ndim//2-nlocal//2,ndim//2+nlocal//2,nlocal)-ndim//2
    print(x2)
    kfrloc = kfensmean_rloc[nmid].reshape(ndim,ndim)
    kfbloc = kfensmean_bloc[nmid].reshape(ndim,ndim)
    kfopt2d = kfopt[nmid].reshape(ndim,ndim)
    plt.plot(x,kfrloc[ndim//2],'r',label='K rloc')
    plt.plot(x,kfbloc[ndim//2],'b',label='K bloc')
    #plt.plot(x,kfnoloc[ndim//2],'k:',label='K noloc')
    plt.plot(x,kfopt2d[ndim//2],'k',label='K true')
    plt.xlim(x2.min(), x2.max())
    plt.title('localization scale = %s beta = %s' % (lscale, cov_param))
    plt.legend()
    plt.savefig('gains.png')
    plt.show()

meankferr_rloc = np.sqrt(meankferr_rloc); meankferr_bloc = np.sqrt(meankferr_bloc)

if verbose:
    diff_rloc = kfensmean_rloc-kfopt
    diff_bloc = kfensmean_bloc-kfopt
    kferr_rloc = (diff_rloc**2).sum()
    kferr_bloc = (diff_bloc**2).sum()
    ensmeankferr_rloc = kferr_rloc/kfopt_frobnorm
    ensmeankferr_bloc = kferr_bloc/kfopt_frobnorm
    ensmeankferr_rloc = np.sqrt(ensmeankferr_rloc); ensmeankferr_bloc = np.sqrt(ensmeankferr_bloc)
    # print out mean error
    print("lscale = %s Kerr_Rloc = %6.4f Kerr_Bloc = %6.4f Kmerr_Rloc = %6.4f Kmerr_Bloc = %6.4f" % (lscale,meankferr_rloc,meankferr_bloc,ensmeankferr_rloc,ensmeankferr_bloc))
else:
    print("lscale = %s Kerr_Rloc = %6.4f Kerr_Bloc = %6.4f" % (lscale,meankferr_rloc,meankferr_bloc))