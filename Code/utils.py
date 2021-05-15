import numpy as np
import scipy.stats as ss
from matplotlib import pyplot as plt
from math import gamma as GammaFN
from scipy.optimize import brentq as root
from scipy.optimize import fsolve

def findMoments(data):
    xbar = np.mean(data)
    std = np.std(data, ddof=1)
    skew = ss.skew(data,bias=False)
    
    return xbar, std, skew

def fitNormal(data, method):
    assert method == 'MLE' or method == 'MOM',"method must = 'MLE' or 'MOM'"
    
    xbar, std, skew = findMoments(data)
    
    if method == 'MLE':
        mu, sigma = ss.norm.fit(data)
        
    elif method == 'MOM':
        mu = xbar
        sigma = std
    
    return mu, sigma

def findNormalReturnPd(mu, sigma, T):
    
    q_T = ss.norm.ppf(1-1/T, mu, sigma)
    
    return q_T

def plotNormal(data, mu, sigma, title, figname):
    
    plt.hist(data, density=True)
    x = np.arange(np.min(data),np.max(data),(np.max(data)-np.min(data))/100.0)
    f_x = ss.norm.pdf(x, mu, sigma)
    plt.plot(x,f_x)
    plt.title(title)
    plt.xlabel('Flow (m^3/s)')
    plt.ylabel('Probability Density')
    plt.savefig(figname)
    plt.clf()
    
def NormalPPCT(data, mu, sigma, title, figname):
    x_sorted = np.sort(data)
    p_observed = np.arange(1,len(data)+1,1)/(len(data)+1)
    x_fitted = ss.norm.ppf(p_observed, mu, sigma)
    rho = np.corrcoef(x_sorted, x_fitted)[0,1]
    
    plt.scatter(x_sorted,x_fitted,color='b')
    plt.plot(x_sorted,x_sorted,color='r')
    plt.xlabel('Observations')
    plt.ylabel('Fitted Values')
    plt.title(title)
    plt.savefig(figname)
    plt.show()
    # plt.clf()
    
    # generate M synthetic samples of n observations
    rhoVector = np.zeros(10000)
    for i in range(10000):
        x = ss.norm.rvs(mu, sigma, len(data))
        rhoVector[i] = np.corrcoef(np.sort(x), x_fitted)[0,1]
        
    count = 0
    for i in range(len(rhoVector)):
        if rho < rhoVector[i]:
            count = count + 1
            
    p_value = 1 - count/10000
    
    return rho, p_value

# fit 2 or 3-parameter Gamma using MLE or MoM
def fitGamma(data, method, npars):
    assert method == 'MLE' or method == 'MOM',"method must = 'MLE' or 'MOM'"
    assert npars == 2 or npars == 3,"npars must = 2 or 3"
    
    xbar, std, skew = findMoments(data)

    if method == 'MLE':
        if npars == 2:
            shape, loc, scale = ss.gamma.fit(data,floc=0)
        elif npars == 3:
            shape, loc, scale = ss.gamma.fit(data)
            
        alpha = shape
        beta = 1/scale
        xi = loc
    elif method == 'MOM':
        if npars == 2:
            alpha = xbar**2/std**2
            beta = xbar/std**2
            xi = 0
        elif npars == 3:
            alpha = 4/skew**2
            beta = np.sqrt(alpha)/std
            xi = xbar - alpha/beta
    
    return alpha, xi, beta

def findGammaReturnPd(alpha, xi, beta, T):
    
    q_T = ss.gamma.ppf(1-1/T, alpha, xi, 1/beta)
    
    return q_T

def GammaPPCT(data, alpha, beta, xi):
    x_sorted = np.sort(data)
    p_observed = np.arange(1,len(data)+1,1)/(len(data)+1)
    x_fitted = ss.gamma.ppf(p_observed, alpha, xi, 1/beta)
    rho = np.corrcoef(x_sorted, x_fitted)[0,1]
    
    plt.scatter(x_sorted,x_fitted,color='b')
    plt.plot(x_sorted,x_sorted,color='r')
    
    # generate M synthetic samples of n observations
    rhoVector = np.zeros(10000)
    for i in range(10000):
        x = ss.gamma.rvs(alpha, xi, 1/beta, len(data))
        rhoVector[i] = np.corrcoef(np.sort(x), x_fitted)[0,1]
        
    count = 0
    for i in range(len(rhoVector)):
        if rho < rhoVector[i]:
            count = count + 1
            
    p_value = 1 - count/10000
    
    return rho, p_value

def fitGEV(data, method):
    assert method == 'MLE' or method == 'MOM',"method must = 'MLE' or 'MOM'"
    
    xbar, std, skew = findMoments(data)
    
    if method == 'MLE':
        kappa = root(lambda x: (np.abs(x)/x) * (-GammaFN(1+3*x) + 3*GammaFN(1+x)*GammaFN(1+2*x) - 2*GammaFN(1+x)**3) / (GammaFN(1+2*x)-GammaFN(1+x)**2)**(3/2) - skew, -1.0, 0.32)
        
        kappa, xi, alpha = ss.genextreme.fit(data, kappa)
    elif method == 'MOM':
        kappa = root(lambda x: (np.abs(x)/x) * (-GammaFN(1+3*x) + 3*GammaFN(1+x)*GammaFN(1+2*x) - 2*GammaFN(1+x)**3) / (GammaFN(1+2*x)-GammaFN(1+x)**2)**(3/2) - skew, -1.0, 0.32)
        alpha = np.sqrt(kappa**2 * std**2 / (GammaFN(1+2*kappa) - GammaFN(1+kappa)**2))
        print(xbar)
        print(alpha)
        print(kappa)
        xi = xbar - (alpha/kappa)*(1-GammaFN(1+kappa))
    
    return kappa, xi, alpha

def findGEVreturnPd(kappa, xi, alpha, T):
    
    q_T = ss.genextreme.ppf(1-1/T, kappa, xi, alpha)
    
    return q_T

def fitGPD(data, x0, method, initialize=False):
    assert method == 'MLE' or method == 'MOM',"method must = 'MLE' or 'MOM'"
    
    xbar, std, skew = findMoments(data)
    
    if method == 'MLE':
        if initialize == True:
            kappa = 0.5*(((xbar-x0)/std)**2-1)
            kappa, x0, alpha = ss.genpareto.fit(data, kappa, floc=x0)
        else:
            kappa, x0, alpha = ss.genpareto.fit(data, floc=x0)
    elif method == 'MOM':
        kappa = 0.5*(((xbar-x0)/std)**2-1)
        alpha = (1+kappa)*(xbar - x0)
    
    return kappa, x0, alpha

def fitWeibull(data, method, npars):
    assert method == 'MLE' or method == 'MOM',"method must = 'MLE' or 'MOM'"
    assert npars == 2 or npars == 3,"npars must = 2 or 3"
    
    xbar, std, skew = findMoments(data)

    if method == 'MLE':
        if npars == 2:
            kappa, xi, alpha = ss.weibull_min.fit(data,floc=0)
        elif npars == 3:
            kappa, xi, alpha = ss.weibull_min.fit(data)
    elif method == 'MOM':
        if npars == 2:
            kappa = root(lambda x: xbar**2 * (GammaFN(1+2/x)/GammaFN(1+1/x)**2 -1) - std**2, 0.02, 10)
            alpha = xbar / GammaFN(1+1/kappa)
            xi = 0
        elif npars == 3:
            def equations(p):
                kappa, xi, alpha = p
                mu = alpha*GammaFN(1+1/kappa) + xi
                sigma = np.sqrt(alpha**2*(GammaFN(1+2/kappa)-(GammaFN(1+1/kappa))**2))
                gamma = (GammaFN(1+3/kappa)*alpha**3 -3*mu*sigma**2 - mu**3)/(sigma**3)
                return (mu-xbar, sigma-std, gamma-skew)
            
            kappa = root(lambda x: xbar**2 * (GammaFN(1+2/x)/GammaFN(1+1/x)**2 -1) - std**2, 0.02, 10)
            alpha = xbar / GammaFN(1+1/kappa)
            xi = 0
            kappa, xi,alpha = fsolve(equations,(kappa,xi,alpha))
            
    return kappa, xi, alpha