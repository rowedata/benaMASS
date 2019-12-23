# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 10:08:48 2019
sampling

@author: Benazir
"""
import scipy
from scipy.stats import gamma
import numpy as np
import pandas as pd
import os 
from scipy import linalg
import random
from scipy.stats import uniform
from scipy.stats import bernoulli
from scipy.stats import norm


def set_prior_values(gt,ph,pcs):
    n = 980
    p = 5
    m = 10
    M = 3
    sigma_m_sq = 1
    sigma_a_sq = 1
    
    #STEP 1: Sample prior distributions
    G = np.array(gt.iloc[:, 4:985])
    
    #tau
    shape, scale = 2., 2.  # mean=4, std=2*sqrt(2)
    tau = np.random.gamma(shape / 2, scale / 2, 1)
    #mu 
    mu = pd.DataFrame(np.random.normal(loc=0.0, scale=sigma_m_sq / tau, size=n))
    h = np.random.uniform(low = 0.0, high = 1.0, size=1)
    log_pi = np.random.uniform(low=np.log(1 / p), high=np.log(M / p), size = 1)
    pi = np.exp(log_pi)
    
    gamma = np.random.binomial(1, pi, size = p)
    sum(gamma)
    
    #beta
    
    beta = np.random.normal(loc = 0.0, scale = sigma_a_sq / tau, size = sum(gamma))

    gamma_df = pd.DataFrame(data = gamma)
    gamma_df.loc[gamma==1] = beta
    beta=gamma_df
    
    # alpha
    X = np.array(pcs.iloc[0:980, 2:12])
    g_prior_var = n * linalg.inv(X.T.dot(X)) / tau
    def is_pos_def(x):
        return np.all(np.linalg.eigvals(x) > 0)
    is_pos_def(g_prior_var)
    
    mean=np.zeros(10)
    cov =g_prior_var
    alpha = np.random.multivariate_normal(mean, cov)
    
    #priors are set, now
    return tau, mu, h, pi, gamma,beta, alpha
#STEP 2: Proposal and sampling
#propose gamma

#remove
def remove(gamma):
    a = gamma
    remove_index = random.randint(1,sum(a))
    q = list(np.cumsum(a))
    a[q.index(remove_index)] = 0
    print("We add/remove ",q.index(remove_index)+1)
    
    return a
#add
def add(gamma):
    b = np.logical_not(gamma)
    b = remove(b)
    a = np.logical_not(b)
    return a
#swap
def swap(gamma):
    gamma = add(gamma)
    gamma = remove(gamma)
    
    return gamma


def propose(gamma, h, p):
    a = gamma
    gamma_prop = swap(a)
# propose pi
    a = sum(gamma_prop)
    pi_prop = np.random.beta(a, p - a + 1, size = None)
    print("proposed pi" , pi_prop)
# propose h
    h_prop = h  + np.random.uniform(-0.1, 0.1)
    
    return h_prop, pi_prop, gamma_prop

def calculate_r(h_prop, pi_prop, gamma_prop, h, pi, gamma): 

    #h
    uniform.ppf(h)
    uniform.ppf(h_prop)
    
    #pi
    uniform.ppf(pi)
    uniform.ppf(pi_prop)
    
    #gamma
    g = pi**gamma.sum()*(1-pi)**(p-gamma.sum()) 
    g_prop = pi_prop**gamma.sum()*(1-pi_prop)**(p-gamma.sum())
    
    #y current
    
    G = genotype.iloc[:, 4:984]
    
    Q = G.T.dot(beta).reset_index()
    Q = Q.drop(columns = 'index')
    
    mean = np.array(mu + Q).ravel()
    cov = np.identity(980) * (1 / tau)
    ys = np.random.multivariate_normal(mean, cov)
    mylist = norm.pdf(ys)
    y = np.prod(np.array(mylist))

    #proposed
    gamma_df = pd.DataFrame(data = gamma_prop)
    gamma_df.loc[gamma==1] = beta
    beta_prop=gamma_df
    
    
    Q_prop = G.T.dot(beta_prop).reset_index()
    Q_prop = Q_prop.drop(columns = 'index')
    
    mean_prop = np.array(mu + Q_prop).ravel()
    ys= np.random.multivariate_normal(mean_prop, cov)
    mylist = norm.pdf(ys)
    y_prop = np.prod(np.array(mylist))
    
    r = (uniform.ppf(h_prop)*uniform.ppf(pi_prop)*g_prop* )/()
    
    if (u < r): 
        h = h_prop
        pi = p_prop
        gamma = gamma_prop
    
    return h, pi, gamma

    #sample tau
    #select relevant covariates
def sample_tau_beta():    
    G = genotype.iloc[:, 4:984]
    G = G.iloc[gamma,]
    G_gamma = G.loc[G.index==1].T.reset_index(drop=True)
    
    G_square = np.square(G)
    s_j = np.sum(G_square, axis =1)
    sigma_a = h/(1 - h) * (1 / s_j.sum())
    
    omega = ((1/sigma_a**2) * np.identity(G_gamma.shape[1]) + G_gamma.T.dot(G_gamma))**(-1)
    
    Xty = G_gamma.T.dot(phenotype)
    Ytx = phenotype.T.dot(G_gamma)
    tau_scale = phenotype.T.dot(phenotype)-Ytx.dot(omega).dot(Xty)
    
    tau = np.random.gamma(n/2, scale = 0.5*tau_scale, size = None)
    
    #sample beta
    mean_beta = omega.dot(Xty)
    cov_beta = np.multiply((1/tau), omega)
    beta = np.random.multivariate_normal(mean_beta, cov_beta)
    
    return tau, beta

    



if __name__ == '__main__':
   
    #global variables
    n = 980 #subjects
    p = 5 #snps
    m = 10 #pc covariates
    M = 3
    
    #READ DATA 840 people by 5 rsIDs
    genotype = pd.read_csv('data/test.mgt.txt', sep = " ", header = None)
    phenotype = pd.read_csv('data/test.ph.txt', sep = " ", header = None)
    principal_components = pd.read_csv('data/qcvcf.eigenvec', sep = " ", header = None)
    
    #set priors
    tau, mu, h, pi, gamma, beta, alpha = set_prior_values(genotype,phenotype,principal_components)
    
    iterations = 1000
    #ITERATION START
    #propose h,p,gamma for Metropolis
    for (i in 1: iterations):
        print i
    h_prop, pi_prop, gamma_prop = propose(gamma, h, p)
    
    h, pi, gamma = calculate_r(h_prop, pi_prop, gamma_prop, h, pi, gamma) #still not done
    
    tau,beta = sample_tau_beta( )
    
    #update alpha
    beta = beta.iloc[gamma,]
    beta_gamma = beta.loc[beta.index==1].reset_index(drop=True)
    G_gamma.reset_index(drop=True)
    
    Q = phenotype - G_gamma.dot(beta_gamma) #solve this!
    pref = n / (n + 1)
    inv = linalg.inv(X.T.dot(X))
    Xtq = X.T.dot(Q)
    
    mean_alpha = pref *inv.dot(Xtq)
    cov_alpha = np.multiply((pref / tau) , )
    alpha = np.random.multivariate_normal(mean_alpha, cov_alpha)
    
    #ITERATION END
    
    PIP= X_included / 1000 #may have high sampling variance 

    
  




