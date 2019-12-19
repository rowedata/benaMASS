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


def propose(gamma,h):
    a = gamma
    gamma_prop = swap(a)
# propose pi
    a = sum(gamma_prop)
    pi_prop = np.random.beta(a, p - a + 1, size = None)
    print("proposed pi" , pi_prop)
# propose h
    h_prop = h  + np.random.uniform(-0.1,0.1)
    
    return h_prop, pi_prop, gamma_prop

def calculate_r(): 

    #h
    uniform.ppf(h)
    uniform.ppf(h_prop)
    
    #pi
    uniform.ppf(pi)
    uniform.ppf(pi_prop)
    
    #gamma vect
    gamma_df = pd.DataFrame(data = gamma)
    gamma_df.loc[gamma==0, 'new'] = pi
    gamma_df.loc[gamma==1, 'new'] = 1 - pi
    
    gamma_prop_df =pd.DataFrame(data = gamma_prop)
    gamma_prop_df.loc[gamma==0, 'new'] = pi_prop
    gamma_prop_df.loc[gamma==1, 'new'] = 1 - pi_prop
    
    gamma_prop_df.new
    gamma_df.new
    
    #y
    
    G = genotype.iloc[:,4:984]
    Q = G.T.dot(beta)
    mean = (mu + Q).T
    cov = np.identity(980) * (1 / tau)
    y = np.random.multivariate_normal(mean, cov)
    norm.pdf(y)



if __name__ == '__main__':
    #READ DATA 840 people by 5 rsIDs

    genotype = pd.read_csv('test.mgt.txt', sep = " ", header = None)
    phenotype = pd.read_csv('test.ph.txt', sep = " ", header = None)
    principal_components = pd.read_csv('qcvcf.eigenvec', sep = " ", header = None)
    #set priors
    tau, mu, h, pi, gamma, beta, alpha = set_prior_values(genotype,phenotype,principal_components)
    #propose h,p,gamma for Metropolis
    
    h_prop, pi_prop, gamma_prop = propose(gamma, h)
    r =calculate_r()
    
    u=np.random.uniform(0,1)
    
    if (u<r): 
        h=h_prop
        pi=p_prop
        gamma=gamma_prop
        
        
    #sample tau
    
    #sample beta
    
    
    
    



#add







