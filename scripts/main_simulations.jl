using DrWatson
@quickactivate "moth_infestation_jl"

using LinearAlgebra
using Distributions
using DataFrames
using DelimitedFiles
using CSV
using Statistics
using Plots
using Random

BLAS.set_num_threads(1)

include(srcdir("main_functions.jl"))
include(srcdir("aux_functions.jl"))

#Reading distance matrices

D_paths=readdir(srcdir("distance_matrices"); join=true)
D_matrices=readdlm.(D_paths, '\t')
sites=string.(readdlm(srcdir("sites_names.txt"), '\t'))
sites=dropdims(sites, dims=2)

#Creating matrices of probability of site j to infect site i

P_matrices=prob_matrix.(D_matrices, 1000, 1.0)
P_logs=P_log.(P_matrices)
T_matrices=T_prob.(P_matrices)

#Running simulations
#The function infest_sim uses as input the following arguments:

#A matrix of direct probability of infections between sites (P)
#A matrix with the natural logarithm of the probability of infections (facilitate computations)
#A matrix with the total, direct and indirect probability of infections (computations of the direct and indirect probability of infections between sites)
#The number of time steps to run the model
#Parameter α of the model
#Number of simulations to perform

#The function return a dataframe containing the following variables:

#area = the total number of sites_names
#alpha = value of the parameter α used in the simulation
#site = ID of each site
#inf_time = Time at which the site was infected
#avg_inftime = Average infestation time of all sites
#prop_inf = Proportion of infested sites at the end of the simulation
#infest_focus = ID of the site that was the initial focus of infestation
#T_in = Total (direct and indirect) probability for a site to be infected by the initial focus of infestation
#T_out = Total (direct and indirect) probability of the inicial focus of infestation to infest all other sites
#P_in = Direct probability for a site to be infected by the initial focus of infestation
#P_out = Direct probability of the inicial focus of infestation to infest all other sites
#nsim = Number of the simulation

infest_sim(P_matrices, P_logs, T_matrices, sites, 30000, 1.0, 2000)