####################
## DEFINE WRAPPER ##
####################

import torch

## Probabilities to rates
def r2p(x):
    return 1 - 1/(x+1)
def p2r(x):
    return x/(1-x)

## Define wrapper
def wrapper_YEM_V1_14(ttheta, t_max, model_r, model_p, Dt0, Bt0, Pt0, Gt0, NNt0, NNNzt0, zzmax, zzz_len):
    """
    Wrapper function that takes a scaling parameter vector and a reference model and applies scaling to each parameter of 
    the reference model to generate a candidate model.
    """
    ## Initiate idx
    idx = 0
    
    ## Initial conditions unstructured states
    Dt0_p = Dt0 * ttheta[idx] * 1; idx += 1
    Bt0_p = Bt0 * ttheta[idx] * 1; idx += 1
    Pt0_p = Pt0 * ttheta[idx] * 1; idx += 1
    Gt0_p = Gt0 * ttheta[idx] * 1; idx += 1
    ## Parameter index: 4
    
    ## Initial conditions structured states
    NNt0_p = NNt0.clone()
    NNNzt0_p = NNNzt0.clone()        
    for i in range(len(NNt0_p)):
        NNt0_p[i] = NNt0[i] * ttheta[idx]
        NNNzt0_p[i] = NNNzt0[i] * ttheta[idx]
        idx += 1
    ## Parameter index: 4 + 6

    ## Unstructured variables 
    model_p.phi_BD = r2p(p2r(model_r.phi_BD) * ttheta[idx]); idx += 1 # [0,1] 
    model_p.phi_DB = r2p(p2r(model_r.phi_DB) * ttheta[idx]); idx += 1 # [0,1] 
    model_p.rho_B = model_r.rho_B * ttheta[idx]; idx += 1 # R+ 
    model_p.phi_DP = r2p(p2r(model_r.phi_DP) * ttheta[idx]); idx += 1 # [0,1]         
    model_p.phi_GD = r2p(p2r(model_r.phi_GD) * ttheta[idx]); idx += 1 # [0,1]
    model_p.phi_PG = r2p(p2r(model_r.phi_PG) * ttheta[idx]); idx += 1 # [0,1]
    model_p.rho_G = model_r.rho_G * ttheta[idx]; idx += 1 # R+
    model_p.gamma_3_G = model_r.gamma_3_G * ttheta[idx]; idx += 1 # R+    

    ## Diet preference (I * 4) 
    model_p.lllambda[4][2] = r2p(p2r(model_r.lllambda[4][2]) * ttheta[idx]); idx += 1 # [0,1]
    model_p.lllambda[5][2] = r2p(p2r(model_r.lllambda[5][2]) * ttheta[idx]); idx += 1 # [0,1]
    # Only for scavenging in wolves and cougars
    
    ## Predation
    for i, j in [[0,2], [1,2], [2,4], [2,5], [3,4], [4,4], [5,5]]:        
        model_p.bbbeta_dp0[i][j] = r2p(p2r(model_r.bbbeta_dp0[i][j]) * ttheta[idx]); idx += 1 # [0,1]
        model_p.bbbeta_dp1[i][j] = model_r.bbbeta_dp1[i][j] * ttheta[idx]; idx += 1 # R
        model_p.bbbeta_dp2[i][j] = model_r.bbbeta_dp2[i][j] * ttheta[idx]; idx += 1 # R               
        
    for i in range(len(NNt0_p)):

        ## Nutrient assimilation and grazing (I * 2)
        model_p.bbeta_dh0_P[i] = r2p(p2r(model_r.bbeta_dh0_P[i]) * ttheta[idx]); idx += 1 # [0,1]
        model_p.bbeta_dh0_G[i] = r2p(p2r(model_r.bbeta_dh0_G[i]) * ttheta[idx]); idx += 1 # [0,1]

        ## Natural survival (I * 3)
        model_p.bbeta_sn0[i] = r2p(p2r(model_r.bbeta_sn0[i]) * ttheta[idx]); idx += 1 # [0,1]
        model_p.bbeta_sn1[i] = model_r.bbeta_sn1[i] * ttheta[idx]; idx += 1 # R
        model_p.bbeta_sn2[i] = model_r.bbeta_sn2[i] * ttheta[idx]; idx += 1 # R    

        ## Consumption and growth (I * 4)
        model_p.aalpha_1[i] = r2p(p2r(model_r.aalpha_1[i] * model_r.aalpha_2[i]) * ttheta[idx]); idx += 1 # [0,1]
        model_p.aalpha_2[i] = torch.ones_like(model_p.aalpha_2[i]) # Un-identifiable
        model_p.ggamma_1[i] = model_r.ggamma_1[i] * ttheta[idx]; idx += 1 # R+
        model_p.ggamma_2[i] = r2p(p2r(model_r.ggamma_2[i]) * ttheta[idx]); idx += 1 # [0,1]        

        ## Reproduction (I * 5)
        model_p.bbeta_b0[i] = r2p(p2r(model_r.bbeta_b0[i]) * ttheta[idx]); idx += 1 # [0,1]
        model_p.bbeta_b1[i] = model_r.bbeta_b1[i] * ttheta[idx]; idx += 1
        model_p.bbeta_b2[i] = model_r.bbeta_b2[i] * ttheta[idx]; idx += 1        
        model_p.nnu_o[i] = model_r.nnu_o[i] * ttheta[idx]; idx += 1
        model_p.mmu_o[i] = model_r.mmu_o[i] * ttheta[idx]; idx += 1

        ## Metabolic costs (I * 1)
        model_p.ddelta[i] = model_r.ddelta[i] * ttheta[idx]; idx += 1        

    ## DEBUG (check number of parameters)
    # print(idx)
    
    ## Simulate model
    predictions = model_p.simulate(t_max=t_max, Dt=Dt0_p, Bt=Bt0_p, Pt=Pt0_p, Gt=Gt0_p, NNt=NNt0_p, NNNzt=NNNzt0_p, zzmax=zzmax, zzz_len=zzz_len)

    return predictions

def build_ttheta_index_map(I, J):
    """
    Build a dictionary mapping parameter names to their index in ttheta.
    I = number of species
    J = number of interaction targets per species
    """
    index_map = {}
    idx = 0

    # Initial unstructured state
    for name in ['Dt0', 'Bt0', 'Pt0', 'Gt0']:
        index_map[name] = idx
        idx += 1

    # Initial structured state (NNt0 and NNNzt0)
    for i in range(I):
        index_map[f'NNt0[{i}]'] = idx
        # index_map[f'NNNzt0[{i}]'] = idx # Not used as theta maps to two parameters
        idx += 1

    # Unstructured parameters (9)
    for name in ['phi_BD', 'phi_DB', 'rho_B', 'phi_DP', 'phi_GD', 'phi_PG', 'rho_G', 'gamma_3_G']:
        index_map[name] = idx
        idx += 1        

    ## Diet preference
    index_map[f'lllambda[{4}][{2}]'] = idx; idx +=1
    index_map[f'lllambda[{5}][{2}]'] = idx; idx +=1

    ## Predation
    for i, j in [[0,2], [1,2], [2,4], [2,5], [3,4], [4,4], [5,5]]:        
        index_map[f'bbbeta_dp0[{i}][{j}]'] = idx; idx += 1
        index_map[f'bbbeta_dp1[{i}][{j}]'] = idx; idx += 1
        index_map[f'bbbeta_dp2[{i}][{j}]'] = idx; idx += 1
    
    # Nutrient assimilation & grazing (2 per species)
    for i in range(I):
        
        ## Nutrients and grass acquisition
        index_map[f'bbeta_dh0_P[{i}]'] = idx; idx += 1
        index_map[f'bbeta_dh0_G[{i}]'] = idx; idx += 1

        # Natural survival (3 per species)
        index_map[f'bbeta_sn0[{i}]'] = idx; idx += 1
        index_map[f'bbeta_sn1[{i}]'] = idx; idx += 1
        index_map[f'bbeta_sn2[{i}]'] = idx; idx += 1

        # Consumption and growth (4 per species)
        index_map[f'aalpha_1[{i}]'] = idx; idx += 1
        # index_map[f'aalpha_2[{i}]'] = 'unidentifiable'  # Not used as theta maps to two parameters
        index_map[f'ggamma_1[{i}]'] = idx; idx += 1
        index_map[f'ggamma_2[{i}]'] = idx; idx += 1

        # Reproduction (5 per species)
        index_map[f'bbeta_b0[{i}]'] = idx; idx += 1
        index_map[f'bbeta_b1[{i}]'] = idx; idx += 1
        index_map[f'bbeta_b2[{i}]'] = idx; idx += 1
        index_map[f'nnu_o[{i}]'] = idx; idx += 1
        index_map[f'mmu_o[{i}]'] = idx; idx += 1

        # Metabolic cost (1 per species)
        index_map[f'ddelta[{i}]'] = idx; idx += 1
    
    ## DEBUG (check number of parameters)
    # print(idx)

    return index_map

def get_index_from_name(name, index_map):
    """
    Given a parameter name, return the index in ttheta using the provided map.
    """
    if name not in index_map:
        raise ValueError(f"Parameter name '{name}' not found in index map.")
    return index_map[name]

def get_name_from_index(index, index_map):
    """
    Given an index, return the parameter name using the provided map.
    """
    for name, idx in index_map.items():
        if idx == index:
            return name
    raise ValueError(f"Index '{index}' not found in index map.")
     
#
###