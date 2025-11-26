##############
## INITIATE ##
##############

## Imports
import torch
import numpy as np
import matplotlib.pyplot as plt
import time
import seaborn as sns
import yaml

#
###

#######################
## UTILITY FUNCTIONS ##
#######################

def hh(x, y):
  return (x * y)/(x + y)

def f_normal(z, mu, sigma):
    prefactor = 1 / (sigma * torch.sqrt(torch.tensor(2 * torch.pi, dtype=torch.float64)))
    exponent = -0.5 * ((z - mu) / sigma) ** 2
    density = prefactor * torch.exp(exponent)
    return density

def pdf(x):    
    x_norm = x/torch.clamp(x.sum(0), min=1e-8)
    return x_norm

#
###

########################
## SURVIVAL FUNCTIONS ##
########################

def p_s_i(z, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    p_s_i_ = p_sn_i(z, beta_sn0_i, beta_sn1_i, beta_sn2_i) * p_sp_i(z, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i)
    return p_s_i_

def p_sn_i(z, beta_sn0_i, beta_sn1_i, beta_sn2_i):
    p_sn_i_ = beta_sn0_i * torch.sigmoid((z - beta_sn1_i) * beta_sn2_i)    
    return p_sn_i_

def p_sp_i(z, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    p_sp_i_ = 1.0
    for j in range(len(NNt)):
      p_sp_i_ = p_sp_i_ * p_sp_ij(z, NNt[j], bbeta_dp0_i[j], bbeta_dp1_i[j], bbeta_dp2_i[j])
    return p_sp_i_

def p_sp_ij(z, Nt_j, beta_dp0_ij, beta_dp1_ij, beta_dp2_ij):
    p_sp_ij_ = torch.exp(Nt_j * torch.log(p_sp_ijk(z, beta_dp0_ij, beta_dp1_ij, beta_dp2_ij)))
    return p_sp_ij_

def p_sp_ijk(z, beta_dp0_ij, beta_dp1_ij, beta_dp2_ij):
    p_sp_ijk_ = 1 - beta_dp0_ij * torch.sigmoid((z - beta_dp1_ij) * beta_dp2_ij)
    return p_sp_ijk_

def total_p_dp_i(z, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    p_dp_i_ = 0.0
    for j in range(len(NNt)):
        p_dp_i_ = p_dp_i_ + (1 - p_sp_ij(z, NNt[j], bbeta_dp0_i[j], bbeta_dp1_i[j], bbeta_dp2_i[j]))
    return p_dp_i_

#
###

#########################################
## PARTITIONING STRUCTURED POPULATIONS ##
#########################################

def Zt_d_i(zz, NNt_i, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    Zt_d_i_ = zz * (1 - p_s_i(zz, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i)) * NNt_i
    return Zt_d_i_.sum()

def Zt_dn_i(zz, NNt_i, beta_sn0_i, beta_sn1_i, beta_sn2_i):
    Zt_dn_i_ = zz * (1 - p_sn_i(zz, beta_sn0_i, beta_sn1_i, beta_sn2_i)) * NNt_i
    return Zt_dn_i_.sum()

def ZZt_dp_i(zz, NNt_i, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    p_sp_i_ = p_sp_i(zz, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i)
    total_p_dp_i_ = total_p_dp_i(zz, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i)
    ZZt_dp_i_ = []

    for j in range(len(NNt)):
        if total_p_dp_i_.sum() < 0.001:
            ZZt_dp_i_.append(torch.tensor(0.0, dtype=torch.float64))
        else:
            weight = (
                (1 - p_sp_ij(zz, NNt[j], bbeta_dp0_i[j], bbeta_dp1_i[j], bbeta_dp2_i[j])).sum()
                / total_p_dp_i_.sum()
            )
            term = (
                zz
                * p_sn_i(zz, beta_sn0_i, beta_sn1_i, beta_sn2_i)
                * (1 - p_sp_i_)
                * weight
                * NNt_i
            ).sum()
            ZZt_dp_i_.append(term)

    return torch.stack(ZZt_dp_i_)  # <-- Return as tensor

#
###

#########################################
## PARTITIONING UNSTRUCTURED RESOURCES ##
#########################################

def p_s_0(beta_sn0_0, NNt, bbeta_dh0_0):
    p_sn_0_ = p_sn_0(beta_sn0_0)
    p_sh_0_ = p_sh_0(NNt, bbeta_dh0_0)
    p_s_0_ = p_sn_0_ * p_sh_0_
    return p_s_0_

def p_sn_0(beta_sn0_0):
    return beta_sn0_0

def p_sh_0(NNt, bbeta_dh0_0):
    p_sh_0_ = 1.0
    for j in range(len(NNt)):
      p_sh_0_ = p_sh_0_ * p_sh_0j(NNt[j], bbeta_dh0_0[j])
    return p_sh_0_

def p_sh_0j(Nt_j, beta_dh0_0j):
    p_sh_0j_ = torch.exp(Nt_j * torch.log(1 - beta_dh0_0j))
    return p_sh_0j_

def total_p_dh_0(NNt, bbeta_dh0_0):
    p_dh_0_ = 0.0
    for j in range(len(NNt)):
      p_dh_0_ = p_dh_0_ + (1 - p_sh_0j(NNt[j], bbeta_dh0_0[j]))
    return p_dh_0_

def Gt_d_0(Gt_0, beta_sn0_0, NNt, bbeta_dh0_0):
    Gt_d_0_ = (1 - p_s_0(beta_sn0_0, NNt, bbeta_dh0_0)) * Gt_0
    return Gt_d_0_

def Gt_dn_0(Gt_0, beta_sn0_0):
    Gt_dn_0_ = (1 - p_sn_0(beta_sn0_0)) * Gt_0
    return Gt_dn_0_

def GGt_dh_0(Gt_0, beta_sn0_0, NNt, bbeta_dh0_0):
    p_sh_0_ = p_sh_0(NNt, bbeta_dh0_0)
    total_p_dh_0_ = total_p_dh_0(NNt, bbeta_dh0_0)
    GGt_dh_0_ = []
    for j in range(len(NNt)):
      if total_p_dh_0_ == 0:
        GGt_dh_0_.append(torch.tensor(0.0, dtype=torch.float64))
      else:
        GGt_dh_0_.append((Gt_0 * p_sn_0(beta_sn0_0) * (1 - p_sh_0_) * (1 - p_sh_0j(NNt[j], bbeta_dh0_0[j])) / total_p_dh_0_))
    return torch.stack(GGt_dh_0_)

#
###

##################################
## CONSUMPTION AND ASSIMILATION ##
##################################

def g_i(z, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i):
    g_i_ = alpha_1_i * alpha_2_i * c_i(z, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i)
    return g_i_ * torch.ones_like(z)

def c_i(z, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i):
    epsilon = 1e-8
    c_i_ = gamma_1_i * gamma_2_i * Zt_c_i / (gamma_1_i * Nt_i + gamma_2_i * Zt_c_i + epsilon)
    return c_i_ * torch.ones_like(z)

#
###

##################
## REPRODUCTION ##
##################

def l_i(z, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i):
    l_i_ =  l_m_i(z, delta_i) + l_r_i(z, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i)
    return l_i_

def l_r_i(z, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i):
    l_r_i_ = mu_o_i * p_b_i(z, beta_b0_i, beta_b1_i, beta_b2_i) * n_o_i(z, nu_o_i, mu_o_i, beta_b1_i)
    return l_r_i_

def p_b_i(z, beta_b0_i, beta_b1_i, beta_b2_i):
    p_b_i_ = beta_b0_i * torch.sigmoid((z - beta_b1_i) * beta_b2_i)
    return p_b_i_

def n_o_i(z, nu_o_i, mu_o_i, beta_b1_i):
    n_o_i_ = nu_o_i * torch.ones_like(z) / mu_o_i * beta_b1_i
    # n_o_i_ = torch.tensor(1)
    return n_o_i_ 

def n_o_i_uncapped(z, nu_o_i, mu_o_i):
    n_o_i_ = nu_o_i * z / mu_o_i
    # n_o_i_ = torch.tensor(1)
    return n_o_i_ 

#
###

#################
## MAINTENANCE ##
#################

def l_m_i(z, delta_i):
    l_m_i_ = delta_i * z**0.75
    return l_m_i_

#
###

###############
## EXCRETION ##
###############

def Zt_e_i(zz, NNt_i, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i, delta_i, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    Zt_e_i_ = Zt_c_i - ((g_i(zz, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i) - l_m_i(zz, delta_i)) * p_s_i(zz, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i) * NNt_i).sum()    
    return Zt_e_i_

#
###

#############
## KERNELS ##
#############

def f_s_i(z, zp, sigma_s_i, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i):
    f_s_i_ = f_normal(z, mu_s_i(zp, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i), sigma_s_i)
    return f_s_i_

def mu_s_i(z, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i):
    mu_s_i_ = z + g_i(z, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i) - l_i(z, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i)
    return mu_s_i_

def f_r_i(z, zp, mu_o_i, sigma_o_i):
    f_r_i_ = f_normal(z, mu_o_i, sigma_o_i)
    return f_r_i_

#
###

###########################
## MAIN KERNEL FUNCTIONS ##
###########################

def K_i(zzz, zzzp, mu_o_i, sigma_o_i, sigma_s_i, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i, nu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i):
    K_r_i = pdf(f_r_i(zzz, zzzp, mu_o_i, sigma_o_i)) * n_o_i(zzzp, nu_o_i, mu_o_i, beta_b1_i) * p_b_i(zzzp, beta_b0_i, beta_b1_i, beta_b2_i)
    K_s_i = pdf(f_s_i(zzz, zzzp, sigma_s_i, alpha_1_i, alpha_2_i, gamma_1_i, Nt_i, gamma_2_i, Zt_c_i, nu_o_i, mu_o_i, beta_b0_i, beta_b1_i, beta_b2_i, delta_i))
    K_i_ = (K_r_i + K_s_i) * p_s_i(zzzp, beta_sn0_i, beta_sn1_i, beta_sn2_i, NNt, bbeta_dp0_i, bbeta_dp1_i, bbeta_dp2_i)
    return K_i_

#
###

#########################
## ECOSYSTEM FUNCTIONS ##
#########################

def DBPG_tdt(D_t, B_t, P_t, G_t, phi_BD, phi_GD, phi_DP, phi_DB, rho_B, phi_PG, rho_G, gamma_3_G, P_c, G_c, Z_e_t, Z_d_t):
  D_tdt = D_t + phi_BD * B_t + phi_GD * G_t - hh(phi_DB * D_t, rho_B * B_t) + Z_d_t + Z_e_t
  B_tdt = B_t + (1 - phi_DP) * hh(phi_DB * D_t, rho_B * B_t) - phi_BD * B_t
  P_tdt = P_t + phi_DP * hh(phi_DB * D_t, rho_B * B_t) - hh(phi_PG * P_t, rho_G * G_t) - P_c
  G_tdt = G_t + hh(phi_PG * P_t, rho_G * G_t) * gamma_3_G - phi_GD * G_t - G_c
  return D_tdt, B_tdt, P_tdt, G_tdt

#
###

#####################
## ECOSYSTEM MODEL ##
#####################

class YellowstoneEcosystemModel:

    ################
    ## PARAMETERS ##
    ################
    
    def __init__(self, param_file):
        
        ## Load YAML parameter file
        with open(param_file, 'r') as f:
            raw_params = yaml.safe_load(f)

        ## Convert all values to float64 torch tensors
        self.params = {}
        for key, value in raw_params.items():
            if isinstance(value, list):
                if isinstance(value[0], list):  ## 2D list
                    self.params[key] = torch.tensor(value, dtype=torch.float64)
                else:  ## 1D list
                    self.params[key] = torch.tensor(value, dtype=torch.float64)
            else:  ## scalar
                self.params[key] = torch.tensor(float(value), dtype=torch.float64)

        ## Assign specific parameters as attributes for direct access
        for key in self.params:
            setattr(self, key, self.params[key])

        ## TWEAKS
        ## Metabolism
        SCALING_METABOLISM = 0.125        
        for i in [0, 1, 2, 3, 4, 5]:        
            self.ddelta[i] *= SCALING_METABOLISM        
                        
    #
    ###

    ################
    ## SIMULATION ##
    ################
    
    def simulate(
        self,
        t_max=200,
        Dt=None,
        Bt=None,
        Pt=None,
        Gt=None,
        NNt=None,
        NNNzt=None,
        zzmax=None,
        zzz_len=100,
    ):
    
        ######################
        ## DEFAULTS + SETUP ##
                    
        # Default state values
        Dt = Dt if isinstance(Dt, torch.Tensor) else torch.tensor(1e9, dtype=torch.float64)
        Bt = Bt if isinstance(Bt, torch.Tensor) else torch.tensor(1e9, dtype=torch.float64)
        Pt = Pt if isinstance(Pt, torch.Tensor) else torch.tensor(1e9, dtype=torch.float64)
        Gt = Gt if isinstance(Gt, torch.Tensor) else torch.tensor(1e9, dtype=torch.float64)
        
        # Population sizes
        if NNt is None:
            NNt = torch.tensor([100000, 100000, 15000, 750, 30, 30], dtype=torch.float64)
        elif not isinstance(NNt, torch.Tensor):
            NNt = torch.tensor(NNt, dtype=torch.float64)
        
        # Max phenotype values
        if zzmax is None:
            zzmax = torch.tensor([2000, 2000, 500, 1000, 200, 200], dtype=torch.float64)
        elif not isinstance(zzmax, torch.Tensor):
            zzmax = torch.tensor(zzmax, dtype=torch.float64)
        
        # Placeholder for phenotype distributions
        if NNNzt is None:
            NNNzt = None  # Will be created below
        
        #####################
        ## INITIATE STATES ##
                
        # Time steps
        tt = range(0, t_max + 1)
        
        # Phenotype space — stacked tensor of shape (n_species, zzz_len)
        zzz = torch.stack([
            torch.linspace(0.01, zzmax[i].item(), zzz_len, dtype=torch.float64)
            for i in range(len(NNt))
        ])  # shape: (n_species, zzz_len)
        
        # Transposed phenotype grids for 2D use
        zzz_mesh = torch.stack([torch.meshgrid(zzz[i], zzz[i], indexing="ij")[0] for i in range(len(NNt))])
        zzz_mesh_T = zzz_mesh.transpose(1, 2)  # shape: (n_species, zzz_len, zzz_len)
        
        # Initial phenotype distributions
        if NNNzt is None:
            Zt0_mean = zzmax / 3.0           # shape: (n_species,)
            Zt0_std = Zt0_mean / 10.0        # shape: (n_species,)
        
            # Compute phenotype distribution (e.g., normal distribution scaled by NNt)
            NNNzt = torch.stack([
                NNt[i] * pdf(f_normal(zzz[i], Zt0_mean[i].item(), Zt0_std[i].item()))
                for i in range(len(NNt))
            ])  # shape: (n_species, zzz_len)
        else:
            NNNzt = torch.stack([z.clone().detach().to(dtype=torch.float64) for z in NNNzt])
        
        # Initial condition tensors (shape: [1] for scalars, or [n_species, zzz_len])
        DDt = Dt.unsqueeze(0)  # shape: (1,)
        BBt = Bt.unsqueeze(0)
        PPt = Pt.unsqueeze(0)
        GGt = Gt.unsqueeze(0)
        NNNt = NNt.unsqueeze(0)
        NNNNzt = NNNzt.unsqueeze(0)  # shape: (1, n_species, zzz_len)
        
        # Compute initial total biomass: sum across species and phenotypes
        Zt_tot = (zzz * NNNzt).sum(dim=1).sum()  # first sum: per species, second: total
        
        # Total biomass (tensors, shape: [1])
        ttotal_biomass = Dt + Bt + Pt + Gt + Zt_tot
        ttotal_biomass = ttotal_biomass.unsqueeze(0)  # for time tracking
        
        #####################
        ## MAIN SIMULATION ##

        ## Step
        for t in range(1, t_max+1):
            
            ############################
            ## ACQUISITION OF BIOMASS ##
                    
            # Available biomass through nutrients and herbivory
            Pt_d = Gt_d_0(Pt, 1 - self.phi_PG, NNt, self.bbeta_dh0_P)
            Pt_dn = Gt_dn_0(Pt, 1 - self.phi_PG)
            PPt_dh = GGt_dh_0(Pt, 1 - self.phi_PG, NNt, self.bbeta_dh0_P)  # shape: (n_species,)
            Gt_d = Gt_d_0(Gt, 1 - self.phi_GD, NNt, self.bbeta_dh0_G)
            Gt_dn = Gt_dn_0(Gt, 1 - self.phi_GD)
            GGt_dh = GGt_dh_0(Gt, 1 - self.phi_GD, NNt, self.bbeta_dh0_G)  # shape: (n_species,)
            Pt_dh = PPt_dh.sum()
            Gt_dh = GGt_dh.sum()
            
            # Scavenging: vectorized
            ZZt_dn_raw = torch.stack([
                Zt_dn_i(zzz[i], NNNzt[i], self.bbeta_sn0[i], self.bbeta_sn1[i], self.bbeta_sn2[i])
                for i in range(len(NNt))
            ])  # shape: (n_species,)
            ZZt_dn = torch.matmul(self.bbbeta_dn0.T, ZZt_dn_raw.unsqueeze(1)).flatten()  # shape: (n_species,)
            
            # Predation: vectorized
            ZZt_dp_raw = torch.stack([
                ZZt_dp_i(
                    zzz[i], NNNzt[i], self.bbeta_sn0[i], self.bbeta_sn1[i], self.bbeta_sn2[i],
                    NNt, self.bbbeta_dp0[i], self.bbbeta_dp1[i], self.bbbeta_dp2[i]
                )
                for i in range(len(NNt))
            ])  # shape: (n_species, zzz_len)
            ZZt_dp = ZZt_dp_raw.sum(dim=0)  # sum across predators
            
            # Diet preference tensor (convert from list if needed)
            if not isinstance(self.lllambda, torch.Tensor):
                self.lllambda = torch.tensor(self.lllambda, dtype=torch.float64)
            
            # Acquired biomass per species
            ZZt_c = (
                self.lllambda[:, 0] * PPt_dh +
                self.lllambda[:, 1] * GGt_dh +
                self.lllambda[:, 2] * ZZt_dn +
                self.lllambda[:, 3] * ZZt_dp
            )  # shape: (n_species,)
            
            # Not acquired biomass
            Zt_notc = (
                (1 - self.lllambda[:, 0]) * PPt_dh +
                (1 - self.lllambda[:, 1]) * GGt_dh +
                (1 - self.lllambda[:, 2]) * ZZt_dn +
                (1 - self.lllambda[:, 3]) * ZZt_dp
            ).sum()  # scalar
            
            # Excretion: still species-specific (keep as loop)
            Zt_e = torch.tensor(0.0, dtype=torch.float64)
            for i in range(len(NNt)):
                Zt_e += Zt_e_i(
                    zzz[i], NNNzt[i], self.aalpha_1[i], self.aalpha_2[i], self.ggamma_1[i],
                    NNt[i], self.ggamma_2[i], ZZt_c[i], self.ddelta[i], self.bbeta_sn0[i],
                    self.bbeta_sn1[i], self.bbeta_sn2[i], NNt,
                    self.bbbeta_dp0[i], self.bbbeta_dp1[i], self.bbbeta_dp2[i]
                )

            #####################
            ## POPULATION STEP ##
                    
            # Compute population at next time step
            updated_NNNzt = []
            
            for i in range(len(NNt)):
                if NNt[i] >= 1:
                    updated_dist = torch.matmul(
                        K_i(
                            zzz_mesh[i], zzz_mesh_T[i], self.mmu_o[i], self.ssigma_o[i],
                            self.ssigma_s[i], self.aalpha_1[i], self.aalpha_2[i],
                            self.ggamma_1[i], NNt[i], self.ggamma_2[i], ZZt_c[i],
                            self.nnu_o[i], self.bbeta_b0[i], self.bbeta_b1[i],
                            self.bbeta_b2[i], self.ddelta[i], self.bbeta_sn0[i],
                            self.bbeta_sn1[i], self.bbeta_sn2[i], NNt,
                            self.bbbeta_dp0[i], self.bbbeta_dp1[i], self.bbbeta_dp2[i]
                        ),
                        NNNzt[i]
                    )
                    updated_NNNzt.append(updated_dist)
                else:
                    updated_NNNzt.append(torch.zeros_like(NNNzt[i]))
            
            # Stack new population tensor
            NNNzt = torch.stack(updated_NNNzt)  # shape: (n_species, zzz_len)
            
            # Update total population per species
            NNt = NNNzt.sum(dim=1)  # shape: (n_species,)
            
            # Accumulate over time
            NNNNzt = torch.cat([NNNNzt, NNNzt.unsqueeze(0)], dim=0)  # (t, n_species, zzz_len)
            NNNt = torch.cat([NNNt, NNt.unsqueeze(0)], dim=0)        # (t, n_species)
            
            ####################
            ## ECOSYSTEM STEP ##
            
            # Next state of ecosystem compartments
            Dt, Bt, Pt, Gt = DBPG_tdt(
                Dt, Bt, Pt, Gt,
                self.phi_BD, self.phi_GD, self.phi_DP, self.phi_DB,
                self.rho_B, self.phi_PG, self.rho_G, self.gamma_3_G,
                P_c=Pt_dh, G_c=Gt_dh, Z_e_t=Zt_e, Z_d_t=Zt_notc
            )

            ## Check total biomass and correct D accordingly
            Zt_tot = ((zzz * NNNzt).sum(dim=1)).sum()  # shape: scalar
            biomass_t = Dt + Bt + Pt + Gt + Zt_tot
            Dt += ttotal_biomass[-1] - biomass_t
            
            # Accumulate ecosystem compartments
            DDt = torch.cat([DDt, Dt.unsqueeze(0)], dim=0)
            BBt = torch.cat([BBt, Bt.unsqueeze(0)], dim=0)
            PPt = torch.cat([PPt, Pt.unsqueeze(0)], dim=0)
            GGt = torch.cat([GGt, Gt.unsqueeze(0)], dim=0)
            
            # Total biomass: compute and accumulate
            Zt_tot = ((zzz * NNNzt).sum(dim=1)).sum()  # shape: scalar
            biomass_t = Dt + Bt + Pt + Gt + Zt_tot
            ttotal_biomass = torch.cat([ttotal_biomass, biomass_t.unsqueeze(0)], dim=0)
            
        return zzz, tt, DDt, BBt, PPt, GGt, NNNt, NNNNzt, ttotal_biomass
        
    #
    ###

#
###

#############################
## VISUALISATION FUNCTIONS ##
#############################

def plot_model(predictions):

    ## Format output
    zzz, tt, DDt, BBt, PPt, GGt, NNNt, NNNNzt, ttotal_biomass = predictions
    NNNNzt_new = [] # (time, species, phenotype) --> (species, time, phenotype)
    for i in range(6):
        NNNzt_new = []
        for NNNzt in NNNNzt:
            NNNzt_new += [NNNzt[i].clone()]
        NNNNzt_new += [torch.stack(NNNzt_new)]
    NNNNzt = torch.stack(NNNNzt_new)

    ## FIGURE 1
    plt.figure(figsize=(12, 4))

    ## Scales
    scales = [100000, 100000, 10, 10, 1, 1]
    labels = [f"Willow (x{scales[0]})",
              f"Aspen (x{scales[1]})",
              f"Elk (x{scales[2]})",
              f"Bison (x{scales[3]})",
              f"Wolves (x{scales[4]})",
              f"Cougars (x{scales[5]})",
             ]
    #
    plt.subplot(1,3,1)
    for i in range(6):
        plt.plot(NNNNzt[i,:].sum(1) / scales[i], label=labels[i])
    #
    # plt.title('$N^{(i)}(t)$')
    plt.xlabel('Time (years)')
    plt.ylabel('Number of individuals')
    #
    plt.legend(loc='upper right', frameon=False)    
    # plt.yscale('log')
    # plt.show()

    ## Define color vector
    color_vector = {
      'Total_biomass': 'black',
      'Dead organic matter': 'brown',
      'Decomposers': 'blue',
      'Nutrients': 'orange',
      'Grass': 'green'
    }

    ## Ecosystem dynamics
    plt.subplot(1,3,2)
    plt.plot(tt, ttotal_biomass, label='Total_biomass', color=color_vector['Total_biomass'])
    plt.plot(tt, DDt, label='Dead organic matter', color=color_vector['Dead organic matter'])
    plt.plot(tt, BBt, label='Decomposers', color=color_vector['Decomposers'])
    plt.plot(tt, PPt, label='Nutrients', color=color_vector['Nutrients'])
    plt.plot(tt, GGt, label='Grass', color=color_vector['Grass'])
    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('Time')
    plt.ylabel('Biomass (kg)')
    # plt.grid(c='white')

    ## Calculate proportions
    total_biomass = DDt.clone().detach() + BBt.clone().detach() + PPt.clone().detach() + GGt.clone().detach()
    D_prop = DDt.clone().detach() / total_biomass
    B_prop = BBt.clone().detach() / total_biomass
    P_prop = PPt.clone().detach() / total_biomass
    G_prop = GGt.clone().detach() / total_biomass

    ## Convert tensors to numpy for plotting
    proportions = torch.stack([D_prop, B_prop, P_prop, G_prop], dim=0).numpy()

    ## Stacked area plot
    plt.subplot(1,3,3)
    labels = ['Dead organic matter', 'Decomposers', 'Nutrients', 'Grass']
    colors = [color_vector[label] for label in labels]  # Extract colors for each label
    plt.stackplot(tt, proportions, labels=labels, colors=colors, alpha=0.7)

    ## Remove gap between plot and y-axis
    plt.xlim(0, np.max(tt)-1)
    plt.ylim(0,1)

    ## Customize plot
    plt.legend(loc='lower right', frameon=False)
    plt.xlabel('Time')
    plt.ylabel('Proportion of Total Biomass')
    # plt.title('Change in Biomass Proportions Over Time')
    # plt.grid(c='white')

    ## End figure
    plt.tight_layout()
    plt.show()


    ## FIGURE 2
    exponents = [.25, 0.25, .25, .25, .25, .25]
    for i in range(6):

        ## Visualise matrices
        plt.figure(figsize=(12, 4))
        #
        plt.subplot(1, 3, 1)
        plt.imshow(NNNNzt[i].T**exponents[i], origin='lower', aspect='auto')
        # plt.colorbar()
        plt.title('$N^{(i)}(z,t)$')
        plt.xlabel('Time (years)')
        plt.ylabel('Biomass Bins')
        #
        plt.subplot(1, 3, 2)
        for t, NNzt in enumerate(NNNNzt[i]):
            plt.plot(zzz[i], NNzt, alpha=0.2, color='black')
        plt.plot(zzz[i], NNzt, color='salmon')
        # plt.colorbar()
        plt.title('$N^{(i)}(z,t)$')
        plt.xlabel('Body mass (kg)')
        plt.ylabel('Density of individuals')
        #
        plt.subplot(1, 3, 3)
        plt.plot(NNNNzt[i,:].sum(1), color='salmon')
        # plt.colorbar()
        plt.title('$N^{(i)}(t)$')
        plt.xlabel('Time (years)')
        plt.ylabel('Number of individuals')
        #
        plt.tight_layout()
        plt.show()

#
###