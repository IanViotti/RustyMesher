import sympy as sp

# Constantes da iteração n
alpha, omega = sp.symbols('alpha omega')
A_ij, B_ij, C_ij = sp.symbols('A_ij B_ij C_ij')
L_x = sp.Symbol('L_x')

# Variáveis Delta intermediárias (Sweep xi) e finais (Sweep eta)
dx_star_im1 = sp.Symbol('dx*_i-1,j')
dx_star_ij  = sp.Symbol('dx*_i,j')
dx_star_ip1 = sp.Symbol('dx*_i+1,j')

dx_jm1 = sp.Symbol('dx_i,j-1')
dx_ij  = sp.Symbol('dx_i,j')
dx_jp1 = sp.Symbol('dx_i,j+1')

# ==========================================
# Passo 1: Sweep na direção Xi (Gera dx_star)
# Eq: (alpha - A_ij * delta_xixi) * dx_star = alpha * omega * L_x
# ==========================================

LHS_xi = alpha * dx_star_ij - A_ij * (dx_star_ip1 - 2*dx_star_ij + dx_star_im1)
RHS_xi = alpha * omega * L_x

# O "Pulo do Gato": Expandir a equação antes de pedir os coeficientes
LHS_xi_exp = sp.expand(LHS_xi)

print("=== COEFICIENTES DE THOMAS (Sweep Xi) ===")
# Usando .coeff(), o SymPy vai cirurgicamente no termo exato
print("a[i] (i-1) =", LHS_xi_exp.coeff(dx_star_im1))
print("b[i] (i)   =", LHS_xi_exp.coeff(dx_star_ij))
print("c[i] (i+1) =", LHS_xi_exp.coeff(dx_star_ip1))
print("d[i] (Dir) =", RHS_xi)


# ==========================================
# Passo 2: Sweep na direção Eta (Gera o dx final)
# Eq: (alpha - C_ij * delta_etaeta) * dx = alpha * dx_star
# ==========================================

LHS_eta = alpha * dx_ij - C_ij * (dx_jp1 - 2*dx_ij + dx_jm1)
RHS_eta = alpha * dx_star_ij

LHS_eta_exp = sp.expand(LHS_eta)

print("\n=== COEFICIENTES DE THOMAS (Sweep Eta) ===")
print("a[j] (j-1) =", LHS_eta_exp.coeff(dx_jm1))
print("b[j] (j)   =", LHS_eta_exp.coeff(dx_ij))
print("c[j] (j+1) =", LHS_eta_exp.coeff(dx_jp1))
print("d[j] (Dir) =", RHS_eta)