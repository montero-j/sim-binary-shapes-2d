# Correlacion entre avalanchas consecutivas (R_a = 0.1)

86 corridas (Cuad-Disco: 26, Pent-Disco: 26, Tri-Disco: 22, Disco-Disco: 12), N entre 37 y 281 avalanchas.

- Spearman (s_i, s_i+1) con p < 0.05: 5 de 86 corridas (por azar se esperan ~4.3).
    - Tri-Disco chi = 0.02: r_S = +0.144, p = 0.0344, p(BH) = 0.592, SEM(8)/SEM(1) = 1.07 +- 0.15
    - Tri-Disco chi = 0.15: r_S = +0.220, p = 0.00321, p(BH) = 0.276, SEM(8)/SEM(1) = 1.18 +- 0.18
    - Cuad-Disco chi = 0.01: r_S = -0.139, p = 0.0309, p(BH) = 0.592, SEM(8)/SEM(1) = 0.88 +- 0.12
    - Cuad-Disco chi = 0.02: r_S = -0.150, p = 0.0219, p(BH) = 0.592, SEM(8)/SEM(1) = 0.91 +- 0.12
    - Pent-Disco chi = 0.06: r_S = +0.159, p = 0.0228, p(BH) = 0.592, SEM(8)/SEM(1) = 1.18 +- 0.17
- Despues de corregir por comparaciones multiples (Benjamini-Hochberg sobre las 86 corridas): 0 con p(BH) < 0.05.
- |rho(1)| > 1.96/sqrt(N): 4 de 86: Cuad-Disco chi = 0.01 (rho1 = -0.164, banda 0.126), Cuad-Disco chi = 0.35 (rho1 = -0.168, banda 0.154), Pent-Disco chi = 0.06 (rho1 = +0.178, banda 0.136), Pent-Disco chi = 0.07 (rho1 = +0.194, banda 0.137)
- Lags k = 1..20 fuera de la banda: 83 de 1698 (por azar se esperan ~85, 5 %).
- SEM(8)/SEM(1): mediana 0.99, rango 0.62-1.30; ruido esperado para avalanchas independientes ~ +-0.15 (1 sigma). Corridas con cociente > 1 + 2 sigma: ninguna.
- Disco-Disco chi = 0.12 (borrador: r_S = -0.11): r_S = -0.109, p = 0.141 (p(BH) = 0.722), banda 1.96/sqrt(N) = 0.144 con N = 186; rho(1) = -0.072; SEM(8)/SEM(1) = 0.95 +- 0.14.
