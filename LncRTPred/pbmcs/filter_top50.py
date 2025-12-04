#!/usr/bin/env python
import pandas as pd

# === Cargar ===
df = pd.read_csv("Resultado_Final.csv")

# === Calcular promedio y nº de modelos >=80% ===
cols = ["Decision_Tree(%)","K-Nearest-Neighbours(%)","Random_Forest(%)","LightGBM(%)"]
df["avg_score"] = df[cols].mean(axis=1)
df["models_ge80"] = (df[cols] >= 80).sum(axis=1)

# === Filtrado ===
# Mantener pares con al menos 3 de 4 modelos >=80% y promedio >=85
filtered = df[(df["avg_score"] >= 85) & (df["models_ge80"] >= 3)]

# Orden descendente por promedio
filtered = filtered.sort_values("avg_score", ascending=False)

# Top 50
top50 = filtered.head(50)

# === Guardar ===
out = "Resultado_Final.TOP50_consensus.csv"
top50.to_csv(out, index=False)
print(f"[OK] Guardado Top 50 en {out} ({len(top50)} filas).")
