# v8_split_11/kinetics/models/yeast_acetate_contois_o2_v1.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple
import numpy as np
from ..core import ParameterSpec, Environment

@dataclass(frozen=True)
class YeastAcetateContoisO2V1:
    """Contois growth on acetate with oxygen limitation. DO is reactor-owned (env.DO_mol_L)."""
    name: str = "yeast_acetate_contois_o2_v1"
    state_ids: Tuple[str, ...] = ("X_g_L", "SAc_g_L")
    param_schema: Tuple[ParameterSpec, ...] = (
        ParameterSpec("mu_max_h", default=0.30, lower=0.0, units="1/h"),
        ParameterSpec("Kc_gS_gX", default=0.50, lower=1e-12, units="gS/gX"),
        ParameterSpec("Yxs_gX_gS", default=0.50, lower=1e-12, units="gX/gS"),
        ParameterSpec("K_DO_mol_L", default=1e-4, lower=0.0, units="mol/L"),
        ParameterSpec("YxO2_gX_per_molO2", default=10.0, lower=1e-12, units="gX/molO2"),
        ParameterSpec("S_min_g_L", default=0.0, lower=0.0, units="g/L"),
        ParameterSpec("X_min_g_L", default=1e-12, lower=0.0, units="g/L"),
    )

    def rhs(self, t_h: float, y: np.ndarray, env: Environment, params: Mapping[str, float]):
        X = float(y[0]); S = float(y[1])
        mu_max=float(params["mu_max_h"]); Kc=float(params["Kc_gS_gX"]); Yxs=float(params["Yxs_gX_gS"])
        K_DO=float(params["K_DO_mol_L"]); YxO2=float(params["YxO2_gX_per_molO2"])
        S_min=float(params["S_min_g_L"]); X_min=float(params["X_min_g_L"])
        X_eff=max(X,X_min); S_eff=max(S,S_min)
        denomS=Kc*X_eff+S_eff
        fS=(S_eff/denomS) if denomS>0 else 0.0
        DO=float(getattr(env,"DO_mol_L",0.0) or 0.0)
        denomDO=K_DO+DO
        fDO=(DO/denomDO) if denomDO>0 else 0.0
        mu=mu_max*fS*fDO
        dX=mu*X_eff
        dS=-(mu/max(Yxs,1e-30))*X_eff
        if S<=S_min and dS<0: dS=0.0
        qO2=(mu/max(YxO2,1e-30))
        return np.array([dX,dS],dtype=float), {"mu_h":mu,"qO2_mol_gX_h":qO2}
