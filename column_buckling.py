import numpy as np
from scipy.optimize import bisect


def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """

    # פונקציית העזר למציאת השורש
    def f(P):

        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        sigma_max = (P / A) * (
            1 + (e * c / (r ** 2)) * (1 / np.cos(theta))
        )

        return sigma_max - sigma_allow

    # עומס אוילר הקריטי
    euler_load = (np.pi**2 * E * (A * r**2)) / (L**2)

    # גבולות לשיטת החצייה
    P_min = 1e-6
    P_max = 0.999 * euler_load

    # פתרון נומרי
    critical_load = bisect(f, P_min, P_max, xtol=1e-6)

    return critical_load
