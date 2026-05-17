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

    # פונקציית המאמץ לפי נוסחת הסקנט
    def stress_function(P):

        # חישוב הארגומנט של sec
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        # sec(x) = 1 / cos(x)
        sec_theta = 1 / np.cos(theta)

        # נוסחת הסקנט
        sigma_max = (P / A) * (
            1 + ((e * c) / (r ** 2)) * sec_theta
        )

        # פונקציית השורש
        return sigma_max - sigma_allow

    # גבולות התחלה לשיטת החצייה
    P_min = 1.0
    P_max = E * A * 0.99

    # מציאת השורש
    critical_load = bisect(
        stress_function,
        P_min,
        P_max,
        xtol=1e-6
    )

    return critical_load


if __name__ == "__main__":

    # דוגמת בדיקה
    L = 3000       # mm
    E = 200000     # MPa
    A = 5000       # mm^2
    r = 80         # mm
    c = 100        # mm
    e = 20         # mm
    sigma_allow = 250   # MPa

    P = find_critical_load(L, E, A, r, c, e, sigma_allow)

    print(f"Critical Load = {P:.3f} N")
